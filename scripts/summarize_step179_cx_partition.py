#!/usr/bin/env python3
"""Summarize Cx partition of the mixed-C drift."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_harris_scan_matrix_module(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_window_fractions(raw: str):
    vals = [float(x.strip()) for x in raw.split(",") if x.strip()]
    if len(vals) < 2:
        raise ValueError("Need at least two fractions")
    if abs(vals[0]) > 1.0e-12 or abs(vals[-1] - 1.0) > 1.0e-12:
        raise ValueError("Fractions must start at 0 and end at 1")
    for i in range(1, len(vals)):
        if vals[i] <= vals[i - 1]:
            raise ValueError("Fractions must be strictly increasing")
    return vals


def derivative_series(times, values):
    if values is None or len(times) < 2 or len(values) != len(times):
        return ([], [])
    out_t = []
    out_v = []
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        out_t.append(times[i])
        out_v.append((values[i] - values[i - 1]) / dt)
    return (out_t, out_v)


def pearson(xs, ys):
    if len(xs) != len(ys) or len(xs) < 2:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) * (x - mx) for x in xs)
    vy = sum((y - my) * (y - my) for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    max_abs = max(abs(v) for v in finite)
    rms_abs = math.sqrt(sum(v * v for v in finite) / len(finite))
    final = finite[-1]
    return (max_abs, rms_abs, final)


def select_window(times, values, t_lo, t_hi, include_right=False):
    return [v for t, v in zip(times, values) if t >= t_lo and (t <= t_hi if include_right else t < t_hi)]


def fmt(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def compute_rows(case_dir, cols, fractions):
    times = cols["time"]
    t_r, cx_rate = derivative_series(times, cols["m4d_mixed_cx2"])
    _, bulk_rate = derivative_series(times, cols["m4d_em_u_bulk"])
    if not t_r or len(cx_rate) != len(bulk_rate):
        return []

    component_rates = {
        "d_mixed_cx_daw_dx2_dt": derivative_series(times, cols["m4d_mixed_cx_daw_dx2"])[1],
        "d_mixed_cx_dw_ax2_dt": derivative_series(times, cols["m4d_mixed_cx_dw_ax2"])[1],
        "d_mixed_cx_cross_dt": derivative_series(times, cols["m4d_mixed_cx_cross"])[1],
    }
    if any(len(v) != len(cx_rate) for v in component_rates.values()):
        return []

    start = t_r[0]
    end = t_r[-1]
    span = end - start
    if span <= 0.0:
        return []

    windows = ["early", "mid", "late"]
    rows = []
    for idx in range(len(fractions) - 1):
        f_lo = fractions[idx]
        f_hi = fractions[idx + 1]
        window = windows[idx] if idx < len(windows) else f"w{idx}"
        t_lo = start + f_lo * span
        t_hi = start + f_hi * span
        include_right = idx == len(fractions) - 2

        bulk_window = select_window(t_r, bulk_rate, t_lo, t_hi, include_right)
        cx_window = select_window(t_r, cx_rate, t_lo, t_hi, include_right)
        bulk_max, bulk_rms, _ = finite_stats(bulk_window)
        cx_max, cx_rms, cx_final = finite_stats(cx_window)
        if not cx_window:
            continue

        rows.append({
            "case_dir": str(case_dir),
            "case_name": case_dir.name,
            "window": window,
            "component": "d_mixed_cx2_dt",
            "max_abs_rate": cx_max,
            "rms_abs_rate": cx_rms,
            "final_rate": cx_final,
            "share_of_cx_rms": 1.0,
            "share_of_bulk_rms": cx_rms / bulk_rms if bulk_rms > 0.0 else math.nan,
            "corr_with_cx": 1.0,
            "corr_with_bulk": pearson(cx_window, bulk_window),
            "cx_rms_abs_rate": cx_rms,
            "bulk_rms_abs_rate": bulk_rms,
        })

        for name, series in component_rates.items():
            window_vals = select_window(t_r, series, t_lo, t_hi, include_right)
            max_abs, rms_abs, final_rate = finite_stats(window_vals)
            rows.append({
                "case_dir": str(case_dir),
                "case_name": case_dir.name,
                "window": window,
                "component": name,
                "max_abs_rate": max_abs,
                "rms_abs_rate": rms_abs,
                "final_rate": final_rate,
                "share_of_cx_rms": rms_abs / cx_rms if cx_rms > 0.0 else math.nan,
                "share_of_bulk_rms": rms_abs / bulk_rms if bulk_rms > 0.0 else math.nan,
                "corr_with_cx": pearson(window_vals, cx_window),
                "corr_with_bulk": pearson(window_vals, bulk_window),
                "cx_rms_abs_rate": cx_rms,
                "bulk_rms_abs_rate": bulk_rms,
            })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True)
    parser.add_argument("--window-fractions", default="0,0.333333,0.666667,1")
    parser.add_argument("--focus-window", default="late", choices=["early", "mid", "late", "all"])
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)
    dirs = [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]
    if not dirs:
        raise SystemExit(f"No directories matched: {args.glob}")

    rows = []
    for d in dirs:
        full_hst = d / "harris_full.out1.hst"
        if not full_hst.exists():
            continue
        cols = hsm.parse_hst(full_hst)
        required = [
            "time",
            "m4d_em_u_bulk",
            "m4d_mixed_cx2",
            "m4d_mixed_cx_daw_dx2",
            "m4d_mixed_cx_dw_ax2",
            "m4d_mixed_cx_cross",
        ]
        if not all(col in cols for col in required):
            continue
        rows.extend(compute_rows(d, cols, fractions))

    if not rows:
        raise SystemExit("No Cx partition rows were computed.")

    window_order = {"early": 0, "mid": 1, "late": 2}
    rows.sort(key=lambda r: (r["case_name"], window_order.get(r["window"], 99), r["component"]))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    focus_windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
    grouped = {}
    for row in rows:
        if row["window"] in focus_windows:
            grouped.setdefault((row["case_name"], row["window"]), []).append(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-179 Cx Partition Audit\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for (case_name, window) in sorted(grouped.keys(), key=lambda x: (x[0], window_order.get(x[1], 99))):
            grp = grouped[(case_name, window)]
            cx = next((r for r in grp if r["component"] == "d_mixed_cx2_dt"), None)
            f.write(f"## {case_name} window={window}\n\n")
            if cx is not None:
                f.write(f"- mixed_cx_rms_abs_rate: `{fmt(cx['rms_abs_rate'])}`\n")
                f.write(f"- bulk_rms_abs_rate: `{fmt(cx['bulk_rms_abs_rate'])}`\n\n")
            ranked = [r for r in grp if r["component"] != "d_mixed_cx2_dt"]
            ranked.sort(key=lambda r: abs(r["share_of_cx_rms"]) if math.isfinite(r["share_of_cx_rms"]) else -1.0, reverse=True)
            f.write("| component | rms_abs_rate | share_of_cx_rms | share_of_bulk_rms | corr_with_cx | corr_with_bulk | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in ranked:
                f.write(
                    f"| {row['component']} | {fmt(row['rms_abs_rate'])} | {fmt(row['share_of_cx_rms'])} | {fmt(row['share_of_bulk_rms'])} | {fmt(row['corr_with_cx'])} | {fmt(row['corr_with_bulk'])} | {fmt(row['max_abs_rate'])} | {fmt(row['final_rate'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
