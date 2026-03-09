#!/usr/bin/env python3
"""Summarize per-mode A_w bookkeeping against diagonal d_z a_w energy drift."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_hsm(script_dir: Path):
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


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    return f"{float(v):.6e}"


def compute_rows(case_dir, cols, fractions):
    times = cols["time"]
    windows = ["early", "mid", "late"]
    rows = []

    for mode in (1, 2):
        target_key = f"m4d_mixed_cz_daw_dz2_mode{mode}"
        t_r, target = derivative_series(times, cols[target_key])
        if not t_r:
            continue
        start = t_r[0]
        end = t_r[-1]
        span = end - start
        if span <= 0.0:
            continue

        series_keys = {
            f"aw_mode{mode}_pi_drive": f"m4d_int_aw_mode{mode}_pi_drive",
            f"aw_mode{mode}_rhs": f"m4d_int_aw_mode{mode}_rhs",
            f"aw_mode{mode}_mix_pi": f"m4d_int_aw_mode{mode}_mix_pi",
            f"aw_mode{mode}_mix_a": f"m4d_int_aw_mode{mode}_mix_a",
            f"aw_mode{mode}_da_dt": f"m4d_int_aw_mode{mode}_da_dt",
        }
        components = {name: derivative_series(times, cols[key])[1] for name, key in series_keys.items()}
        if any(len(v) != len(target) for v in components.values()):
            continue

        for idx in range(len(fractions) - 1):
            f_lo = fractions[idx]
            f_hi = fractions[idx + 1]
            window = windows[idx] if idx < len(windows) else f"w{idx}"
            t_lo = start + f_lo * span
            t_hi = start + f_hi * span
            include_right = idx == len(fractions) - 2

            target_window = select_window(t_r, target, t_lo, t_hi, include_right)
            target_max, target_rms, target_final = finite_stats(target_window)
            if not target_window:
                continue

            rows.append({
                "case_dir": str(case_dir),
                "case_name": case_dir.name,
                "mode": mode,
                "window": window,
                "component": f"d_mixed_cz_daw_dz2_mode{mode}_dt",
                "max_abs_rate": target_max,
                "rms_abs_rate": target_rms,
                "final_rate": target_final,
                "share_of_target_rms": 1.0,
                "corr_with_target": 1.0,
            })

            for name, series in components.items():
                window_vals = select_window(t_r, series, t_lo, t_hi, include_right)
                max_abs, rms_abs, final_rate = finite_stats(window_vals)
                rows.append({
                    "case_dir": str(case_dir),
                    "case_name": case_dir.name,
                    "mode": mode,
                    "window": window,
                    "component": name,
                    "max_abs_rate": max_abs,
                    "rms_abs_rate": rms_abs,
                    "final_rate": final_rate,
                    "share_of_target_rms": rms_abs / target_rms if target_rms > 0.0 else math.nan,
                    "corr_with_target": pearson(window_vals, target_window),
                })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True)
    parser.add_argument("--window-fractions", default="0,0.333333,0.666667,1")
    parser.add_argument("--focus-window", default="late", choices=["early", "mid", "late", "all"])
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]

    rows = []
    for case_dir in case_dirs:
        hst = case_dir / "harris_full.out1.hst"
        if not hst.exists():
            continue
        cols = hsm.parse_hst(hst)
        required = [
            "time",
            "m4d_mixed_cz_daw_dz2_mode1",
            "m4d_mixed_cz_daw_dz2_mode2",
            "m4d_int_aw_mode1_pi_drive",
            "m4d_int_aw_mode1_rhs",
            "m4d_int_aw_mode1_mix_pi",
            "m4d_int_aw_mode1_mix_a",
            "m4d_int_aw_mode1_da_dt",
            "m4d_int_aw_mode2_pi_drive",
            "m4d_int_aw_mode2_rhs",
            "m4d_int_aw_mode2_mix_pi",
            "m4d_int_aw_mode2_mix_a",
            "m4d_int_aw_mode2_da_dt",
        ]
        if not all(k in cols for k in required):
            continue
        rows.extend(compute_rows(case_dir, cols, fractions))

    if not rows:
        raise SystemExit("No A_w bookkeeping rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    focus_windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
    window_order = {"early": 0, "mid": 1, "late": 2}
    grouped = {}
    for row in rows:
        if row["window"] in focus_windows:
            grouped.setdefault((row["case_name"], row["mode"], row["window"]), []).append(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-164 A_w Mode Bookkeeping Audit\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1], window_order.get(x[2], 99))):
            case_name, mode, window = key
            grp = grouped[key]
            target = next((r for r in grp if r["component"] == f"d_mixed_cz_daw_dz2_mode{mode}_dt"), None)
            f.write(f"## {case_name} mode={mode} window={window}\n\n")
            if target:
                f.write(f"- target_rms_abs_rate: `{fmt(target['rms_abs_rate'])}`\n\n")
            ranked = [r for r in grp if not r["component"].startswith("d_mixed_cz_daw_dz2_mode")]
            ranked.sort(key=lambda r: abs(r["share_of_target_rms"]) if math.isfinite(r["share_of_target_rms"]) else -1.0, reverse=True)
            f.write("| component | rms_abs_rate | share_of_target_rms | corr_with_target | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|\n")
            for row in ranked:
                f.write(
                    f"| {row['component']} | {fmt(row['rms_abs_rate'])} | {fmt(row['share_of_target_rms'])} | {fmt(row['corr_with_target'])} | {fmt(row['max_abs_rate'])} | {fmt(row['final_rate'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
