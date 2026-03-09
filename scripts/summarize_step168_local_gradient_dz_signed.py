#!/usr/bin/env python3
"""Summarize Step-168 local-gradient dz signed diagnostics."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_hsm(script_dir: Path):
    path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("hsm", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def parse_window_fractions(text: str):
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    if len(vals) < 2 or vals[0] != 0.0 or vals[-1] != 1.0:
        raise ValueError("window fractions must start at 0 and end at 1")
    return vals


def derivative_series(times, values):
    out_t, out_v = [], []
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        out_t.append(times[i])
        out_v.append((values[i] - values[i - 1]) / dt)
    return out_t, out_v


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
        full_target_key = f"m4d_int_aw_mode{mode}_gradz_power_mix_w0_local_gradient"
        dz_target_key = f"m4d_int_aw_mode{mode}_gradz_power_mix_w0_local_dz"
        t_r, full_target = derivative_series(times, cols[full_target_key])
        t_dz, dz_target = derivative_series(times, cols[dz_target_key])
        if not t_r or t_r != t_dz:
            continue
        start = t_r[0]
        end = t_r[-1]
        span = end - start
        if span <= 0.0:
            continue
        series_keys = {
            f"aw_mode{mode}_gradz_power_mix_w0_local_dz_signed": f"m4d_int_aw_mode{mode}_gradz_power_mix_w0_local_dz_signed",
            f"aw_mode{mode}_gradz_power_mix_w0_local_dz_signed_mode1": f"m4d_int_aw_mode{mode}_gradz_power_mix_w0_local_dz_signed_mode1",
            f"aw_mode{mode}_gradz_power_mix_w0_local_dz_signed_mode2": f"m4d_int_aw_mode{mode}_gradz_power_mix_w0_local_dz_signed_mode2",
        }
        components = {name: derivative_series(times, cols[key])[1] for name, key in series_keys.items()}
        if any(len(v) != len(full_target) for v in components.values()):
            continue
        for idx in range(len(fractions) - 1):
            f_lo = fractions[idx]
            f_hi = fractions[idx + 1]
            window = windows[idx] if idx < len(windows) else f"w{idx}"
            t_lo = start + f_lo * span
            t_hi = start + f_hi * span
            include_right = idx == len(fractions) - 2
            full_target_window = select_window(t_r, full_target, t_lo, t_hi, include_right)
            dz_target_window = select_window(t_r, dz_target, t_lo, t_hi, include_right)
            if not full_target_window or not dz_target_window:
                continue
            full_max, full_rms, full_final = finite_stats(full_target_window)
            dz_max, dz_rms, dz_final = finite_stats(dz_target_window)
            rows.append({
                "case_dir": str(case_dir),
                "case_name": case_dir.name,
                "mode": mode,
                "window": window,
                "component": f"d_aw_mode{mode}_gradz_power_mix_w0_local_gradient_dt",
                "max_abs_rate": full_max,
                "rms_abs_rate": full_rms,
                "final_rate": full_final,
                "share_of_dz_target_rms": full_rms / dz_rms if dz_rms > 0.0 else math.nan,
                "share_of_full_target_rms": 1.0,
                "corr_with_dz_target": pearson(full_target_window, dz_target_window),
                "corr_with_full_target": 1.0,
            })
            rows.append({
                "case_dir": str(case_dir),
                "case_name": case_dir.name,
                "mode": mode,
                "window": window,
                "component": f"d_aw_mode{mode}_gradz_power_mix_w0_local_dz_dt",
                "max_abs_rate": dz_max,
                "rms_abs_rate": dz_rms,
                "final_rate": dz_final,
                "share_of_dz_target_rms": 1.0,
                "share_of_full_target_rms": dz_rms / full_rms if full_rms > 0.0 else math.nan,
                "corr_with_dz_target": 1.0,
                "corr_with_full_target": pearson(dz_target_window, full_target_window),
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
                    "share_of_dz_target_rms": rms_abs / dz_rms if dz_rms > 0.0 else math.nan,
                    "share_of_full_target_rms": rms_abs / full_rms if full_rms > 0.0 else math.nan,
                    "corr_with_dz_target": pearson(window_vals, dz_target_window),
                    "corr_with_full_target": pearson(window_vals, full_target_window),
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
            "m4d_int_aw_mode1_gradz_power_mix_w0_local_gradient",
            "m4d_int_aw_mode1_gradz_power_mix_w0_local_dz",
            "m4d_int_aw_mode1_gradz_power_mix_w0_local_dz_signed",
            "m4d_int_aw_mode1_gradz_power_mix_w0_local_dz_signed_mode1",
            "m4d_int_aw_mode1_gradz_power_mix_w0_local_dz_signed_mode2",
            "m4d_int_aw_mode2_gradz_power_mix_w0_local_gradient",
            "m4d_int_aw_mode2_gradz_power_mix_w0_local_dz",
            "m4d_int_aw_mode2_gradz_power_mix_w0_local_dz_signed",
            "m4d_int_aw_mode2_gradz_power_mix_w0_local_dz_signed_mode1",
            "m4d_int_aw_mode2_gradz_power_mix_w0_local_dz_signed_mode2",
        ]
        if not all(k in cols for k in required):
            continue
        rows.extend(compute_rows(case_dir, cols, fractions))

    if not rows:
        raise SystemExit("No local-gradient dz signed rows were computed.")

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
        f.write("# Step-168 Local-Gradient dz Signed Audit\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1], window_order.get(x[2], 99))):
            case_name, mode, window = key
            grp = grouped[key]
            full_target = next((r for r in grp if r["component"] == f"d_aw_mode{mode}_gradz_power_mix_w0_local_gradient_dt"), None)
            dz_target = next((r for r in grp if r["component"] == f"d_aw_mode{mode}_gradz_power_mix_w0_local_dz_dt"), None)
            f.write(f"## {case_name} mode={mode} window={window}\n\n")
            if full_target:
                f.write(f"- full_target_rms_abs_rate: `{fmt(full_target['rms_abs_rate'])}`\n")
            if dz_target:
                f.write(f"- dz_target_rms_abs_rate: `{fmt(dz_target['rms_abs_rate'])}`\n")
            f.write("\n")
            ranked = [r for r in grp if r["component"] not in {
                f"d_aw_mode{mode}_gradz_power_mix_w0_local_gradient_dt",
                f"d_aw_mode{mode}_gradz_power_mix_w0_local_dz_dt",
            }]
            ranked.sort(key=lambda r: abs(r["share_of_dz_target_rms"]) if math.isfinite(r["share_of_dz_target_rms"]) else -1.0, reverse=True)
            f.write("| component | rms_abs_rate | share_of_dz_target_rms | share_of_full_target_rms | corr_with_dz_target | corr_with_full_target | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in ranked:
                f.write(
                    f"| {row['component']} | {fmt(row['rms_abs_rate'])} | {fmt(row['share_of_dz_target_rms'])} | {fmt(row['share_of_full_target_rms'])} | {fmt(row['corr_with_dz_target'])} | {fmt(row['corr_with_full_target'])} | {fmt(row['max_abs_rate'])} | {fmt(row['final_rate'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
