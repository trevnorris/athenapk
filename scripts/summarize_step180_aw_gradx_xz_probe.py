#!/usr/bin/env python3
"""Compare discrete A_w grad-x update-power proxies and combined x+z sums."""

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


def select_window(times, values, t_lo, t_hi, include_right=False):
    return [
        v
        for t, v in zip(times, values)
        if t >= t_lo and (t <= t_hi if include_right else t < t_hi)
    ]


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    return (
        max(abs(v) for v in finite),
        math.sqrt(sum(v * v for v in finite) / len(finite)),
        finite[-1],
    )


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


def dot(xs, ys):
    return sum(x * y for x, y in zip(xs, ys))


def fit_single(y, x):
    xx = dot(x, x)
    if xx <= 0.0:
        return math.nan
    return dot(y, x) / xx


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    return f"{float(v):.6e}"


def sum_series(*series_list):
    return [sum(vals) for vals in zip(*series_list)]


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
            "m4d_em_u_bulk",
            "m4d_int_ja_ea",
            "m4d_int_jw_ew",
            "m4d_mixed_cx_daw_dx2_mode0",
            "m4d_mixed_cx_daw_dx2_mode1",
            "m4d_mixed_cx_daw_dx2_mode2",
            "m4d_mixed_cx_daw_dx2_mode3",
            "m4d_mixed_cx_daw_dx2_mode_diag",
            "m4d_mixed_cz_daw_dz2_mode_diag",
            "m4d_int_aw_mode0_gradx_power_da_dt_discrete",
            "m4d_int_aw_mode1_gradx_power_da_dt_discrete",
            "m4d_int_aw_mode2_gradx_power_da_dt_discrete",
            "m4d_int_aw_mode3_gradx_power_da_dt_discrete",
            "m4d_int_aw_mode0_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode1_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode2_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode3_gradz_power_da_dt_discrete",
        ]
        if not all(k in cols for k in required):
            continue

        times = cols["time"]
        t_r, du = derivative_series(times, cols["m4d_em_u_bulk"])
        _, dja = derivative_series(times, cols["m4d_int_ja_ea"])
        _, djw = derivative_series(times, cols["m4d_int_jw_ew"])
        residual = [u + a + w for u, a, w in zip(du, dja, djw)]
        if not t_r:
            continue

        gradx_components = {
            f"aw_mode{mode}_gradx_power_da_dt_discrete": derivative_series(
                times, cols[f"m4d_int_aw_mode{mode}_gradx_power_da_dt_discrete"]
            )[1]
            for mode in range(4)
        }
        gradz_components = {
            f"aw_mode{mode}_gradz_power_da_dt_discrete": derivative_series(
                times, cols[f"m4d_int_aw_mode{mode}_gradz_power_da_dt_discrete"]
            )[1]
            for mode in range(4)
        }
        target_x = {
            f"d_mixed_cx_daw_dx2_mode{mode}_dt": derivative_series(
                times, cols[f"m4d_mixed_cx_daw_dx2_mode{mode}"]
            )[1]
            for mode in range(4)
        }
        target_x_diag = derivative_series(times, cols["m4d_mixed_cx_daw_dx2_mode_diag"])[1]
        target_z_diag = derivative_series(times, cols["m4d_mixed_cz_daw_dz2_mode_diag"])[1]
        if any(len(v) != len(residual) for v in list(gradx_components.values()) + list(gradz_components.values()) + list(target_x.values()) + [target_x_diag, target_z_diag]):
            continue

        aggregate = {
            "aw_modes0123_gradx_power_da_dt_discrete": sum_series(*[gradx_components[f"aw_mode{mode}_gradx_power_da_dt_discrete"] for mode in range(4)]),
            "aw_modes0123_gradz_power_da_dt_discrete": sum_series(*[gradz_components[f"aw_mode{mode}_gradz_power_da_dt_discrete"] for mode in range(4)]),
        }
        aggregate["aw_modes0123_gradxz_power_da_dt_discrete"] = sum_series(
            aggregate["aw_modes0123_gradx_power_da_dt_discrete"],
            aggregate["aw_modes0123_gradz_power_da_dt_discrete"],
        )
        aggregate["d_mixed_cx_daw_dx2_mode_diag"] = target_x_diag
        aggregate["d_mixed_cz_daw_dz2_mode_diag"] = target_z_diag

        start = t_r[0]
        end = t_r[-1]
        span = end - start
        if span <= 0.0:
            continue
        windows = ["early", "mid", "late"]

        for idx in range(len(fractions) - 1):
            window = windows[idx] if idx < len(windows) else f"w{idx}"
            if args.focus_window != "all" and window != args.focus_window:
                continue
            t_lo = start + fractions[idx] * span
            t_hi = start + fractions[idx + 1] * span
            include_right = idx == len(fractions) - 2

            residual_window = select_window(t_r, residual, t_lo, t_hi, include_right)
            if not residual_window:
                continue
            res_max, res_rms, res_final = finite_stats(residual_window)
            rows.append({
                "case_name": case_dir.name,
                "window": window,
                "metric_group": "bulk",
                "component": "em_bulk_ledger_residual",
                "max_abs_rate": res_max,
                "rms_abs_rate": res_rms,
                "final_rate": res_final,
                "share_of_target_rms": 1.0,
                "corr_with_target": 1.0,
                "fit_coeff": 1.0,
                "corrected_rms_abs_rate": 0.0,
                "rms_reduction_factor": math.inf,
            })
            for name, series in aggregate.items():
                vals = select_window(t_r, series, t_lo, t_hi, include_right)
                max_abs, rms_abs, final_rate = finite_stats(vals)
                alpha = fit_single(residual_window, vals)
                corrected = [y - alpha * x for y, x in zip(residual_window, vals)] if math.isfinite(alpha) else []
                _, corrected_rms, _ = finite_stats(corrected)
                reduction = res_rms / corrected_rms if corrected_rms > 0.0 else math.inf
                rows.append({
                    "case_name": case_dir.name,
                    "window": window,
                    "metric_group": "bulk",
                    "component": name,
                    "max_abs_rate": max_abs,
                    "rms_abs_rate": rms_abs,
                    "final_rate": final_rate,
                    "share_of_target_rms": rms_abs / res_rms if res_rms > 0.0 else math.nan,
                    "corr_with_target": pearson(vals, residual_window),
                    "fit_coeff": alpha,
                    "corrected_rms_abs_rate": corrected_rms,
                    "rms_reduction_factor": reduction,
                })

            for mode in range(4):
                target_vals = select_window(t_r, target_x[f"d_mixed_cx_daw_dx2_mode{mode}_dt"], t_lo, t_hi, include_right)
                if not target_vals:
                    continue
                target_max, target_rms, target_final = finite_stats(target_vals)
                rows.append({
                    "case_name": case_dir.name,
                    "window": window,
                    "metric_group": f"mode{mode}_x",
                    "component": f"d_mixed_cx_daw_dx2_mode{mode}_dt",
                    "max_abs_rate": target_max,
                    "rms_abs_rate": target_rms,
                    "final_rate": target_final,
                    "share_of_target_rms": 1.0,
                    "corr_with_target": 1.0,
                    "fit_coeff": 1.0,
                    "corrected_rms_abs_rate": 0.0,
                    "rms_reduction_factor": math.inf,
                })
                comp_name = f"aw_mode{mode}_gradx_power_da_dt_discrete"
                comp_vals = select_window(t_r, gradx_components[comp_name], t_lo, t_hi, include_right)
                max_abs, rms_abs, final_rate = finite_stats(comp_vals)
                alpha = fit_single(target_vals, comp_vals)
                corrected = [y - alpha * x for y, x in zip(target_vals, comp_vals)] if math.isfinite(alpha) else []
                _, corrected_rms, _ = finite_stats(corrected)
                reduction = target_rms / corrected_rms if corrected_rms > 0.0 else math.inf
                rows.append({
                    "case_name": case_dir.name,
                    "window": window,
                    "metric_group": f"mode{mode}_x",
                    "component": comp_name,
                    "max_abs_rate": max_abs,
                    "rms_abs_rate": rms_abs,
                    "final_rate": final_rate,
                    "share_of_target_rms": rms_abs / target_rms if target_rms > 0.0 else math.nan,
                    "corr_with_target": pearson(comp_vals, target_vals),
                    "fit_coeff": alpha,
                    "corrected_rms_abs_rate": corrected_rms,
                    "rms_reduction_factor": reduction,
                })

    if not rows:
        raise SystemExit("No Step-180 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-180 A_w Grad-X / X+Z Discrete Update-Power Probe\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for row in rows:
            pass
        grouped = {}
        order = {"early": 0, "mid": 1, "late": 2}
        for row in rows:
            if args.focus_window == "all" or row["window"] == args.focus_window:
                grouped.setdefault((row["case_name"], row["metric_group"], row["window"]), []).append(row)
        for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1], order.get(x[2], 99))):
            case_name, metric_group, window = key
            grp = grouped[key]
            f.write(f"## {case_name} group={metric_group} window={window}\n\n")
            target = next((r for r in grp if r["component"].startswith("d_mixed") or r["component"] == "em_bulk_ledger_residual"), None)
            if target:
                f.write(f"- target_rms_abs_rate: `{fmt(target['rms_abs_rate'])}`\n\n")
            ranked = [r for r in grp if r is not target]
            ranked.sort(key=lambda r: (r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf, abs(r["corr_with_target"]) if math.isfinite(r["corr_with_target"]) else -1.0), reverse=True)
            f.write("| component | rms_abs_rate | share_of_target_rms | corr_with_target | fit_coeff | corrected_rms_abs_rate | rms_reduction_factor | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in ranked:
                f.write(
                    f"| {r['component']} | {fmt(r['rms_abs_rate'])} | {fmt(r['share_of_target_rms'])} | {fmt(r['corr_with_target'])} | {fmt(r['fit_coeff'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} | {fmt(r['max_abs_rate'])} | {fmt(r['final_rate'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
