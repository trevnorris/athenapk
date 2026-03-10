#!/usr/bin/env python3
"""Compare exact discrete A_w gradient-energy increments against surviving ledger targets."""

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
    return [v for t, v in zip(times, values) if t >= t_lo and (t <= t_hi if include_right else t < t_hi)]


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


def add_rows(rows, case_dir, window, target_name, target_window, components):
    target_max, target_rms, target_final = finite_stats(target_window)
    rows.append({
        "case_dir": str(case_dir),
        "case_name": case_dir.name,
        "window": window,
        "target": target_name,
        "component": target_name,
        "max_abs_rate": target_max,
        "rms_abs_rate": target_rms,
        "final_rate": target_final,
        "share_of_target_rms": 1.0,
        "corr_with_target": 1.0,
        "fit_coeff": 1.0,
        "corrected_rms_abs_rate": 0.0,
        "rms_reduction_factor": math.inf,
    })
    for name, vals in components.items():
        max_abs, rms_abs, final_rate = finite_stats(vals)
        alpha = fit_single(target_window, vals)
        if math.isfinite(alpha):
            corrected = [y - alpha * x for y, x in zip(target_window, vals)]
            _, corrected_rms, _ = finite_stats(corrected)
            reduction = target_rms / corrected_rms if corrected_rms > 0.0 else math.inf
        else:
            corrected_rms = math.nan
            reduction = math.nan
        rows.append({
            "case_dir": str(case_dir),
            "case_name": case_dir.name,
            "window": window,
            "target": target_name,
            "component": name,
            "max_abs_rate": max_abs,
            "rms_abs_rate": rms_abs,
            "final_rate": final_rate,
            "share_of_target_rms": rms_abs / target_rms if target_rms > 0.0 else math.nan,
            "corr_with_target": pearson(vals, target_window),
            "fit_coeff": alpha,
            "corrected_rms_abs_rate": corrected_rms,
            "rms_reduction_factor": reduction,
        })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True)
    parser.add_argument("--window-fractions", default="0,0.333333,0.666667,1")
    parser.add_argument("--focus-window", default="late", choices=["early", "mid", "late", "all"])
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
    windows = ["early", "mid", "late"]
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
            "m4d_mixed_cz_daw_dz2_mode1",
            "m4d_mixed_cz_daw_dz2_mode2",
            "m4d_int_aw_mode0_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode0_gradz_energy_delta_exact_discrete",
            "m4d_int_aw_mode0_gradz_energy_delta_quadratic_discrete",
            "m4d_int_aw_mode1_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode1_gradz_energy_delta_exact_discrete",
            "m4d_int_aw_mode1_gradz_energy_delta_quadratic_discrete",
            "m4d_int_aw_mode2_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode2_gradz_energy_delta_exact_discrete",
            "m4d_int_aw_mode2_gradz_energy_delta_quadratic_discrete",
            "m4d_int_aw_mode3_gradz_power_da_dt_discrete",
            "m4d_int_aw_mode3_gradz_energy_delta_exact_discrete",
            "m4d_int_aw_mode3_gradz_energy_delta_quadratic_discrete",
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
        start = t_r[0]
        end = t_r[-1]
        span = end - start
        if span <= 0.0:
            continue

        mode_components = {}
        for mode in (1, 2):
            mode_components[mode] = {
                f"aw_mode{mode}_gradz_power_da_dt_discrete": derivative_series(times, cols[f"m4d_int_aw_mode{mode}_gradz_power_da_dt_discrete"])[1],
                f"aw_mode{mode}_gradz_energy_delta_exact_discrete": derivative_series(times, cols[f"m4d_int_aw_mode{mode}_gradz_energy_delta_exact_discrete"])[1],
                f"aw_mode{mode}_gradz_energy_delta_quadratic_discrete": derivative_series(times, cols[f"m4d_int_aw_mode{mode}_gradz_energy_delta_quadratic_discrete"])[1],
            }

        bulk_components = {
            "aw_modes0123_gradz_power_da_dt_discrete": sum_series(
                derivative_series(times, cols["m4d_int_aw_mode0_gradz_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode1_gradz_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode2_gradz_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode3_gradz_power_da_dt_discrete"])[1],
            ),
            "aw_modes0123_gradz_energy_delta_exact_discrete": sum_series(
                derivative_series(times, cols["m4d_int_aw_mode0_gradz_energy_delta_exact_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode1_gradz_energy_delta_exact_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode2_gradz_energy_delta_exact_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode3_gradz_energy_delta_exact_discrete"])[1],
            ),
            "aw_modes0123_gradz_energy_delta_quadratic_discrete": sum_series(
                derivative_series(times, cols["m4d_int_aw_mode0_gradz_energy_delta_quadratic_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode1_gradz_energy_delta_quadratic_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode2_gradz_energy_delta_quadratic_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode3_gradz_energy_delta_quadratic_discrete"])[1],
            ),
        }
        if any(len(v) != len(residual) for v in bulk_components.values()):
            continue

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
            bulk_window_components = {
                name: select_window(t_r, vals, t_lo, t_hi, include_right)
                for name, vals in bulk_components.items()
            }
            add_rows(rows, case_dir, window, "em_bulk_ledger_residual", residual_window, bulk_window_components)

            for mode in (1, 2):
                _, target = derivative_series(times, cols[f"m4d_mixed_cz_daw_dz2_mode{mode}"])
                target_window = select_window(t_r, target, t_lo, t_hi, include_right)
                if not target_window:
                    continue
                comps = {
                    name: select_window(t_r, vals, t_lo, t_hi, include_right)
                    for name, vals in mode_components[mode].items()
                }
                add_rows(rows, case_dir, window, f"d_mixed_cz_daw_dz2_mode{mode}_dt", target_window, comps)

    if not rows:
        raise SystemExit("No Step-177 rows were computed.")

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
            grouped.setdefault((row["case_name"], row["target"], row["window"]), []).append(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-177 A_w Exact Delta Probe\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for key in sorted(grouped.keys(), key=lambda x: (x[0], x[1], window_order.get(x[2], 99))):
            case_name, target_name, window = key
            grp = grouped[key]
            target = next((r for r in grp if r["component"] == target_name), None)
            f.write(f"## {case_name} target={target_name} window={window}\n\n")
            if target:
                f.write(f"- target_rms_abs_rate: `{fmt(target['rms_abs_rate'])}`\n\n")
            ranked = [r for r in grp if r["component"] != target_name]
            ranked.sort(
                key=lambda r: (
                    r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf,
                    abs(r["corr_with_target"]) if math.isfinite(r["corr_with_target"]) else -1.0,
                ),
                reverse=True,
            )
            f.write("| component | rms_abs_rate | share_of_target_rms | corr_with_target | fit_coeff | corrected_rms_abs_rate | rms_reduction_factor | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in ranked:
                f.write(
                    f"| {row['component']} | {fmt(row['rms_abs_rate'])} | {fmt(row['share_of_target_rms'])} | {fmt(row['corr_with_target'])} | {fmt(row['fit_coeff'])} | {fmt(row['corrected_rms_abs_rate'])} | {fmt(row['rms_reduction_factor'])} | {fmt(row['max_abs_rate'])} | {fmt(row['final_rate'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
