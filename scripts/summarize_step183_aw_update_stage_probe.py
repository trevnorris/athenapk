#!/usr/bin/env python3
"""Compare bulk residual against discrete A_w update-stage sums."""

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
    return out_t, out_v


def select_window(times, values, t_lo, t_hi, include_right=False):
    return [
        v for t, v in zip(times, values)
        if t >= t_lo and (t <= t_hi if include_right else t < t_hi)
    ]


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return math.nan, math.nan, math.nan
    return max(abs(v) for v in finite), math.sqrt(sum(v * v for v in finite) / len(finite)), finite[-1]


def pearson(xs, ys):
    if len(xs) != len(ys) or len(xs) < 2:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def fit_single(y, x):
    xx = sum(v * v for v in x)
    if xx <= 0.0:
        return math.nan
    return sum(a * b for a, b in zip(y, x)) / xx


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    if v == math.inf:
        return "inf"
    return f"{float(v):.6e}"


def sum_series(*series_list):
    return [sum(vals) for vals in zip(*series_list)]


def sub_series(a, b):
    return [x - y for x, y in zip(a, b)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]
    rows = []

    for case_dir in case_dirs:
        hst = case_dir / "harris_full.out1.hst"
        if not hst.exists():
            continue
        cols = hsm.parse_hst(hst)
        required = [
            "time", "m4d_em_u_bulk", "m4d_int_ja_ea", "m4d_int_jw_ew",
            "m4d_int_aw_gradx_power_pi_drive_discrete_sum",
            "m4d_int_aw_gradz_power_pi_drive_discrete_sum",
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
        start, end = t_r[0], t_r[-1]
        t_lo = start + (end - start) * 2.0 / 3.0
        y = select_window(t_r, residual, t_lo, end, include_right=True)
        _, res_rms, _ = finite_stats(y)

        pi_drive_x = select_window(
            t_r, derivative_series(times, cols["m4d_int_aw_gradx_power_pi_drive_discrete_sum"])[1],
            t_lo, end, include_right=True)
        pi_drive_z = select_window(
            t_r, derivative_series(times, cols["m4d_int_aw_gradz_power_pi_drive_discrete_sum"])[1],
            t_lo, end, include_right=True)
        pi_drive_xz = sum_series(pi_drive_x, pi_drive_z)
        full_da_dt_xz = select_window(
            t_r,
            sum_series(
                derivative_series(times, cols["m4d_int_aw_mode0_gradx_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode1_gradx_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode2_gradx_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode3_gradx_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode0_gradz_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode1_gradz_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode2_gradz_power_da_dt_discrete"])[1],
                derivative_series(times, cols["m4d_int_aw_mode3_gradz_power_da_dt_discrete"])[1],
            ),
            t_lo, end, include_right=True,
        )
        mix_remainder_xz = sub_series(full_da_dt_xz, pi_drive_xz)

        components = {
            "pi_drive_x_sum": pi_drive_x,
            "pi_drive_z_sum": pi_drive_z,
            "pi_drive_xz_sum": pi_drive_xz,
            "mix_remainder_xz_sum": mix_remainder_xz,
            "full_da_dt_xz_sum": full_da_dt_xz,
        }

        rows.append({
            "case_name": case_dir.name,
            "component": "em_bulk_ledger_residual",
            "rms_abs_rate": res_rms,
            "corr_with_residual": 1.0,
            "fit_coeff": 1.0,
            "corrected_rms_abs_rate": 0.0,
            "rms_reduction_factor": math.inf,
        })
        for name, vals in components.items():
            _, rms_abs, _ = finite_stats(vals)
            alpha = fit_single(y, vals)
            corrected = [yy - alpha * vv for yy, vv in zip(y, vals)] if math.isfinite(alpha) else []
            _, corrected_rms, _ = finite_stats(corrected)
            reduction = res_rms / corrected_rms if corrected_rms > 0 else math.inf
            rows.append({
                "case_name": case_dir.name,
                "component": name,
                "rms_abs_rate": rms_abs,
                "corr_with_residual": pearson(vals, y),
                "fit_coeff": alpha,
                "corrected_rms_abs_rate": corrected_rms,
                "rms_reduction_factor": reduction,
            })

    if not rows:
        raise SystemExit("No Step-183 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-183 A_w Update Stage Probe\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            f.write(f"## {case_name} window=late\n\n")
            target = next(r for r in rows if r["case_name"] == case_name and r["component"] == "em_bulk_ledger_residual")
            f.write(f"- residual_rms_abs_rate: `{fmt(target['rms_abs_rate'])}`\n\n")
            f.write("| component | rms_abs_rate | corr_with_residual | fit_coeff | corrected_rms_abs_rate | rms_reduction_factor |\n")
            f.write("|:---|---:|---:|---:|---:|---:|\n")
            ranked = [r for r in rows if r["case_name"] == case_name and r["component"] != "em_bulk_ledger_residual"]
            ranked.sort(key=lambda r: r["rms_reduction_factor"], reverse=True)
            for r in ranked:
                f.write(
                    f"| {r['component']} | {fmt(r['rms_abs_rate'])} | {fmt(r['corr_with_residual'])} | "
                    f"{fmt(r['fit_coeff'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
