#!/usr/bin/env python3
"""Validate fixed post-fix z-block correction rules across time windows."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


Z_COEFF_FIXED = 8.304215e-01
X_COEFF_FIXED = 8.288144e-02


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
    out_t, out_v = [], []
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        out_t.append(times[i])
        out_v.append((values[i] - values[i - 1]) / dt)
    return out_t, out_v


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


def sum_series(*series):
    return [sum(vals) for vals in zip(*series)]


def fit_single(y, x):
    xx = sum(v * v for v in x)
    if xx <= 0.0:
        return math.nan
    return sum(a * b for a, b in zip(y, x)) / xx


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


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    if v == math.inf:
        return "inf"
    return f"{float(v):.6e}"


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
            "time",
            "m4d_em_u_bulk",
            "m4d_int_ja_ea",
            "m4d_int_jw_ew",
            "m4d_int_cx_energy_delta_exact_discrete_sum",
            "m4d_int_cz_energy_delta_exact_discrete_sum",
        ] + [f"m4d_int_aw_mode{n}_gradx_power_da_dt_discrete" for n in range(4)] + [
            f"m4d_int_aw_mode{n}_gradz_power_da_dt_discrete" for n in range(4)
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
        span = end - start
        windows = {
            "early": (start, start + span / 3.0),
            "mid": (start + span / 3.0, start + 2.0 * span / 3.0),
            "late": (start + 2.0 * span / 3.0, end),
        }

        _, cx = derivative_series(times, cols["m4d_int_cx_energy_delta_exact_discrete_sum"])
        _, cz = derivative_series(times, cols["m4d_int_cz_energy_delta_exact_discrete_sum"])
        x_da = []
        z_da = []
        for n in range(4):
            _, dx = derivative_series(times, cols[f"m4d_int_aw_mode{n}_gradx_power_da_dt_discrete"])
            _, dz = derivative_series(times, cols[f"m4d_int_aw_mode{n}_gradz_power_da_dt_discrete"])
            x_da.append(dx)
            z_da.append(dz)
        x_block = sum_series(cx, *x_da)
        z_block = sum_series(cz, *z_da)
        fixed_z = [Z_COEFF_FIXED * v for v in z_block]
        fixed_xz = [X_COEFF_FIXED * x + Z_COEFF_FIXED * z for x, z in zip(x_block, z_block)]

        for window_name, (t_lo, t_hi) in windows.items():
            y = select_window(t_r, residual, t_lo, t_hi, include_right=(window_name == "late"))
            z_vals = select_window(t_r, fixed_z, t_lo, t_hi, include_right=(window_name == "late"))
            xz_vals = select_window(t_r, fixed_xz, t_lo, t_hi, include_right=(window_name == "late"))
            _, y_rms, _ = finite_stats(y)
            for label, vals in {
                "z_block_fixed": z_vals,
                "xz_block_fixed": xz_vals,
            }.items():
                corrected = [yy - vv for yy, vv in zip(y, vals)]
                _, corr_rms, _ = finite_stats(corrected)
                reduction = y_rms / corr_rms if corr_rms > 0 else math.inf
                rows.append(
                    {
                        "case_name": case_dir.name,
                        "window": window_name,
                        "candidate": label,
                        "residual_rms_abs_rate": y_rms,
                        "corr_with_residual": pearson(vals, y),
                        "implied_scale": fit_single(y, vals),
                        "corrected_rms_abs_rate": corr_rms,
                        "rms_reduction_factor": reduction,
                        "z_coeff_fixed": Z_COEFF_FIXED,
                        "x_coeff_fixed": X_COEFF_FIXED if label == "xz_block_fixed" else 0.0,
                    }
                )

    if not rows:
        raise SystemExit("No Step-201 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-201 Post-fix Z-Block Fixed-Rule Validation\n\n")
        f.write(f"- fixed_z_coeff: `{fmt(Z_COEFF_FIXED)}`\n")
        f.write(f"- fixed_x_coeff: `{fmt(X_COEFF_FIXED)}`\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            f.write(f"## {case_name}\n\n")
            f.write("| window | candidate | residual_rms | corr | implied_scale | corrected_rms | reduction |\n")
            f.write("|:---|:---|---:|---:|---:|---:|---:|\n")
            ranked = [r for r in rows if r["case_name"] == case_name]
            for r in ranked:
                f.write(
                    f"| {r['window']} | {r['candidate']} | {fmt(r['residual_rms_abs_rate'])} | "
                    f"{fmt(r['corr_with_residual'])} | {fmt(r['implied_scale'])} | {fmt(r['corrected_rms_abs_rate'])} | "
                    f"{fmt(r['rms_reduction_factor'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
