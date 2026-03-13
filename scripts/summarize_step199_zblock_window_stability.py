#!/usr/bin/env python3
"""Assess stability of the post-fix z-block fit across time windows."""
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


def finite_rms(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return math.nan
    return math.sqrt(sum(v * v for v in finite) / len(finite))


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


def solve_linear_system(mat, rhs):
    n = len(rhs)
    a = [row[:] + [rhs_i] for row, rhs_i in zip(mat, rhs)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-30:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        for j in range(col, n + 1):
            a[col][j] /= scale
        for row in range(n):
            if row == col:
                continue
            factor = a[row][col]
            if factor == 0.0:
                continue
            for j in range(col, n + 1):
                a[row][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def fit_single(y, x):
    xx = sum(v * v for v in x)
    alpha = sum(a * b for a, b in zip(y, x)) / xx if xx > 0 else math.nan
    fitted = [alpha * v for v in x] if math.isfinite(alpha) else []
    return alpha, fitted


def fit_multi(y, comps):
    n = len(comps)
    gram = [[0.0 for _ in range(n)] for _ in range(n)]
    rhs = [0.0 for _ in range(n)]
    for i in range(n):
        for j in range(n):
            gram[i][j] = sum(a * b for a, b in zip(comps[i], comps[j]))
        rhs[i] = sum(a * b for a, b in zip(comps[i], y))
    coeffs = solve_linear_system(gram, rhs)
    if coeffs is None:
        return None, []
    fitted = [sum(c * comps[k][idx] for k, c in enumerate(coeffs)) for idx in range(len(y))]
    return coeffs, fitted


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    if v == math.inf:
        return "inf"
    return f"{float(v):.6e}"


def series_sum(*series):
    return [sum(vals) for vals in zip(*series)]


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

        late_z_alpha = None
        for window_name, (t_lo, t_hi) in windows.items():
            y = select_window(t_r, residual, t_lo, t_hi, include_right=True)
            res_rms = finite_rms(y)
            cx_w = select_window(t_r, cx, t_lo, t_hi, include_right=True)
            cz_w = select_window(t_r, cz, t_lo, t_hi, include_right=True)
            x_da_w = [select_window(t_r, d, t_lo, t_hi, include_right=True) for d in x_da]
            z_da_w = [select_window(t_r, d, t_lo, t_hi, include_right=True) for d in z_da]
            x_da_sum = series_sum(*x_da_w)
            z_da_sum = series_sum(*z_da_w)
            x_block = series_sum(cx_w, x_da_sum)
            z_block = series_sum(cz_w, z_da_sum)

            z_alpha, z_fit = fit_single(y, z_block)
            z_corr_rms = finite_rms([yy - ff for yy, ff in zip(y, z_fit)])
            if window_name == "late":
                late_z_alpha = z_alpha

            coeffs, fit = fit_multi(y, [x_block, z_block])
            fit_rms = finite_rms([yy - ff for yy, ff in zip(y, fit)]) if fit else math.nan

            fixed_late_rms = math.nan
            fixed_late_corr = math.nan
            if late_z_alpha is not None:
                fixed = [late_z_alpha * v for v in z_block]
                fixed_late_rms = finite_rms([yy - ff for yy, ff in zip(y, fixed)])
                fixed_late_corr = pearson(fixed, y)

            rows.append(
                {
                    "case_name": case_dir.name,
                    "window": window_name,
                    "residual_rms": res_rms,
                    "z_block_corr": pearson(z_block, y),
                    "z_block_fit_coeff": z_alpha,
                    "z_block_corrected_rms": z_corr_rms,
                    "xz_twofit_corr": pearson(fit, y) if fit else math.nan,
                    "xz_twofit_coeff_x": coeffs[0] if coeffs else math.nan,
                    "xz_twofit_coeff_z": coeffs[1] if coeffs else math.nan,
                    "xz_twofit_corrected_rms": fit_rms,
                    "late_z_coeff_applied_corr": fixed_late_corr,
                    "late_z_coeff_applied_rms": fixed_late_rms,
                }
            )

    if not rows:
        raise SystemExit("No Step-199 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r.keys()}))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-199 Z-block Window Stability\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            f.write(f"## {case_name}\n\n")
            f.write("| window | residual rms | z corr | z coeff | z corrected rms | xz corr | x coeff | z coeff | xz corrected rms | late-z corr | late-z corrected rms |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in rows:
                if r["case_name"] != case_name:
                    continue
                f.write(
                    f"| {r['window']} | {fmt(r['residual_rms'])} | {fmt(r['z_block_corr'])} | {fmt(r['z_block_fit_coeff'])} | {fmt(r['z_block_corrected_rms'])} | {fmt(r['xz_twofit_corr'])} | {fmt(r['xz_twofit_coeff_x'])} | {fmt(r['xz_twofit_coeff_z'])} | {fmt(r['xz_twofit_corrected_rms'])} | {fmt(r['late_z_coeff_applied_corr'])} | {fmt(r['late_z_coeff_applied_rms'])} |\n"
                )


if __name__ == "__main__":
    main()
