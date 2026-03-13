#!/usr/bin/env python3
"""Fit the post-fix bulk residual against exact primitive mixed-block pieces."""
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


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    return (max(abs(v) for v in finite), math.sqrt(sum(v * v for v in finite) / len(finite)), finite[-1])


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


def sum_series(*series_list):
    return [sum(vals) for vals in zip(*series_list)]


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
    block_rows = []

    primitive_cols = {
        "cx_daw_dx2": "m4d_int_cx_daw_dx2_energy_delta_exact_discrete_sum",
        "cx_dw_ax2": "m4d_int_cx_dw_ax2_energy_delta_exact_discrete_sum",
        "cx_cross": "m4d_int_cx_cross_energy_delta_exact_discrete_sum",
        "cz_daw_dz2": "m4d_int_cz_daw_dz2_energy_delta_exact_discrete_sum",
        "cz_dw_az2": "m4d_int_cz_dw_az2_energy_delta_exact_discrete_sum",
        "cz_cross": "m4d_int_cz_cross_energy_delta_exact_discrete_sum",
    }

    for case_dir in case_dirs:
        hst = case_dir / "harris_full.out1.hst"
        if not hst.exists():
            continue
        cols = hsm.parse_hst(hst)
        required = ["time", "m4d_em_u_bulk", "m4d_int_ja_ea", "m4d_int_jw_ew", *primitive_cols.values()]
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

        comp_series = {}
        comp_rms = {}
        for name, col in primitive_cols.items():
            _, deriv = derivative_series(times, cols[col])
            comp_series[name] = select_window(t_r, deriv, t_lo, end, include_right=True)
            _, rms_abs, _ = finite_stats(comp_series[name])
            comp_rms[name] = rms_abs
            alpha = fit_single(y, comp_series[name])
            corrected = [yy - alpha * vv for yy, vv in zip(y, comp_series[name])] if math.isfinite(alpha) else []
            _, corrected_rms, _ = finite_stats(corrected)
            reduction = res_rms / corrected_rms if corrected_rms > 0 else math.inf
            rows.append({
                "case_name": case_dir.name,
                "candidate": name,
                "corr_with_residual": pearson(comp_series[name], y),
                "fit_coeff": alpha,
                "rms_abs_rate": rms_abs,
                "corrected_rms_abs_rate": corrected_rms,
                "rms_reduction_factor": reduction,
            })

        blocks = {
            "cx_block_exact": ["cx_daw_dx2", "cx_dw_ax2", "cx_cross"],
            "cz_block_exact": ["cz_daw_dz2", "cz_dw_az2", "cz_cross"],
            "cxz_block_exact": [
                "cx_daw_dx2", "cx_dw_ax2", "cx_cross",
                "cz_daw_dz2", "cz_dw_az2", "cz_cross",
            ],
        }
        for block_name, labels in blocks.items():
            active_labels = [
                label for label in labels if math.isfinite(comp_rms[label]) and comp_rms[label] > 1.0e-12
            ]
            coeffs = None
            fitted = []
            if active_labels:
                coeffs, fitted = fit_multi(y, [comp_series[label] for label in active_labels])
            _, corrected_rms, _ = finite_stats([yy - ff for yy, ff in zip(y, fitted)])
            reduction = res_rms / corrected_rms if corrected_rms > 0 else math.inf
            row = {
                "case_name": case_dir.name,
                "block": block_name,
                "corr_with_residual": pearson(fitted, y) if fitted else math.nan,
                "corrected_rms_abs_rate": corrected_rms,
                "rms_reduction_factor": reduction,
            }
            if coeffs is None:
                row["fit_status"] = "singular"
            else:
                row["fit_status"] = "fit"
                for label, coeff in zip(active_labels, coeffs):
                    row[f"coeff_{label}"] = coeff
            block_rows.append(row)

    if not rows or not block_rows:
        raise SystemExit("No Step-196 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    block_csv = out_csv.with_name(out_csv.stem + "_blocks.csv")
    with block_csv.open("w", newline="", encoding="utf-8") as f:
        fieldnames = sorted({k for row in block_rows for k in row.keys()})
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(block_rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-196 Mixed-Block Exact Probe\n\n")
        for case_name in sorted({r['case_name'] for r in rows}):
            f.write(f"## {case_name} window=late\n\n")
            f.write("### Primitive fits\n\n")
            f.write("| candidate | corr | fit coeff | corrected rms | reduction |\n")
            f.write("|:---|---:|---:|---:|---:|\n")
            ranked = [r for r in rows if r['case_name'] == case_name]
            ranked.sort(key=lambda r: r['rms_reduction_factor'], reverse=True)
            for r in ranked:
                f.write(f"| {r['candidate']} | {fmt(r['corr_with_residual'])} | {fmt(r['fit_coeff'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} |\n")
            f.write("\n### Block fits\n\n")
            for r in [row for row in block_rows if row['case_name'] == case_name]:
                f.write(f"- `{r['block']}` corr=`{fmt(r['corr_with_residual'])}` corrected_rms=`{fmt(r['corrected_rms_abs_rate'])}` reduction=`{fmt(r['rms_reduction_factor'])}`\n")
                coeff_items = sorted((k, v) for k, v in r.items() if k.startswith('coeff_'))
                for k, v in coeff_items:
                    f.write(f"  - `{k}` = `{fmt(v)}`\n")
                f.write("\n")
        f.write(f"Additional CSV: `{block_csv}`\n")

if __name__ == "__main__":
    main()
