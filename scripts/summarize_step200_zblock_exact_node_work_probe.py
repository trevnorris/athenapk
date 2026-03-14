#!/usr/bin/env python3
"""Compare exact Cx/Cz energy deltas against exact-node A_w gradient work channels."""

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
            "m4d_int_aw_gradx_work_exact_node_sum",
            "m4d_int_aw_gradz_work_exact_node_sum",
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

        _, cx = derivative_series(times, cols["m4d_int_cx_energy_delta_exact_discrete_sum"])
        _, cz = derivative_series(times, cols["m4d_int_cz_energy_delta_exact_discrete_sum"])
        _, xwork = derivative_series(times, cols["m4d_int_aw_gradx_work_exact_node_sum"])
        _, zwork = derivative_series(times, cols["m4d_int_aw_gradz_work_exact_node_sum"])

        cx = select_window(t_r, cx, t_lo, end, include_right=True)
        cz = select_window(t_r, cz, t_lo, end, include_right=True)
        xwork = select_window(t_r, xwork, t_lo, end, include_right=True)
        zwork = select_window(t_r, zwork, t_lo, end, include_right=True)

        xz_exact = sum_series(cx, cz)
        xz_work = sum_series(xwork, zwork)

        candidates = {
            "cx_exact_sum": cx,
            "cz_exact_sum": cz,
            "cxz_exact_sum": xz_exact,
            "xwork_exact_node_sum": xwork,
            "zwork_exact_node_sum": zwork,
            "xzwork_exact_node_sum": xz_work,
        }

        rows.append(
            {
                "case_name": case_dir.name,
                "candidate": "em_bulk_ledger_residual",
                "rms_abs_rate": res_rms,
                "corr_with_residual": 1.0,
                "corrected_rms_abs_rate": 0.0,
                "rms_reduction_factor": math.inf,
            }
        )
        for name, vals in candidates.items():
            _, rms_abs, _ = finite_stats(vals)
            alpha, fitted = fit_single(y, vals)
            corrected = [yy - ff for yy, ff in zip(y, fitted)] if fitted else []
            _, corrected_rms, _ = finite_stats(corrected)
            reduction = res_rms / corrected_rms if corrected_rms > 0 else math.inf
            rows.append(
                {
                    "case_name": case_dir.name,
                    "candidate": name,
                    "rms_abs_rate": rms_abs,
                    "corr_with_residual": pearson(vals, y),
                    "fit_coeff": alpha,
                    "corrected_rms_abs_rate": corrected_rms,
                    "rms_reduction_factor": reduction,
                }
            )

        combos = {
            "x_exact_plus_xwork_exact_node": [cx, xwork],
            "z_exact_plus_zwork_exact_node": [cz, zwork],
            "xz_exact_plus_xzwork_exact_node": [xz_exact, xz_work],
        }
        for name, comps in combos.items():
            coeffs, fitted = fit_multi(y, comps)
            corrected = [yy - ff for yy, ff in zip(y, fitted)] if fitted else []
            _, corrected_rms, _ = finite_stats(corrected)
            reduction = res_rms / corrected_rms if corrected_rms > 0 else math.inf
            row = {
                "case_name": case_dir.name,
                "candidate": name,
                "corr_with_residual": pearson(fitted, y) if fitted else math.nan,
                "corrected_rms_abs_rate": corrected_rms,
                "rms_reduction_factor": reduction,
            }
            if coeffs is not None:
                row["coeff_0"] = coeffs[0]
                row["coeff_1"] = coeffs[1]
            rows.append(row)

    if not rows:
        raise SystemExit("No Step-200 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r.keys()}))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-200 Z-Block Exact-Node Work Probe\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            target = next(
                r
                for r in rows
                if r["case_name"] == case_name and r["candidate"] == "em_bulk_ledger_residual"
            )
            f.write(f"## {case_name} window=late\n\n")
            f.write(f"- residual_rms_abs_rate: `{fmt(target['rms_abs_rate'])}`\n\n")
            f.write("| candidate | corr_with_residual | fit coeff / coeffs | corrected_rms_abs_rate | rms_reduction_factor |\n")
            f.write("|:---|---:|:---|---:|---:|\n")
            ranked = [r for r in rows if r["case_name"] == case_name and r["candidate"] != "em_bulk_ledger_residual"]
            ranked.sort(key=lambda r: r["rms_reduction_factor"], reverse=True)
            for r in ranked:
                coeff = ", ".join(
                    f"{k}={fmt(v)}" for k, v in sorted(r.items()) if k.startswith("coeff_")
                ) or fmt(r.get("fit_coeff"))
                f.write(
                    f"| {r['candidate']} | {fmt(r['corr_with_residual'])} | {coeff} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
