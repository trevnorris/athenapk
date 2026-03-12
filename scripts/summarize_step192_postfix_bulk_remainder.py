#!/usr/bin/env python3
"""Decompose the post-fix bulk residual remainder after exact Cx/Cz subtraction."""

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
    return [v for t, v in zip(times, values) if t >= t_lo and (t <= t_hi if include_right else t < t_hi)]


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return math.nan, math.nan, math.nan
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


def fit_two(y, x1, x2):
    s11 = sum(a * a for a in x1)
    s22 = sum(b * b for b in x2)
    s12 = sum(a * b for a, b in zip(x1, x2))
    t1 = sum(a * yy for a, yy in zip(x1, y))
    t2 = sum(b * yy for b, yy in zip(x2, y))
    det = s11 * s22 - s12 * s12
    if abs(det) <= 1.0e-30:
        return math.nan, math.nan
    a = (t1 * s22 - t2 * s12) / det
    b = (s11 * t2 - s12 * t1) / det
    return a, b


def fit_single(y, x):
    denom = sum(v * v for v in x)
    if denom <= 1.0e-30:
        return math.nan
    return sum(a * b for a, b in zip(y, x)) / denom


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
            "m4d_mixed_ew2",
            "m4d_mixed_c2",
            "m4d_mixed_cx2",
            "m4d_mixed_cy2",
            "m4d_mixed_cz2",
            "m4d_em_u_resolved",
        ]
        if not all(k in cols for k in required):
            continue

        t_r, du = derivative_series(cols["time"], cols["m4d_em_u_bulk"])
        _, dja = derivative_series(cols["time"], cols["m4d_int_ja_ea"])
        _, djw = derivative_series(cols["time"], cols["m4d_int_jw_ew"])
        _, dcx = derivative_series(cols["time"], cols["m4d_int_cx_energy_delta_exact_discrete_sum"])
        _, dcz = derivative_series(cols["time"], cols["m4d_int_cz_energy_delta_exact_discrete_sum"])
        residual = [u + a + w for u, a, w in zip(du, dja, djw)]
        if not t_r:
            continue

        start, end = t_r[0], t_r[-1]
        t_lo = start + (end - start) * 2.0 / 3.0
        t_hi = end
        y = select_window(t_r, residual, t_lo, t_hi, include_right=True)
        x = select_window(t_r, dcx, t_lo, t_hi, include_right=True)
        z = select_window(t_r, dcz, t_lo, t_hi, include_right=True)
        coeff_x, coeff_z = fit_two(y, x, z)
        remainder = [yy - coeff_x * xx - coeff_z * zz for yy, xx, zz in zip(y, x, z)]
        rem_max, rem_rms, rem_final = finite_stats(remainder)
        rows.append(
            {
                "case_name": case_dir.name,
                "component": "postfit_remainder",
                "fit_coeff": 1.0,
                "corr_with_remainder": 1.0,
                "rms_abs_rate": rem_rms,
                "max_abs_rate": rem_max,
                "final_rate": rem_final,
                "corrected_rms_abs_rate": 0.0,
                "rms_reduction_factor": math.inf,
                "coeff_x": coeff_x,
                "coeff_z": coeff_z,
            }
        )

        candidate_map = {
            "d_mixed_ew2_dt": derivative_series(cols["time"], cols["m4d_mixed_ew2"])[1],
            "d_mixed_c2_dt": derivative_series(cols["time"], cols["m4d_mixed_c2"])[1],
            "d_mixed_cx2_dt": derivative_series(cols["time"], cols["m4d_mixed_cx2"])[1],
            "d_mixed_cy2_dt": derivative_series(cols["time"], cols["m4d_mixed_cy2"])[1],
            "d_mixed_cz2_dt": derivative_series(cols["time"], cols["m4d_mixed_cz2"])[1],
            "d_em_u_resolved_dt": derivative_series(cols["time"], cols["m4d_em_u_resolved"])[1],
        }
        for name, series in candidate_map.items():
            window_vals = select_window(t_r, series, t_lo, t_hi, include_right=True)
            coeff = fit_single(remainder, window_vals)
            corrected = [yy - coeff * xx for yy, xx in zip(remainder, window_vals)]
            _, corrected_rms, _ = finite_stats(corrected)
            reduction = rem_rms / corrected_rms if corrected_rms > 0 else math.inf
            max_abs, rms_abs, final_rate = finite_stats(window_vals)
            rows.append(
                {
                    "case_name": case_dir.name,
                    "component": name,
                    "fit_coeff": coeff,
                    "corr_with_remainder": pearson(remainder, window_vals),
                    "rms_abs_rate": rms_abs,
                    "max_abs_rate": max_abs,
                    "final_rate": final_rate,
                    "corrected_rms_abs_rate": corrected_rms,
                    "rms_reduction_factor": reduction,
                    "coeff_x": coeff_x,
                    "coeff_z": coeff_z,
                }
            )

    if not rows:
        raise SystemExit("No Step-192 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-192 Post-Fix Bulk Remainder Decomposition\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for case_name in sorted({r['case_name'] for r in rows}):
            grp = [r for r in rows if r["case_name"] == case_name]
            rem = next(r for r in grp if r["component"] == "postfit_remainder")
            f.write(f"## {case_name} window=late\n\n")
            f.write(f"- postfit_remainder_rms_abs_rate: `{fmt(rem['rms_abs_rate'])}`\n")
            f.write(f"- fitted_exact_cx_coeff: `{fmt(rem['coeff_x'])}`\n")
            f.write(f"- fitted_exact_cz_coeff: `{fmt(rem['coeff_z'])}`\n\n")
            f.write("| component | corr_with_remainder | fit_coeff | corrected_rms_abs_rate | rms_reduction_factor | rms_abs_rate | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|\n")
            ranked = [r for r in grp if r["component"] != "postfit_remainder"]
            ranked.sort(key=lambda r: r["rms_reduction_factor"], reverse=True)
            for r in ranked:
                f.write(
                    f"| {r['component']} | {fmt(r['corr_with_remainder'])} | {fmt(r['fit_coeff'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} | {fmt(r['rms_abs_rate'])} | {fmt(r['max_abs_rate'])} | {fmt(r['final_rate'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
