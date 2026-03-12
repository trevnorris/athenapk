#!/usr/bin/env python3
"""Compare continuum mixed Cx/Cz bulk-energy bookkeeping against exact-current host bookkeeping."""

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
            "m4d_em_u_bulk_exact_current_cxz",
            "m4d_int_ja_ea",
            "m4d_int_jw_ew",
            "m4d_mixed_cx2",
            "m4d_mixed_cz2",
            "m4d_int_cx_energy_exact_current_sum",
            "m4d_int_cz_energy_exact_current_sum",
        ]
        if not all(k in cols for k in required):
            continue

        times = cols["time"]
        t_r, du_bulk = derivative_series(times, cols["m4d_em_u_bulk"])
        _, du_bulk_exact = derivative_series(times, cols["m4d_em_u_bulk_exact_current_cxz"])
        _, dja = derivative_series(times, cols["m4d_int_ja_ea"])
        _, djw = derivative_series(times, cols["m4d_int_jw_ew"])
        _, d_cx_mixed = derivative_series(times, cols["m4d_mixed_cx2"])
        _, d_cz_mixed = derivative_series(times, cols["m4d_mixed_cz2"])
        _, d_cx_exact = derivative_series(times, cols["m4d_int_cx_energy_exact_current_sum"])
        _, d_cz_exact = derivative_series(times, cols["m4d_int_cz_energy_exact_current_sum"])
        if not t_r:
            continue

        residual = [u + a + w for u, a, w in zip(du_bulk, dja, djw)]
        residual_exact_bulk = [u + a + w for u, a, w in zip(du_bulk_exact, dja, djw)]
        bulk_gap = [new - old for old, new in zip(cols["m4d_em_u_bulk"], cols["m4d_em_u_bulk_exact_current_cxz"])]
        cx_state_gap = [new - old for old, new in zip(cols["m4d_mixed_cx2"], cols["m4d_int_cx_energy_exact_current_sum"])]
        cz_state_gap = [new - old for old, new in zip(cols["m4d_mixed_cz2"], cols["m4d_int_cz_energy_exact_current_sum"])]
        _, d_bulk_gap = derivative_series(times, bulk_gap)
        _, d_cx_gap = derivative_series(times, cx_state_gap)
        _, d_cz_gap = derivative_series(times, cz_state_gap)

        start, end = t_r[0], t_r[-1]
        t_lo = start + (end - start) * 2.0 / 3.0
        include_right = True
        y = select_window(t_r, residual, t_lo, end, include_right)
        y_exact = select_window(t_r, residual_exact_bulk, t_lo, end, include_right)
        bulk_gap_w = select_window(times, bulk_gap, t_lo, end, include_right)
        cx_gap_w = select_window(times, cx_state_gap, t_lo, end, include_right)
        cz_gap_w = select_window(times, cz_state_gap, t_lo, end, include_right)
        d_bulk_gap_w = select_window(t_r, d_bulk_gap, t_lo, end, include_right)
        d_cx_gap_w = select_window(t_r, d_cx_gap, t_lo, end, include_right)
        d_cz_gap_w = select_window(t_r, d_cz_gap, t_lo, end, include_right)
        d_cx_exact_w = select_window(t_r, d_cx_exact, t_lo, end, include_right)
        d_cz_exact_w = select_window(t_r, d_cz_exact, t_lo, end, include_right)
        d_cx_mixed_w = select_window(t_r, d_cx_mixed, t_lo, end, include_right)
        d_cz_mixed_w = select_window(t_r, d_cz_mixed, t_lo, end, include_right)

        base_max, base_rms, base_final = finite_stats(y)
        exact_max, exact_rms, exact_final = finite_stats(y_exact)
        rows.extend([
            {
                "case_name": case_dir.name,
                "component": "baseline_residual",
                "max_abs": base_max,
                "rms": base_rms,
                "final": base_final,
                "corr_with_baseline": 1.0,
            },
            {
                "case_name": case_dir.name,
                "component": "exact_current_bulk_residual",
                "max_abs": exact_max,
                "rms": exact_rms,
                "final": exact_final,
                "corr_with_baseline": pearson(y_exact, y),
            },
            {
                "case_name": case_dir.name,
                "component": "bulk_state_gap",
                "max_abs": finite_stats(bulk_gap_w)[0],
                "rms": finite_stats(bulk_gap_w)[1],
                "final": finite_stats(bulk_gap_w)[2],
                "corr_with_baseline": math.nan,
            },
            {
                "case_name": case_dir.name,
                "component": "cx_state_gap",
                "max_abs": finite_stats(cx_gap_w)[0],
                "rms": finite_stats(cx_gap_w)[1],
                "final": finite_stats(cx_gap_w)[2],
                "corr_with_baseline": math.nan,
            },
            {
                "case_name": case_dir.name,
                "component": "cz_state_gap",
                "max_abs": finite_stats(cz_gap_w)[0],
                "rms": finite_stats(cz_gap_w)[1],
                "final": finite_stats(cz_gap_w)[2],
                "corr_with_baseline": math.nan,
            },
            {
                "case_name": case_dir.name,
                "component": "d_bulk_gap_dt",
                "max_abs": finite_stats(d_bulk_gap_w)[0],
                "rms": finite_stats(d_bulk_gap_w)[1],
                "final": finite_stats(d_bulk_gap_w)[2],
                "corr_with_baseline": pearson(d_bulk_gap_w, y),
            },
            {
                "case_name": case_dir.name,
                "component": "d_cx_gap_dt",
                "max_abs": finite_stats(d_cx_gap_w)[0],
                "rms": finite_stats(d_cx_gap_w)[1],
                "final": finite_stats(d_cx_gap_w)[2],
                "corr_with_baseline": pearson(d_cx_gap_w, y),
            },
            {
                "case_name": case_dir.name,
                "component": "d_cz_gap_dt",
                "max_abs": finite_stats(d_cz_gap_w)[0],
                "rms": finite_stats(d_cz_gap_w)[1],
                "final": finite_stats(d_cz_gap_w)[2],
                "corr_with_baseline": pearson(d_cz_gap_w, y),
            },
            {
                "case_name": case_dir.name,
                "component": "d_cx_exact_current_dt",
                "max_abs": finite_stats(d_cx_exact_w)[0],
                "rms": finite_stats(d_cx_exact_w)[1],
                "final": finite_stats(d_cx_exact_w)[2],
                "corr_with_baseline": pearson(d_cx_exact_w, y),
            },
            {
                "case_name": case_dir.name,
                "component": "d_cz_exact_current_dt",
                "max_abs": finite_stats(d_cz_exact_w)[0],
                "rms": finite_stats(d_cz_exact_w)[1],
                "final": finite_stats(d_cz_exact_w)[2],
                "corr_with_baseline": pearson(d_cz_exact_w, y),
            },
            {
                "case_name": case_dir.name,
                "component": "d_cx_mixed_dt",
                "max_abs": finite_stats(d_cx_mixed_w)[0],
                "rms": finite_stats(d_cx_mixed_w)[1],
                "final": finite_stats(d_cx_mixed_w)[2],
                "corr_with_baseline": pearson(d_cx_mixed_w, y),
            },
            {
                "case_name": case_dir.name,
                "component": "d_cz_mixed_dt",
                "max_abs": finite_stats(d_cz_mixed_w)[0],
                "rms": finite_stats(d_cz_mixed_w)[1],
                "final": finite_stats(d_cz_mixed_w)[2],
                "corr_with_baseline": pearson(d_cz_mixed_w, y),
            },
        ])

    if not rows:
        raise SystemExit("No Step-193 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-193 Exact-Current Cx/Cz Equivalence Probe\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for case_name in sorted({r['case_name'] for r in rows}):
            grp = [r for r in rows if r['case_name'] == case_name]
            f.write(f"## {case_name} window=late\n\n")
            f.write("| component | rms | max_abs | final | corr_with_baseline |\n")
            f.write("|:---|---:|---:|---:|---:|\n")
            order = {
                "baseline_residual": 0,
                "exact_current_bulk_residual": 1,
                "bulk_state_gap": 2,
                "cx_state_gap": 3,
                "cz_state_gap": 4,
                "d_bulk_gap_dt": 5,
                "d_cx_gap_dt": 6,
                "d_cz_gap_dt": 7,
                "d_cx_exact_current_dt": 8,
                "d_cz_exact_current_dt": 9,
                "d_cx_mixed_dt": 10,
                "d_cz_mixed_dt": 11,
            }
            grp.sort(key=lambda r: order.get(r["component"], 99))
            for r in grp:
                f.write(
                    f"| {r['component']} | {fmt(r['rms'])} | {fmt(r['max_abs'])} | {fmt(r['final'])} | {fmt(r['corr_with_baseline'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
