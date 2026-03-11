#!/usr/bin/env python3
"""Probe whether the surviving residual is primarily a sqrt(pi)*lambda weighting of the z channel."""

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


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    max_abs = max(abs(v) for v in finite)
    rms_abs = math.sqrt(sum(v * v for v in finite) / len(finite))
    final = finite[-1]
    return (max_abs, rms_abs, final)


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


def safe_ratio(num, den):
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def parse_case_rows(run_log: Path):
    rows = []
    header = None
    with run_log.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("case,"):
                header = line.split(",")
                continue
            if header and (line.startswith("controlled,") or line.startswith("full,")):
                vals = line.split(",")
                if len(vals) == len(header):
                    rows.append(dict(zip(header, vals)))
    by_case = {row["case"]: row for row in rows}
    return by_case.get("controlled", {}), by_case.get("full", {})


def safe_float(v: str) -> float:
    try:
        return float(v)
    except Exception:
        return math.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", required=True)
    ap.add_argument("--baseline-case", default="baseline_quadrature")
    ap.add_argument("--focus-window", default="late", choices=["early", "mid", "late"])
    ap.add_argument("--lambda", dest="lambda_value", type=float, default=1.0)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]
    case_dir = next((d for d in case_dirs if d.name == args.baseline_case), None)
    if case_dir is None:
        raise SystemExit(f"Baseline case not found: {args.baseline_case}")

    run_log = case_dir / "run.log"
    hst = case_dir / "harris_full.out1.hst"
    if not run_log.exists() or not hst.exists():
        raise SystemExit(f"Missing run.log or hst for {case_dir}")

    controlled, full = parse_case_rows(run_log)
    cols = hsm.parse_hst(hst)
    times = cols["time"]

    _, du = derivative_series(times, cols.get("m4d_em_u_bulk"))
    t_r, dja = derivative_series(times, cols.get("m4d_int_ja_ea"))
    _, djw = derivative_series(times, cols.get("m4d_int_jw_ew"))
    if not t_r or len(du) != len(dja) or len(du) != len(djw):
        raise SystemExit("Could not construct ledger residual series")
    residual = [u + a + w for u, a, w in zip(du, dja, djw)]

    x_series = [
        sum(vs) for vs in zip(
            derivative_series(times, cols["m4d_int_aw_mode0_gradx_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode1_gradx_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode2_gradx_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode3_gradx_power_da_dt_discrete"])[1],
        )
    ]
    z_series = [
        sum(vs) for vs in zip(
            derivative_series(times, cols["m4d_int_aw_mode0_gradz_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode1_gradz_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode2_gradz_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode3_gradz_power_da_dt_discrete"])[1],
        )
    ]

    start, end = t_r[0], t_r[-1]
    if args.focus_window == "early":
        t_lo, t_hi = start, start + (end - start) / 3.0
    elif args.focus_window == "mid":
        t_lo, t_hi = start + (end - start) / 3.0, start + 2.0 * (end - start) / 3.0
    else:
        t_lo, t_hi = start + 2.0 * (end - start) / 3.0, end

    idx = [i for i, t in enumerate(t_r) if t >= t_lo and t <= t_hi]
    y = [residual[i] for i in idx]
    x = [x_series[i] for i in idx]
    z = [z_series[i] for i in idx]
    z_int = args.lambda_value * math.sqrt(math.pi)

    candidates = {
        "z_current": z,
        "z_current_scaled_zint": [z_int * v for v in z],
        "z_current_scaled_zint_minus_x_current": [(z_int * zv) - xv for zv, xv in zip(z, x)],
        "z_current_minus_x_current": [zv - xv for zv, xv in zip(z, x)],
    }

    rows = []
    res_max, res_rms, _ = finite_stats(y)
    for name, series in candidates.items():
        cand_max, cand_rms, cand_final = finite_stats(series)
        alpha = fit_single(y, series)
        corrected = [yy - alpha * xx for yy, xx in zip(y, series)] if math.isfinite(alpha) else []
        corr_max, corr_rms, corr_final = finite_stats(corrected)
        rows.append({
            "case_name": case_dir.name,
            "window": args.focus_window,
            "candidate": name,
            "z_int": z_int,
            "candidate_rms_abs_rate": cand_rms,
            "corr_with_residual": pearson(series, y),
            "fit_coeff": alpha,
            "corrected_rms_abs_rate": corr_rms,
            "rms_reduction_factor": safe_ratio(res_rms, corr_rms),
        })

    rows.sort(key=lambda r: -(r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-186 Z-Channel Normalization Probe\n\n")
        f.write(f"Baseline case: `{case_dir.name}`.\n\n")
        f.write(f"Controlled psi0 span: `{safe_float(controlled.get('final_psi0_span', 'nan')):.6e}`\n")
        f.write(f"Full psi0 span: `{safe_float(full.get('final_psi0_span', 'nan')):.6e}`\n")
        f.write(f"z_int = lambda * sqrt(pi): `{fmt(z_int)}`\n\n")
        f.write(f"## window={args.focus_window}\n\n")
        f.write(f"- residual_rms_abs_rate: `{fmt(res_rms)}`\n")
        f.write(f"- residual_max_abs_rate: `{fmt(res_max)}`\n\n")
        f.write("| candidate | corr_with_residual | fit_coeff | corrected_rms_abs_rate | rms_reduction_factor |\n")
        f.write("|:---|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(
                f"| {r['candidate']} | {fmt(r['corr_with_residual'])} | {fmt(r['fit_coeff'])} | "
                f"{fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} |\n"
            )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
