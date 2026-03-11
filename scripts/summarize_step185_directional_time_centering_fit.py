#!/usr/bin/env python3
"""Fit directional x/z A_w update-power channels with independent temporal alignment choices."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from itertools import product
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


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    if v == math.inf:
        return "inf"
    return f"{float(v):.6e}"


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


def safe_ratio(num, den):
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def fit_multi(y, xs):
    mat = [[sum(a * b for a, b in zip(xi, xj)) for xj in xs] for xi in xs]
    vec = [sum(a * b for a, b in zip(y, xi)) for xi in xs]
    n = len(vec)
    aug = [row[:] + [vec[i]] for i, row in enumerate(mat)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1.0e-30:
            return None
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        piv = aug[col][col]
        for j in range(col, n + 1):
            aug[col][j] /= piv
        for r in range(n):
            if r == col:
                continue
            fac = aug[r][col]
            if fac == 0.0:
                continue
            for j in range(col, n + 1):
                aug[r][j] -= fac * aug[col][j]
    return [aug[i][n] for i in range(n)]


def window_bounds(times, fractions):
    start = times[0]
    end = times[-1]
    span = end - start
    if span <= 0.0:
        return []
    names = ["early", "mid", "late"]
    out = []
    for idx in range(len(fractions) - 1):
        f_lo = fractions[idx]
        f_hi = fractions[idx + 1]
        name = names[idx] if idx < len(names) else f"w{idx}"
        t_lo = start + f_lo * span
        t_hi = start + f_hi * span
        include_right = idx == len(fractions) - 2
        out.append((name, t_lo, t_hi, include_right))
    return out


def shift_prev(values):
    return [math.nan] + values[:-1]


def shift_next(values):
    return values[1:] + [math.nan]


def midpoint_prev(values):
    out = [math.nan]
    out.extend(0.5 * (values[i - 1] + values[i]) for i in range(1, len(values)))
    return out


def midpoint_next(values):
    out = [0.5 * (values[i] + values[i + 1]) for i in range(len(values) - 1)]
    out.append(math.nan)
    return out


def paired_window(times, target, x_series, z_series, t_lo, t_hi, include_right=False):
    y_out, x_out, z_out = [], [], []
    for t, y, x, z in zip(times, target, x_series, z_series):
        if not (t >= t_lo and (t <= t_hi if include_right else t < t_hi)):
            continue
        if not (math.isfinite(y) and math.isfinite(x) and math.isfinite(z)):
            continue
        y_out.append(y)
        x_out.append(x)
        z_out.append(z)
    return y_out, x_out, z_out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", required=True)
    ap.add_argument("--baseline-case", default="baseline_quadrature")
    ap.add_argument("--window-fractions", default="0,0.333333,0.666667,1")
    ap.add_argument("--focus-window", default="late", choices=["early", "mid", "late", "all"])
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
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

    required = [
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
        raise SystemExit("Missing Step-183 A_w discrete update-power channels")

    x_current = [
        sum(vs) for vs in zip(
            derivative_series(times, cols["m4d_int_aw_mode0_gradx_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode1_gradx_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode2_gradx_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode3_gradx_power_da_dt_discrete"])[1],
        )
    ]
    z_current = [
        sum(vs) for vs in zip(
            derivative_series(times, cols["m4d_int_aw_mode0_gradz_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode1_gradz_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode2_gradz_power_da_dt_discrete"])[1],
            derivative_series(times, cols["m4d_int_aw_mode3_gradz_power_da_dt_discrete"])[1],
        )
    ]

    if len(x_current) != len(residual) or len(z_current) != len(residual):
        raise SystemExit("Series length mismatch")

    x_variants = {
        "x_current": x_current,
        "x_prev": shift_prev(x_current),
        "x_next": shift_next(x_current),
        "x_mid_prev": midpoint_prev(x_current),
        "x_mid_next": midpoint_next(x_current),
    }
    z_variants = {
        "z_current": z_current,
        "z_prev": shift_prev(z_current),
        "z_next": shift_next(z_current),
        "z_mid_prev": midpoint_prev(z_current),
        "z_mid_next": midpoint_next(z_current),
    }

    rows = []
    for window, t_lo, t_hi, include_right in window_bounds(t_r, fractions):
        if args.focus_window != "all" and window != args.focus_window:
            continue
        res_window = [
            v for t, v in zip(t_r, residual)
            if t >= t_lo and (t <= t_hi if include_right else t < t_hi) and math.isfinite(v)
        ]
        if not res_window:
            continue
        res_max, res_rms, _ = finite_stats(res_window)
        for x_name, x_series in x_variants.items():
            for z_name, z_series in z_variants.items():
                y_win, x_win, z_win = paired_window(t_r, residual, x_series, z_series, t_lo, t_hi, include_right)
                if len(y_win) < 2:
                    continue
                coeffs = fit_multi(y_win, [z_win, x_win])
                if coeffs is None:
                    continue
                pred = [coeffs[0] * z + coeffs[1] * x for z, x in zip(z_win, x_win)]
                corrected = [y - p for y, p in zip(y_win, pred)]
                _, corrected_rms, corrected_final = finite_stats(corrected)
                pred_max, pred_rms, pred_final = finite_stats(pred)
                rows.append({
                    "case_name": case_dir.name,
                    "window": window,
                    "z_variant": z_name,
                    "x_variant": x_name,
                    "samples": len(y_win),
                    "residual_rms_abs_rate": res_rms,
                    "residual_max_abs_rate": res_max,
                    "predicted_rms_abs_rate": pred_rms,
                    "predicted_max_abs_rate": pred_max,
                    "corr_with_residual": pearson(pred, y_win),
                    "coeff_z": coeffs[0],
                    "coeff_x": coeffs[1],
                    "corrected_rms_abs_rate": corrected_rms,
                    "corrected_final_rate": corrected_final,
                    "rms_reduction_factor": safe_ratio(res_rms, corrected_rms),
                })

    if not rows:
        raise SystemExit("No Step-185 rows were computed")

    rows.sort(key=lambda r: (r["window"], -(r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf), r["z_variant"], r["x_variant"]))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-185 Directional Time-Centering Fit\n\n")
        f.write(f"Baseline case: `{case_dir.name}`.\n\n")
        f.write(f"Controlled psi0 span: `{safe_float(controlled.get('final_psi0_span', 'nan')):.6e}`\n")
        f.write(f"Full psi0 span: `{safe_float(full.get('final_psi0_span', 'nan')):.6e}`\n\n")
        windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
        for window in windows:
            grp = [r for r in rows if r["window"] == window]
            if not grp:
                continue
            f.write(f"## window={window}\n\n")
            f.write(f"- residual_rms_abs_rate: `{fmt(grp[0]['residual_rms_abs_rate'])}`\n")
            f.write(f"- residual_max_abs_rate: `{fmt(grp[0]['residual_max_abs_rate'])}`\n\n")
            f.write("| z_variant | x_variant | samples | corr_with_residual | coeff_z | coeff_x | corrected_rms_abs_rate | rms_reduction_factor |\n")
            f.write("|:---|:---|---:|---:|---:|---:|---:|---:|\n")
            for r in grp[:15]:
                f.write(
                    f"| {r['z_variant']} | {r['x_variant']} | {r['samples']} | {fmt(r['corr_with_residual'])} | "
                    f"{fmt(r['coeff_z'])} | {fmt(r['coeff_x'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
