#!/usr/bin/env python3
"""Fit the surviving quadrature-baseline ledger residual against diagonal (d_z a_w)^2 channels."""

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
    out = []
    for t, v in zip(times, values):
        if t >= t_lo and (t <= t_hi if include_right else t < t_hi):
            out.append(v)
    return out


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


def fmt(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def safe_ratio(num, den):
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


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


def dot(xs, ys):
    return sum(x * y for x, y in zip(xs, ys))


def solve_linear_system(mat, vec):
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


def fit_single(y, x):
    xx = dot(x, x)
    if xx <= 0.0:
        return math.nan
    return dot(y, x) / xx


def fit_multi(y, xs):
    mat = [[dot(xi, xj) for xj in xs] for xi in xs]
    vec = [dot(y, xi) for xi in xs]
    return solve_linear_system(mat, vec)


def combine(coeffs, xs):
    return [sum(c * x[i] for c, x in zip(coeffs, xs)) for i in range(len(xs[0]))]


def parse_case_rows(run_log: Path):
    rows = []
    with run_log.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("case,"):
                header = line.split(",")
                continue
            if line.startswith("controlled,") or line.startswith("full,"):
                vals = line.split(",")
                if len(vals) == len(header):
                    rows.append(dict(zip(header, vals)))
    by_case = {row["case"]: row for row in rows}
    return by_case.get("controlled", {}), by_case.get("full", {})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True)
    parser.add_argument("--baseline-case", default="baseline_quadrature")
    parser.add_argument("--window-fractions", default="0,0.333333,0.666667,1")
    parser.add_argument("--focus-window", default="late", choices=["early", "mid", "late", "all"])
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

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
    residual = [u + a + w for u, a, w in zip(du, dja, djw)]

    series = {
        "mode_diag": derivative_series(times, cols.get("m4d_mixed_cz_daw_dz2_mode_diag"))[1],
        "mode0": derivative_series(times, cols.get("m4d_mixed_cz_daw_dz2_mode0"))[1],
        "mode1": derivative_series(times, cols.get("m4d_mixed_cz_daw_dz2_mode1"))[1],
        "mode2": derivative_series(times, cols.get("m4d_mixed_cz_daw_dz2_mode2"))[1],
        "mode3": derivative_series(times, cols.get("m4d_mixed_cz_daw_dz2_mode3"))[1],
    }
    if not t_r or any(len(v) != len(residual) for v in series.values()):
        raise SystemExit("Series length mismatch")

    rows = []
    for window, t_lo, t_hi, include_right in window_bounds(t_r, fractions):
        if args.focus_window != "all" and window != args.focus_window:
            continue
        y = select_window(t_r, residual, t_lo, t_hi, include_right)
        if not y:
            continue
        r_max, r_rms, r_final = finite_stats(y)
        xdiag = select_window(t_r, series["mode_diag"], t_lo, t_hi, include_right)
        x0 = select_window(t_r, series["mode0"], t_lo, t_hi, include_right)
        x1 = select_window(t_r, series["mode1"], t_lo, t_hi, include_right)
        x2 = select_window(t_r, series["mode2"], t_lo, t_hi, include_right)
        x3 = select_window(t_r, series["mode3"], t_lo, t_hi, include_right)

        candidates = []
        # fixed unit-coefficient probes
        candidates.append(("diag_unit", [1.0], xdiag))
        candidates.append(("modes12_unit", [1.0, 1.0], [x1, x2]))
        # least-squares fits
        a_diag = fit_single(y, xdiag)
        candidates.append(("diag_fit", [a_diag], xdiag))
        a12 = fit_multi(y, [x1, x2])
        if a12 is not None:
            candidates.append(("modes12_fit", a12, [x1, x2]))
        a4 = fit_multi(y, [x0, x1, x2, x3])
        if a4 is not None:
            candidates.append(("modes0123_fit", a4, [x0, x1, x2, x3]))

        for name, coeffs, xs in candidates:
            xs_list = xs if isinstance(xs, list) and xs and isinstance(xs[0], list) else [xs]
            pred = combine(coeffs, xs_list)
            corrected = [yv - pv for yv, pv in zip(y, pred)]
            p_max, p_rms, p_final = finite_stats(pred)
            c_max, c_rms, c_final = finite_stats(corrected)
            row = {
                "case_name": case_dir.name,
                "window": window,
                "candidate": name,
                "residual_rms_abs_rate": r_rms,
                "residual_max_abs_rate": r_max,
                "predicted_rms_abs_rate": p_rms,
                "predicted_max_abs_rate": p_max,
                "corr_with_residual": pearson(pred, y),
                "corrected_rms_abs_rate": c_rms,
                "corrected_max_abs_rate": c_max,
                "corrected_final_rate": c_final,
                "rms_reduction_factor": safe_ratio(r_rms, c_rms),
                "predicted_share_of_residual_rms": safe_ratio(p_rms, r_rms),
                "coeff_0": coeffs[0] if len(coeffs) > 0 else math.nan,
                "coeff_1": coeffs[1] if len(coeffs) > 1 else math.nan,
                "coeff_2": coeffs[2] if len(coeffs) > 2 else math.nan,
                "coeff_3": coeffs[3] if len(coeffs) > 3 else math.nan,
            }
            rows.append(row)

    if not rows:
        raise SystemExit("No rows computed")

    rows.sort(key=lambda r: (r["window"], -(r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf), r["candidate"]))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-172 d_z a_w^2 Bookkeeping Fit Probe\n\n")
        f.write(f"Baseline case: `{case_dir.name}`.\n\n")
        f.write(f"Controlled psi0 span: `{float(controlled.get('final_psi0_span', 'nan')):.6e}`\n")
        f.write(f"Full psi0 span: `{float(full.get('final_psi0_span', 'nan')):.6e}`\n\n")
        windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
        for window in windows:
            grp = [r for r in rows if r["window"] == window]
            if not grp:
                continue
            f.write(f"## window={window}\n\n")
            f.write(f"- residual_rms_abs_rate: `{fmt(grp[0]['residual_rms_abs_rate'])}`\n")
            f.write(f"- residual_max_abs_rate: `{fmt(grp[0]['residual_max_abs_rate'])}`\n\n")
            f.write("| candidate | corr_with_residual | predicted_share_of_residual_rms | coeff_0 | coeff_1 | coeff_2 | coeff_3 | corrected_rms_abs_rate | rms_reduction_factor | corrected_max_abs_rate | corrected_final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in grp:
                f.write(
                    f"| {r['candidate']} | {fmt(r['corr_with_residual'])} | {fmt(r['predicted_share_of_residual_rms'])} | {fmt(r['coeff_0'])} | {fmt(r['coeff_1'])} | {fmt(r['coeff_2'])} | {fmt(r['coeff_3'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} | {fmt(r['corrected_max_abs_rate'])} | {fmt(r['corrected_final_rate'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
