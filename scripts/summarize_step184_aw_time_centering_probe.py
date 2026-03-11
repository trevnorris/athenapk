#!/usr/bin/env python3
"""Test time-centering / stage-shift hypotheses for the full discrete A_w update-power sum."""

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


def fit_single(y, x):
    xx = sum(v * v for v in x)
    if xx <= 0.0:
        return math.nan
    return sum(a * b for a, b in zip(y, x)) / xx


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


def sum_series(*series_list):
    return [sum(vals) for vals in zip(*series_list)]


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


def paired_window(times, target, candidate, t_lo, t_hi, include_right=False):
    t_out = []
    y_out = []
    x_out = []
    for t, y, x in zip(times, target, candidate):
        if not (t >= t_lo and (t <= t_hi if include_right else t < t_hi)):
            continue
        if not (math.isfinite(y) and math.isfinite(x)):
            continue
        t_out.append(t)
        y_out.append(y)
        x_out.append(x)
    return t_out, y_out, x_out


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

    full_da_dt_xz = sum_series(
        derivative_series(times, cols["m4d_int_aw_mode0_gradx_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode1_gradx_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode2_gradx_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode3_gradx_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode0_gradz_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode1_gradz_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode2_gradz_power_da_dt_discrete"])[1],
        derivative_series(times, cols["m4d_int_aw_mode3_gradz_power_da_dt_discrete"])[1],
    )
    if len(full_da_dt_xz) != len(residual):
        raise SystemExit("Series length mismatch for full_da_dt_xz")

    candidates = {
        "full_da_dt_xz_current": full_da_dt_xz,
        "full_da_dt_xz_prev": shift_prev(full_da_dt_xz),
        "full_da_dt_xz_next": shift_next(full_da_dt_xz),
        "full_da_dt_xz_mid_prev": midpoint_prev(full_da_dt_xz),
        "full_da_dt_xz_mid_next": midpoint_next(full_da_dt_xz),
    }

    rows = []
    for window, t_lo, t_hi, include_right in window_bounds(t_r, fractions):
        if args.focus_window != "all" and window != args.focus_window:
            continue
        residual_window = [
            v for t, v in zip(t_r, residual)
            if t >= t_lo and (t <= t_hi if include_right else t < t_hi) and math.isfinite(v)
        ]
        if not residual_window:
            continue
        res_max, res_rms, res_final = finite_stats(residual_window)
        rows.append({
            "case_name": case_dir.name,
            "window": window,
            "candidate": "em_bulk_ledger_residual",
            "samples": len(residual_window),
            "candidate_rms_abs_rate": res_rms,
            "candidate_max_abs_rate": res_max,
            "corr_with_residual": 1.0,
            "fit_coeff": 1.0,
            "corrected_rms_abs_rate": 0.0,
            "corrected_max_abs_rate": 0.0,
            "corrected_final_rate": 0.0,
            "rms_reduction_factor": math.inf,
            "candidate_share_of_residual_rms": 1.0,
        })
        for name, series in candidates.items():
            _, y_win, x_win = paired_window(t_r, residual, series, t_lo, t_hi, include_right)
            if len(y_win) < 2:
                continue
            x_max, x_rms, x_final = finite_stats(x_win)
            alpha = fit_single(y_win, x_win)
            corrected = [y - alpha * x for y, x in zip(y_win, x_win)] if math.isfinite(alpha) else []
            c_max, c_rms, c_final = finite_stats(corrected)
            rows.append({
                "case_name": case_dir.name,
                "window": window,
                "candidate": name,
                "samples": len(y_win),
                "candidate_rms_abs_rate": x_rms,
                "candidate_max_abs_rate": x_max,
                "corr_with_residual": pearson(x_win, y_win),
                "fit_coeff": alpha,
                "corrected_rms_abs_rate": c_rms,
                "corrected_max_abs_rate": c_max,
                "corrected_final_rate": c_final,
                "rms_reduction_factor": safe_ratio(res_rms, c_rms),
                "candidate_share_of_residual_rms": safe_ratio(x_rms, res_rms),
            })

    if not rows:
        raise SystemExit("No Step-184 rows were computed")

    rows.sort(key=lambda r: (r["window"], -(r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf), r["candidate"]))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-184 A_w Time-Centering Probe\n\n")
        f.write(f"Baseline case: `{case_dir.name}`.\n\n")
        f.write(f"Controlled psi0 span: `{safe_float(controlled.get('final_psi0_span', 'nan')):.6e}`\n")
        f.write(f"Full psi0 span: `{safe_float(full.get('final_psi0_span', 'nan')):.6e}`\n\n")
        windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
        for window in windows:
            grp = [r for r in rows if r["window"] == window]
            if not grp:
                continue
            target = next(r for r in grp if r["candidate"] == "em_bulk_ledger_residual")
            ranked = [r for r in grp if r["candidate"] != "em_bulk_ledger_residual"]
            f.write(f"## window={window}\n\n")
            f.write(f"- residual_rms_abs_rate: `{fmt(target['candidate_rms_abs_rate'])}`\n")
            f.write(f"- residual_max_abs_rate: `{fmt(target['candidate_max_abs_rate'])}`\n\n")
            f.write("| candidate | samples | corr_with_residual | fit_coeff | candidate_share_of_residual_rms | corrected_rms_abs_rate | rms_reduction_factor |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|\n")
            for r in ranked:
                f.write(
                    f"| {r['candidate']} | {r['samples']} | {fmt(r['corr_with_residual'])} | {fmt(r['fit_coeff'])} | "
                    f"{fmt(r['candidate_share_of_residual_rms'])} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
