#!/usr/bin/env python3
"""Test ledger-residual corrections built from A_w gradient-power diagnostics."""

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
    return (out_t, out_v)


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


def select_window(times, values, t_lo, t_hi, include_right=False):
    return [v for t, v in zip(times, values) if t >= t_lo and (t <= t_hi if include_right else t < t_hi)]


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    return f"{float(v):.6e}"


def parse_case_rows(run_log: Path):
    lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()
    header = controlled = full = None
    for i, line in enumerate(lines):
        if not line.startswith("case,final_psi0_span,"):
            continue
        header = line
        for j in range(i + 1, min(i + 24, len(lines))):
            if lines[j].startswith("controlled,"):
                controlled = lines[j]
            elif lines[j].startswith("full,"):
                full = lines[j]
        if header and controlled and full:
            break
    if header is None or controlled is None or full is None:
        raise RuntimeError(f"Could not find controlled/full rows in {run_log}")
    keys = next(csv.reader([header]))
    controlled_vals = next(csv.reader([controlled]))
    full_vals = next(csv.reader([full]))
    return dict(zip(keys, controlled_vals)), dict(zip(keys, full_vals))


def safe_float(v: str) -> float:
    try:
        return float(v)
    except Exception:
        return math.nan


def safe_ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def window_bounds(times, fractions):
    start = times[0]
    end = times[-1]
    span = end - start
    if span <= 0.0:
        return []
    windows = ["early", "mid", "late"]
    bounds = []
    for idx in range(len(fractions) - 1):
        f_lo = fractions[idx]
        f_hi = fractions[idx + 1]
        name = windows[idx] if idx < len(windows) else f"w{idx}"
        t_lo = start + f_lo * span
        t_hi = start + f_hi * span
        include_right = idx == len(fractions) - 2
        bounds.append((name, t_lo, t_hi, include_right))
    return bounds


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
    if not t_r or len(du) != len(dja) or len(du) != len(djw):
        raise SystemExit("Could not construct ledger residual series")
    residual = [u + a + w for u, a, w in zip(du, dja, djw)]

    def ds(name: str):
        return derivative_series(times, cols.get(name))[1]

    candidates = {
        "gradz_power_da_dt_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_da_dt"), ds("m4d_int_aw_mode2_gradz_power_da_dt"))],
        "gradz_power_mix_a_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_mix_a"), ds("m4d_int_aw_mode2_gradz_power_mix_a"))],
        "gradz_power_mix_w0_local_gradient_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_mix_w0_local_gradient"), ds("m4d_int_aw_mode2_gradz_power_mix_w0_local_gradient"))],
        "gradz_power_mix_w0_local_dz_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_mix_w0_local_dz"), ds("m4d_int_aw_mode2_gradz_power_mix_w0_local_dz"))],
        "gradz_power_mix_w0_local_mode2_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_mix_w0_local_mode2"), ds("m4d_int_aw_mode2_gradz_power_mix_w0_local_mode2"))],
        "gradz_power_mix_w0_local_dz_signed_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_mix_w0_local_dz_signed"), ds("m4d_int_aw_mode2_gradz_power_mix_w0_local_dz_signed"))],
        "gradz_power_mix_w0_local_dz_signed_mode2_modes12": [a + b for a, b in zip(ds("m4d_int_aw_mode1_gradz_power_mix_w0_local_dz_signed_mode2"), ds("m4d_int_aw_mode2_gradz_power_mix_w0_local_dz_signed_mode2"))],
        "mixed_cz_daw_dz2_modes12": [a + b for a, b in zip(ds("m4d_mixed_cz_daw_dz2_mode1"), ds("m4d_mixed_cz_daw_dz2_mode2"))],
        "mixed_cz_daw_dz2_mode_diag": ds("m4d_mixed_cz_daw_dz2_mode_diag"),
    }

    if any(len(v) != len(residual) for v in candidates.values()):
        raise SystemExit("Candidate series length mismatch")

    rows = []
    for window, t_lo, t_hi, include_right in window_bounds(t_r, fractions):
        if args.focus_window != "all" and window != args.focus_window:
            continue
        r_win = select_window(t_r, residual, t_lo, t_hi, include_right)
        if not r_win:
            continue
        r_max, r_rms, r_final = finite_stats(r_win)
        for name, series in candidates.items():
            s_win = select_window(t_r, series, t_lo, t_hi, include_right)
            s_max, s_rms, s_final = finite_stats(s_win)
            minus = [rv - sv for rv, sv in zip(r_win, s_win)]
            plus = [rv + sv for rv, sv in zip(r_win, s_win)]
            m_max, m_rms, m_final = finite_stats(minus)
            p_max, p_rms, p_final = finite_stats(plus)
            best_sign = "-" if (math.isfinite(m_rms) and (not math.isfinite(p_rms) or m_rms <= p_rms)) else "+"
            best_rms = m_rms if best_sign == "-" else p_rms
            best_max = m_max if best_sign == "-" else p_max
            best_final = m_final if best_sign == "-" else p_final
            rows.append({
                "case_name": case_dir.name,
                "window": window,
                "candidate": name,
                "residual_rms_abs_rate": r_rms,
                "residual_max_abs_rate": r_max,
                "candidate_rms_abs_rate": s_rms,
                "candidate_max_abs_rate": s_max,
                "corr_with_residual": pearson(s_win, r_win),
                "best_sign": best_sign,
                "corrected_rms_abs_rate": best_rms,
                "corrected_max_abs_rate": best_max,
                "corrected_final_rate": best_final,
                "rms_reduction_factor": safe_ratio(r_rms, best_rms),
                "candidate_share_of_residual_rms": safe_ratio(s_rms, r_rms),
            })

    if not rows:
        raise SystemExit("No correction rows were computed")

    rows.sort(key=lambda r: (r["window"], -(r["rms_reduction_factor"] if math.isfinite(r["rms_reduction_factor"]) else -math.inf), r["candidate"]))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-171 A_w Gradient-Power Correction Probe\n\n")
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
            f.write("| candidate | corr_with_residual | candidate_share_of_residual_rms | best_sign | corrected_rms_abs_rate | rms_reduction_factor | corrected_max_abs_rate | corrected_final_rate |\n")
            f.write("|:---|---:|---:|:---:|---:|---:|---:|---:|\n")
            for r in grp:
                f.write(
                    f"| {r['candidate']} | {fmt(r['corr_with_residual'])} | {fmt(r['candidate_share_of_residual_rms'])} | {r['best_sign']} | {fmt(r['corrected_rms_abs_rate'])} | {fmt(r['rms_reduction_factor'])} | {fmt(r['corrected_max_abs_rate'])} | {fmt(r['corrected_final_rate'])} |\n"
                )
            f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
