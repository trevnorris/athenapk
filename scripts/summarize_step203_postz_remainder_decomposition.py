#!/usr/bin/env python3
"""Decompose the small post-z-block ledger remainder into surviving energy channels."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import itertools
import math
from pathlib import Path

Z_COEFF_FIXED = 8.304215e-01


def load_hsm(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    if v == math.inf:
        return "inf"
    return f"{float(v):.6e}"


def derivative_skip_initial_interval(times, values):
    if values is None or len(times) < 3 or len(values) != len(times):
        return ([], [])
    out_t, out_v = [], []
    for i in range(1, len(times)):
        if i == 1:
            continue
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


def fit_single(y, x):
    xx = sum(v * v for v in x)
    if xx <= 0.0:
        return math.nan
    return sum(a * b for a, b in zip(y, x)) / xx


def fit_pair(y, x0, x1):
    s00 = sum(v * v for v in x0)
    s11 = sum(v * v for v in x1)
    s01 = sum(a * b for a, b in zip(x0, x1))
    b0 = sum(a * b for a, b in zip(y, x0))
    b1 = sum(a * b for a, b in zip(y, x1))
    det = s00 * s11 - s01 * s01
    if abs(det) <= 0.0:
        return (math.nan, math.nan)
    c0 = (b0 * s11 - b1 * s01) / det
    c1 = (b1 * s00 - b0 * s01) / det
    return (c0, c1)


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
            "m4d_int_cz_energy_delta_exact_discrete_sum",
            "m4d_mixed_cx2",
            "m4d_mixed_cy2",
            "m4d_mixed_cz2",
            "m4d_mixed_c2",
            "m4d_mixed_ew2",
            "m4d_em_u_resolved",
            "m4d_int_cx_energy_delta_exact_discrete_sum",
        ] + [f"m4d_int_aw_mode{n}_gradz_power_da_dt_discrete" for n in range(4)] + [
            f"m4d_int_aw_mode{n}_gradx_power_da_dt_discrete" for n in range(4)
        ]
        if not all(k in cols for k in required):
            continue

        times = cols["time"]
        z_block = [
            cols["m4d_int_cz_energy_delta_exact_discrete_sum"][i]
            + sum(cols[f"m4d_int_aw_mode{n}_gradz_power_da_dt_discrete"][i] for n in range(4))
            for i in range(len(times))
        ]
        corrected_bulk = [
            cols["m4d_em_u_bulk"][i] - Z_COEFF_FIXED * z_block[i] for i in range(len(times))
        ]
        t_r, corrected_residual = derivative_skip_initial_interval(
            times,
            [
                corrected_bulk[i] + cols["m4d_int_ja_ea"][i] + cols["m4d_int_jw_ew"][i]
                for i in range(len(times))
            ],
        )
        if not t_r:
            continue

        x_block_state = [
            cols["m4d_int_cx_energy_delta_exact_discrete_sum"][i]
            + sum(cols[f"m4d_int_aw_mode{n}_gradx_power_da_dt_discrete"][i] for n in range(4))
            for i in range(len(times))
        ]
        candidate_state = {
            "x_block_sum": x_block_state,
            "mixed_cx2": cols["m4d_mixed_cx2"],
            "mixed_cy2": cols["m4d_mixed_cy2"],
            "mixed_cz2": cols["m4d_mixed_cz2"],
            "mixed_c2": cols["m4d_mixed_c2"],
            "mixed_ew2": cols["m4d_mixed_ew2"],
            "resolved": cols["m4d_em_u_resolved"],
        }
        candidate_rates = {
            name: derivative_skip_initial_interval(times, values)[1]
            for name, values in candidate_state.items()
        }

        postz_max, postz_rms, postz_final = finite_stats(corrected_residual)
        start, end = t_r[0], t_r[-1]
        span = end - start
        windows = {
            "overall": (start, end),
            "early": (start, start + span / 3.0),
            "mid": (start + span / 3.0, start + 2.0 * span / 3.0),
            "late": (start + 2.0 * span / 3.0, end),
        }

        for window_name, (t_lo, t_hi) in windows.items():
            y = select_window(t_r, corrected_residual, t_lo, t_hi, include_right=True if window_name in ("overall", "late") else False)
            if not y:
                continue
            residual_max, residual_rms, residual_final = finite_stats(y)
            single_rows = []
            for name, series in candidate_rates.items():
                x = select_window(
                    t_r, series, t_lo, t_hi, include_right=True if window_name in ("overall", "late") else False
                )
                if len(x) != len(y):
                    continue
                coeff = fit_single(y, x)
                corrected = [yy - coeff * xx for yy, xx in zip(y, x)]
                _, corrected_rms, _ = finite_stats(corrected)
                reduction = residual_rms / corrected_rms if corrected_rms > 0.0 else math.inf
                candidate_max, candidate_rms, candidate_final = finite_stats(x)
                row = {
                    "case_name": case_dir.name,
                    "window": window_name,
                    "fit_kind": "single",
                    "candidate_0": name,
                    "candidate_1": "",
                    "coeff_0": coeff,
                    "coeff_1": math.nan,
                    "corr_with_residual": pearson(x, y),
                    "candidate_0_rms_abs_rate": candidate_rms,
                    "candidate_1_rms_abs_rate": math.nan,
                    "residual_max_abs_rate": residual_max,
                    "residual_rms_abs_rate": residual_rms,
                    "residual_final_rate": residual_final,
                    "corrected_rms_abs_rate": corrected_rms,
                    "rms_reduction_factor": reduction,
                    "postz_max_abs_rate": postz_max,
                    "postz_rms_abs_rate": postz_rms,
                    "postz_final_rate": postz_final,
                    "z_coeff_fixed": Z_COEFF_FIXED,
                }
                rows.append(row)
                single_rows.append(row)

            best_pair_row = None
            for name0, name1 in itertools.combinations(candidate_rates.keys(), 2):
                x0 = select_window(
                    t_r,
                    candidate_rates[name0],
                    t_lo,
                    t_hi,
                    include_right=True if window_name in ("overall", "late") else False,
                )
                x1 = select_window(
                    t_r,
                    candidate_rates[name1],
                    t_lo,
                    t_hi,
                    include_right=True if window_name in ("overall", "late") else False,
                )
                if len(x0) != len(y) or len(x1) != len(y):
                    continue
                c0, c1 = fit_pair(y, x0, x1)
                if not math.isfinite(c0) or not math.isfinite(c1):
                    continue
                corrected = [yy - c0 * v0 - c1 * v1 for yy, v0, v1 in zip(y, x0, x1)]
                _, corrected_rms, _ = finite_stats(corrected)
                reduction = residual_rms / corrected_rms if corrected_rms > 0.0 else math.inf
                pair_row = {
                    "case_name": case_dir.name,
                    "window": window_name,
                    "fit_kind": "pair",
                    "candidate_0": name0,
                    "candidate_1": name1,
                    "coeff_0": c0,
                    "coeff_1": c1,
                    "corr_with_residual": math.nan,
                    "candidate_0_rms_abs_rate": finite_stats(x0)[1],
                    "candidate_1_rms_abs_rate": finite_stats(x1)[1],
                    "residual_max_abs_rate": residual_max,
                    "residual_rms_abs_rate": residual_rms,
                    "residual_final_rate": residual_final,
                    "corrected_rms_abs_rate": corrected_rms,
                    "rms_reduction_factor": reduction,
                    "postz_max_abs_rate": postz_max,
                    "postz_rms_abs_rate": postz_rms,
                    "postz_final_rate": postz_final,
                    "z_coeff_fixed": Z_COEFF_FIXED,
                }
                if (
                    best_pair_row is None
                    or pair_row["corrected_rms_abs_rate"] < best_pair_row["corrected_rms_abs_rate"]
                ):
                    best_pair_row = pair_row
            if best_pair_row is not None:
                rows.append(best_pair_row)

    if not rows:
        raise SystemExit("No Step-203 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-203 Post-Z Remainder Decomposition\n\n")
        f.write(f"- fixed_z_coeff: `{fmt(Z_COEFF_FIXED)}`\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            case_rows = [r for r in rows if r["case_name"] == case_name]
            postz = case_rows[0]
            f.write(f"## {case_name}\n\n")
            f.write(
                f"- post-z ledger max abs rate: `{fmt(postz['postz_max_abs_rate'])}`\n"
                f"- post-z ledger rms abs rate: `{fmt(postz['postz_rms_abs_rate'])}`\n"
                f"- post-z ledger final rate: `{fmt(postz['postz_final_rate'])}`\n\n"
            )
            for window_name in ("overall", "early", "mid", "late"):
                window_rows = [r for r in case_rows if r["window"] == window_name]
                if not window_rows:
                    continue
                f.write(f"### {window_name}\n\n")
                best_pair = next((r for r in window_rows if r["fit_kind"] == "pair"), None)
                if best_pair is not None:
                    f.write(
                        f"- best pair: `{best_pair['candidate_0']} + {best_pair['candidate_1']}` "
                        f"with corrected RMS `{fmt(best_pair['corrected_rms_abs_rate'])}` "
                        f"(coeffs `{fmt(best_pair['coeff_0'])}`, `{fmt(best_pair['coeff_1'])}`)\n\n"
                    )
                singles = [r for r in window_rows if r["fit_kind"] == "single"]
                singles.sort(key=lambda r: r["corrected_rms_abs_rate"])
                f.write("| candidate | corr | coeff | corrected rms | reduction |\n")
                f.write("|:---|---:|---:|---:|---:|\n")
                for row in singles:
                    f.write(
                        f"| {row['candidate_0']} | {fmt(row['corr_with_residual'])} | "
                        f"{fmt(row['coeff_0'])} | {fmt(row['corrected_rms_abs_rate'])} | "
                        f"{fmt(row['rms_reduction_factor'])} |\n"
                    )
                f.write("\n")


if __name__ == "__main__":
    main()
