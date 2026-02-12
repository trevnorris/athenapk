#!/usr/bin/env python3
"""Decompose late-window EM ledger residual into signed history-term rates."""

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_harris_scan_matrix_module(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_gain_tag(tag: str):
    sign = 1.0
    if tag.startswith("m"):
        sign = -1.0
        tag = tag[1:]
    elif tag.startswith("minus"):
        sign = -1.0
        tag = tag[5:]

    if "em" in tag:
        base, exp = tag.split("em", 1)
        if not base or not exp:
            return math.nan
        try:
            return sign * float(base) * (10.0 ** (-int(exp)))
        except ValueError:
            return math.nan
    try:
        return sign * float(tag)
    except ValueError:
        return math.nan


def parse_window_fractions(raw: str):
    vals = [float(x.strip()) for x in raw.split(",") if x.strip()]
    if len(vals) < 2:
        raise ValueError("Need at least two fractions.")
    if abs(vals[0]) > 1.0e-12 or abs(vals[-1] - 1.0) > 1.0e-12:
        raise ValueError("Fractions must start at 0 and end at 1.")
    for i in range(1, len(vals)):
        if vals[i] <= vals[i - 1]:
            raise ValueError("Fractions must be strictly increasing.")
    return vals


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


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


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    max_abs = max(abs(v) for v in finite)
    rms_abs = math.sqrt(sum(v * v for v in finite) / len(finite))
    final = finite[-1]
    return (max_abs, rms_abs, final)


def derivative_series(times, values):
    if values is None:
        return ([], [])
    if len(times) < 2 or len(values) != len(times):
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


def select_window(times, values, t_lo, t_hi, include_left=True, include_right=False):
    out = []
    for t, v in zip(times, values):
        ge = (t >= t_lo) if include_left else (t > t_lo)
        le = (t <= t_hi) if include_right else (t < t_hi)
        if ge and le:
            out.append(v)
    return out


def compute_term_rows(case_dir, case_name, cols, gain_tag, gain, fractions):
    times = cols.get("time")
    if times is None or len(times) < 2:
        return []

    # Primary EM bulk-ledger residual terms.
    _, du = derivative_series(times, cols.get("m4d_em_u_bulk"))
    t_r, dja = derivative_series(times, cols.get("m4d_int_ja_ea"))
    _, djw = derivative_series(times, cols.get("m4d_int_jw_ew"))
    if not t_r or len(du) != len(dja) or len(du) != len(djw):
        return []
    residual = [u + a + w for u, a, w in zip(du, dja, djw)]

    # Candidate signed channels for decomposition.
    term_columns = {
        "d_em_u_bulk_dt": "m4d_em_u_bulk",
        "d_int_ja_ea_dt": "m4d_int_ja_ea",
        "d_int_jw_ew_dt": "m4d_int_jw_ew",
        "d_int_bridge_power_mode0_dt": "m4d_int_bridge_power_mode0",
        "d_int_bridge_power_mode2_dt": "m4d_int_bridge_power_mode2",
        "d_int_bridge_power_sum_dt": "m4d_int_bridge_power_sum",
        "d_int_src_timelike_a0_from_piw_dt": "m4d_int_src_timelike_a0_from_piw",
        "d_int_src_timelike_aw_from_pi0_dt": "m4d_int_src_timelike_aw_from_pi0",
        "d_int_src_em_gauge_dt": "m4d_int_src_em_gauge",
        "d_int_rhs_a0_mode0_lap_dt": "m4d_int_rhs_a0_mode0_lap",
        "d_int_rhs_a0_mode0_mass_dt": "m4d_int_rhs_a0_mode0_mass",
        "d_int_rhs_a0_mode0_current_dt": "m4d_int_rhs_a0_mode0_current",
        "d_int_rhs_a0_mode0_damping_dt": "m4d_int_rhs_a0_mode0_damping",
        "d_int_rhs_a0_mode0_timelike_dt": "m4d_int_rhs_a0_mode0_timelike",
        "d_int_rhs_a0_mode0_gauge_dt": "m4d_int_rhs_a0_mode0_gauge",
        "d_int_rhs_a0_mode0_total_dt": "m4d_int_rhs_a0_mode0_total",
        "d_int_rhs_ay_mode0_lap_dt": "m4d_int_rhs_ay_mode0_lap",
        "d_int_rhs_ay_mode0_mass_dt": "m4d_int_rhs_ay_mode0_mass",
        "d_int_rhs_ay_mode0_current_dt": "m4d_int_rhs_ay_mode0_current",
        "d_int_rhs_ay_mode0_damping_dt": "m4d_int_rhs_ay_mode0_damping",
        "d_int_rhs_ay_mode0_spatial_mixed_dt": "m4d_int_rhs_ay_mode0_spatial_mixed",
        "d_int_rhs_ay_mode0_even_bridge_dt": "m4d_int_rhs_ay_mode0_even_bridge",
        "d_int_rhs_ay_mode0_total_dt": "m4d_int_rhs_ay_mode0_total",
        "d_int_rhs_aw_mode0_lap_dt": "m4d_int_rhs_aw_mode0_lap",
        "d_int_rhs_aw_mode0_current_dt": "m4d_int_rhs_aw_mode0_current",
        "d_int_rhs_aw_mode0_damping_dt": "m4d_int_rhs_aw_mode0_damping",
        "d_int_rhs_aw_mode0_spatial_mixed_dt": "m4d_int_rhs_aw_mode0_spatial_mixed",
        "d_int_rhs_aw_mode0_timelike_dt": "m4d_int_rhs_aw_mode0_timelike",
        "d_int_rhs_aw_mode0_total_dt": "m4d_int_rhs_aw_mode0_total",
        "d_int_rhs_aw_mode2_even_bridge_dt": "m4d_int_rhs_aw_mode2_even_bridge",
    }

    term_rate_series = {}
    for term, col in term_columns.items():
        t_term, v_term = derivative_series(times, cols.get(col))
        if not t_term:
            continue
        if len(t_term) != len(t_r):
            continue
        term_rate_series[term] = v_term

    start = t_r[0]
    end = t_r[-1]
    span = end - start
    if span <= 0.0:
        return []

    rows = []
    for idx in range(len(fractions) - 1):
        f_lo = fractions[idx]
        f_hi = fractions[idx + 1]
        window = "early" if idx == 0 else ("late" if idx == len(fractions) - 2 else "mid")
        t_lo = start + f_lo * span
        t_hi = start + f_hi * span
        include_right = idx == len(fractions) - 2

        residual_window = select_window(t_r, residual, t_lo, t_hi, True, include_right)
        res_max, res_rms, res_final = finite_stats(residual_window)
        if not residual_window:
            continue

        # Baseline summary row.
        rows.append(
            {
                "case_dir": str(case_dir),
                "gain_tag": gain_tag,
                "gain": gain,
                "case": case_name,
                "window": window,
                "term": "residual_em_bulk_baseline",
                "max_abs_rate": res_max,
                "rms_abs_rate": res_rms,
                "final_rate": res_final,
                "corr_with_residual": 1.0,
                "share_of_residual_rms": 1.0,
                "residual_max_abs_rate": res_max,
                "residual_rms_abs_rate": res_rms,
            }
        )

        # Decomposition terms.
        for term_name, term_series in term_rate_series.items():
            term_window = select_window(t_r, term_series, t_lo, t_hi, True, include_right)
            term_max, term_rms, term_final = finite_stats(term_window)
            corr = pearson(term_window, residual_window)
            share = (
                term_rms / res_rms
                if math.isfinite(term_rms) and math.isfinite(res_rms) and res_rms > 0.0
                else math.nan
            )
            rows.append(
                {
                    "case_dir": str(case_dir),
                    "gain_tag": gain_tag,
                    "gain": gain,
                    "case": case_name,
                    "window": window,
                    "term": term_name,
                    "max_abs_rate": term_max,
                    "rms_abs_rate": term_rms,
                    "final_rate": term_final,
                    "corr_with_residual": corr,
                    "share_of_residual_rms": share,
                    "residual_max_abs_rate": res_max,
                    "residual_rms_abs_rate": res_rms,
                }
            )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        default="/projects/fluid-engine/out/step69_paired_window_g*",
        help="Glob for case directories containing controlled/full .hst outputs",
    )
    parser.add_argument(
        "--window-fractions",
        default="0,0.333333,0.666667,1",
        help="Comma-separated monotonically increasing fractional boundaries from 0 to 1",
    )
    parser.add_argument(
        "--focus-window",
        default="late",
        choices=["early", "mid", "late", "all"],
        help="Window to emphasize in markdown ranking tables",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of top terms (by RMS share) to print per focus group",
    )
    parser.add_argument(
        "--out-csv",
        default="/projects/fluid-engine/out/step70_ledger_term_decomp.csv",
        help="Path to write decomposition CSV",
    )
    parser.add_argument(
        "--out-md",
        default="/projects/fluid-engine/out/step70_ledger_term_decomp.md",
        help="Path to write markdown summary",
    )
    args = parser.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
    if not dirs:
        raise SystemExit(f"No directories matched: {args.glob}")

    rows = []
    for d in dirs:
        controlled_hst = d / "harris_controlled.out1.hst"
        full_hst = d / "harris_full.out1.hst"
        if not controlled_hst.exists() or not full_hst.exists():
            continue

        gain_tag = d.name.split("_g", 1)[-1]
        gain = parse_gain_tag(gain_tag)
        controlled_cols = hsm.parse_hst(controlled_hst)
        full_cols = hsm.parse_hst(full_hst)

        rows.extend(compute_term_rows(d, "controlled", controlled_cols, gain_tag, gain, fractions))
        rows.extend(compute_term_rows(d, "full", full_cols, gain_tag, gain, fractions))

    if not rows:
        raise SystemExit("No decomposition rows were computed.")

    rows.sort(
        key=lambda r: (
            (math.inf if not math.isfinite(r["gain"]) else r["gain"]),
            r["case_dir"],
            0 if r["case"] == "controlled" else 1,
            {"early": 0, "mid": 1, "late": 2}.get(r["window"], 9),
            r["term"],
        )
    )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Ledger Term Decomposition\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")

        focus_windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
        focus_rows = [r for r in rows if r["case"] == "full" and r["window"] in focus_windows]
        if not focus_rows:
            f.write("No focus rows found.\n")
        else:
            grouped = {}
            for r in focus_rows:
                grouped.setdefault((r["gain"], r["window"]), []).append(r)

            for (gain, window) in sorted(grouped.keys(), key=lambda x: (x[0], x[1])):
                grp = grouped[(gain, window)]
                residual = [r for r in grp if r["term"] == "residual_em_bulk_baseline"]
                residual_rms = residual[0]["rms_abs_rate"] if residual else math.nan
                residual_max = residual[0]["max_abs_rate"] if residual else math.nan
                f.write(f"## gain={fmt(gain)} window={window}\n\n")
                f.write(
                    f"- residual_rms_abs_rate: `{fmt(residual_rms)}`\n"
                    f"- residual_max_abs_rate: `{fmt(residual_max)}`\n\n"
                )
                ranked = [
                    r
                    for r in grp
                    if r["term"] != "residual_em_bulk_baseline"
                    and math.isfinite(r["share_of_residual_rms"])
                ]
                ranked.sort(key=lambda r: abs(r["share_of_residual_rms"]), reverse=True)
                ranked = ranked[: max(args.top_n, 1)]
                f.write("| term | rms_abs_rate | share_of_residual_rms | corr_with_residual | max_abs_rate | final_rate |\n")
                f.write("|:---|---:|---:|---:|---:|---:|\n")
                for r in ranked:
                    f.write(
                        "| "
                        + " | ".join(
                            [
                                r["term"],
                                fmt(r["rms_abs_rate"]),
                                fmt(r["share_of_residual_rms"]),
                                fmt(r["corr_with_residual"]),
                                fmt(r["max_abs_rate"]),
                                fmt(r["final_rate"]),
                            ]
                        )
                        + " |\n"
                    )
                f.write("\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
