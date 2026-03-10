#!/usr/bin/env python3
"""Compare EM bulk ledger against discrete bulk-energy replacements for the dominant A_w channel."""

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


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    return f"{float(v):.6e}"


def make_window_bounds(times, fractions):
    start = times[0]
    end = times[-1]
    span = end - start
    windows = ["early", "mid", "late"]
    out = []
    for i in range(len(fractions) - 1):
        out.append((
            windows[i] if i < len(windows) else f"w{i}",
            start + fractions[i] * span,
            start + fractions[i + 1] * span,
            i == len(fractions) - 2,
        ))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True)
    parser.add_argument("--window-fractions", default="0,0.333333,0.666667,1")
    parser.add_argument("--focus-window", default="late", choices=["early", "mid", "late", "all"])
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]

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
            "m4d_mixed_cz_daw_dz2_mode_diag",
            "m4d_em_u_bulk_discrete_aw_modes12",
            "m4d_em_u_bulk_discrete_aw_modes0123",
        ]
        if not all(k in cols for k in required):
            continue

        times = cols["time"]
        t_r, du = derivative_series(times, cols["m4d_em_u_bulk"])
        _, dja = derivative_series(times, cols["m4d_int_ja_ea"])
        _, djw = derivative_series(times, cols["m4d_int_jw_ew"])
        _, ddiag = derivative_series(times, cols["m4d_mixed_cz_daw_dz2_mode_diag"])
        _, ddisc12 = derivative_series(times, cols["m4d_em_u_bulk_discrete_aw_modes12"])
        _, ddisc0123 = derivative_series(times, cols["m4d_em_u_bulk_discrete_aw_modes0123"])
        if not t_r:
            continue

        residual = [u + a + w for u, a, w in zip(du, dja, djw)]
        residual_corr12 = [u - d0 + d1 + a + w for u, d0, d1, a, w in zip(du, ddiag, ddisc12, dja, djw)]
        residual_corr0123 = [u - d0 + d1 + a + w for u, d0, d1, a, w in zip(du, ddiag, ddisc0123, dja, djw)]
        series = {
            "baseline_residual": residual,
            "diag_channel_rate": ddiag,
            "discrete_aw_modes12_rate": ddisc12,
            "discrete_aw_modes0123_rate": ddisc0123,
            "corrected_residual_modes12": residual_corr12,
            "corrected_residual_modes0123": residual_corr0123,
        }

        for window, t_lo, t_hi, include_right in make_window_bounds(t_r, fractions):
            if args.focus_window != "all" and window != args.focus_window:
                continue
            base = select_window(t_r, residual, t_lo, t_hi, include_right)
            if not base:
                continue
            base_max, base_rms, base_final = finite_stats(base)
            for name, vals in series.items():
                window_vals = select_window(t_r, vals, t_lo, t_hi, include_right)
                max_abs, rms_abs, final_rate = finite_stats(window_vals)
                rows.append({
                    "case_name": case_dir.name,
                    "window": window,
                    "component": name,
                    "max_abs_rate": max_abs,
                    "rms_abs_rate": rms_abs,
                    "final_rate": final_rate,
                    "share_of_baseline_rms": rms_abs / base_rms if base_rms > 0.0 else math.nan,
                })

    if not rows:
        raise SystemExit("No Step-176 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    focus_windows = ["early", "mid", "late"] if args.focus_window == "all" else [args.focus_window]
    order = {"early": 0, "mid": 1, "late": 2}
    grouped = {}
    for row in rows:
        if row["window"] in focus_windows:
            grouped.setdefault((row["case_name"], row["window"]), []).append(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-176 Discrete Bulk-Energy Equivalence Probe\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for key in sorted(grouped.keys(), key=lambda x: (x[0], order.get(x[1], 99))):
            case_name, window = key
            grp = grouped[key]
            f.write(f"## {case_name} window={window}\n\n")
            f.write("| component | rms_abs_rate | share_of_baseline_rms | max_abs_rate | final_rate |\n")
            f.write("|:---|---:|---:|---:|---:|\n")
            for row in grp:
                f.write(
                    f"| {row['component']} | {fmt(row['rms_abs_rate'])} | "
                    f"{fmt(row['share_of_baseline_rms'])} | {fmt(row['max_abs_rate'])} | "
                    f"{fmt(row['final_rate'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
