#!/usr/bin/env python3
"""Compare EM bulk ledger residuals under alternative source-work time centerings."""

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
        v
        for t, v in zip(times, values)
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

    combos = {
        "legacy_oldnew": ("m4d_int_ja_ea", "m4d_int_jw_ew"),
        "old_old": ("m4d_int_ja_ea_old", "m4d_int_jw_ew_old"),
        "old_mid": ("m4d_int_ja_ea_old", "m4d_int_jw_ew_mid"),
        "mid_mid": ("m4d_int_ja_ea_mid", "m4d_int_jw_ew_mid"),
        "mid_new": ("m4d_int_ja_ea_mid", "m4d_int_jw_ew_new"),
        "new_new": ("m4d_int_ja_ea_new", "m4d_int_jw_ew_new"),
    }

    for case_dir in case_dirs:
        hst = case_dir / "harris_full.out1.hst"
        if not hst.exists():
            continue
        cols = hsm.parse_hst(hst)
        required = ["time", "m4d_em_u_bulk", *{c for pair in combos.values() for c in pair}]
        if not all(k in cols for k in required):
            continue

        times = cols["time"]
        t_r, du = derivative_series(times, cols["m4d_em_u_bulk"])
        if not t_r:
            continue

        start, end = t_r[0], t_r[-1]
        t_lo = start + (end - start) * 2.0 / 3.0

        for combo_name, (ja_name, jw_name) in combos.items():
            _, dja = derivative_series(times, cols[ja_name])
            _, djw = derivative_series(times, cols[jw_name])
            residual = [u + a + w for u, a, w in zip(du, dja, djw)]
            window = select_window(t_r, residual, t_lo, end, include_right=True)
            max_abs, rms_abs, final_rate = finite_stats(window)
            rows.append(
                {
                    "case_name": case_dir.name,
                    "candidate": combo_name,
                    "ja_series": ja_name,
                    "jw_series": jw_name,
                    "max_abs_rate": max_abs,
                    "rms_abs_rate": rms_abs,
                    "final_rate": final_rate,
                }
            )

    if not rows:
        raise SystemExit("No Step-189 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-189 Source-Work Time-Centering Probe\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            f.write(f"## {case_name} window=late\n\n")
            ranked = [r for r in rows if r["case_name"] == case_name]
            ranked.sort(key=lambda r: (math.inf if math.isnan(r["rms_abs_rate"]) else r["rms_abs_rate"]))
            f.write("| candidate | ja_series | jw_series | rms_abs_rate | max_abs_rate | final_rate |\n")
            f.write("|:---|:---|:---|---:|---:|---:|\n")
            for row in ranked:
                f.write(
                    f"| {row['candidate']} | {row['ja_series']} | {row['jw_series']} | "
                    f"{fmt(row['rms_abs_rate'])} | {fmt(row['max_abs_rate'])} | "
                    f"{fmt(row['final_rate'])} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
