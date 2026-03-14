#!/usr/bin/env python3
"""Compute corrected EM bulk-ledger metrics using fixed post-fix z-block rules."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path

Z_COEFF_FIXED = 8.304215e-01
X_COEFF_FIXED = 8.288144e-02
LEDGER_TOL = 1.0


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


def status_from_max(max_abs):
    if not math.isfinite(max_abs):
        return "N/A"
    return "PASS" if max_abs <= LEDGER_TOL else "FAIL"


def sum_series(*series_list):
    return [sum(vals) for vals in zip(*series_list)]


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
            "m4d_int_cx_energy_delta_exact_discrete_sum",
            "m4d_int_cz_energy_delta_exact_discrete_sum",
        ] + [f"m4d_int_aw_mode{n}_gradx_power_da_dt_discrete" for n in range(4)] + [
            f"m4d_int_aw_mode{n}_gradz_power_da_dt_discrete" for n in range(4)
        ]
        if not all(k in cols for k in required):
            continue

        times = cols["time"]
        em = cols["m4d_em_u_bulk"]
        ja = cols["m4d_int_ja_ea"]
        jw = cols["m4d_int_jw_ew"]
        cx = cols["m4d_int_cx_energy_delta_exact_discrete_sum"]
        cz = cols["m4d_int_cz_energy_delta_exact_discrete_sum"]
        x_da = [cols[f"m4d_int_aw_mode{n}_gradx_power_da_dt_discrete"] for n in range(4)]
        z_da = [cols[f"m4d_int_aw_mode{n}_gradz_power_da_dt_discrete"] for n in range(4)]
        x_block = sum_series(cx, *x_da)
        z_block = sum_series(cz, *z_da)

        variants = {
            "legacy_bulk": em,
            "z_block_fixed": [e - Z_COEFF_FIXED * z for e, z in zip(em, z_block)],
            "xz_block_fixed": [
                e - Z_COEFF_FIXED * z - X_COEFF_FIXED * x
                for e, z, x in zip(em, z_block, x_block)
            ],
        }
        for label, em_variant in variants.items():
            max_abs, rms_abs, final_rate = hsm.em_bulk_ledger_abs_rate_metrics(
                times, em_variant, ja, jw
            )
            rows.append(
                {
                    "case_name": case_dir.name,
                    "candidate": label,
                    "ledger_max_abs_rate": max_abs,
                    "ledger_rms_abs_rate": rms_abs,
                    "ledger_final_rate": final_rate,
                    "ledger_status": status_from_max(max_abs),
                    "z_coeff_fixed": Z_COEFF_FIXED if label != "legacy_bulk" else 0.0,
                    "x_coeff_fixed": X_COEFF_FIXED if label == "xz_block_fixed" else 0.0,
                }
            )

    if not rows:
        raise SystemExit("No Step-202 rows were computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-202 Z-Block Corrected Ledger A/B\n\n")
        f.write(f"- ledger_abs_rate_tol: `{fmt(LEDGER_TOL)}`\n")
        f.write(f"- fixed_z_coeff: `{fmt(Z_COEFF_FIXED)}`\n")
        f.write(f"- fixed_x_coeff: `{fmt(X_COEFF_FIXED)}`\n\n")
        for case_name in sorted({r["case_name"] for r in rows}):
            f.write(f"## {case_name}\n\n")
            f.write("| candidate | max abs rate | rms abs rate | final rate | status |\n")
            f.write("|:---|---:|---:|---:|:---|\n")
            ranked = [r for r in rows if r["case_name"] == case_name]
            ranked.sort(key=lambda r: (r["ledger_status"] != "PASS", r["ledger_max_abs_rate"]))
            for r in ranked:
                f.write(
                    f"| {r['candidate']} | {fmt(r['ledger_max_abs_rate'])} | {fmt(r['ledger_rms_abs_rate'])} | "
                    f"{fmt(r['ledger_final_rate'])} | {r['ledger_status']} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
