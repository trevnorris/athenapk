#!/usr/bin/env python3
"""Summarize direct z-block ledger patch A/B using corrected vs legacy bulk channels."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path

LEDGER_TOL = 1.0


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_hsm(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fmt(v):
    try:
        return f"{float(v):.6e}"
    except Exception:
        return str(v)


def ledger_status(max_abs):
    if not math.isfinite(max_abs):
        return "N/A"
    return "PASS" if max_abs <= LEDGER_TOL else "FAIL"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--summary-csv", required=True)
    ap.add_argument("--phase2-csv", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dir = Path(args.case_dir)
    hst = case_dir / "harris_full.out1.hst"
    cols = hsm.parse_hst(hst)

    times = cols["time"]
    ja = cols["m4d_int_ja_ea"]
    jw = cols["m4d_int_jw_ew"]

    channels = {
        "legacy_bulk": cols["m4d_em_u_bulk_legacy"],
        "zblock_corrected_bulk": cols["m4d_em_u_bulk"],
        "zblock_correction_only": cols["m4d_em_u_bulk_zblock_correction"],
    }

    rows = []
    for label, em in channels.items():
        max_abs, rms_abs, final_rate = hsm.em_bulk_ledger_abs_rate_metrics(times, em, ja, jw)
        rows.append(
            {
                "candidate": label,
                "ledger_max_abs_rate": max_abs,
                "ledger_rms_abs_rate": rms_abs,
                "ledger_final_rate": final_rate,
                "ledger_status": ledger_status(max_abs),
            }
        )

    summary = {r["case_name"]: r for r in read_csv(Path(args.summary_csv))}
    phase2 = {r["case_name"]: r for r in read_csv(Path(args.phase2_csv))}
    case_name = case_dir.name
    s = summary.get(case_name, {})
    p = phase2.get(case_name, {})

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-204 Direct Z-Block Ledger Patch A/B\n\n")
        f.write(f"- case: `{case_name}`\n")
        f.write(f"- psi0 ratio full/control: `{fmt(p.get('psi0_ratio_full_over_controlled'))}`\n")
        f.write(f"- projected status: `{p.get('projected_status', 'nan')}`\n")
        f.write(f"- local status: `{p.get('local_status', 'nan')}`\n")
        f.write(f"- transport status: `{p.get('transport_status', 'nan')}`\n")
        f.write(
            f"- top transport: `{s.get('top_abs_component', 'nan')}={fmt(s.get('top_abs_rate'))}`\n\n"
        )
        f.write("| candidate | ledger max abs rate | ledger rms abs rate | final rate | status |\n")
        f.write("|:---|---:|---:|---:|:---|\n")
        rows_sorted = sorted(rows, key=lambda r: (r["ledger_status"] != "PASS", r["ledger_max_abs_rate"]))
        for row in rows_sorted:
            f.write(
                f"| {row['candidate']} | {fmt(row['ledger_max_abs_rate'])} | "
                f"{fmt(row['ledger_rms_abs_rate'])} | {fmt(row['ledger_final_rate'])} | "
                f"{row['ledger_status']} |\n"
            )


if __name__ == "__main__":
    main()
