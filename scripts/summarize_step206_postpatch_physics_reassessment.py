#!/usr/bin/env python3
"""Summarize the patched-branch physics signal after the Step-204 z-block ledger correction."""

from __future__ import annotations

import argparse
import csv
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


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    try:
        return f"{float(v):.6e}"
    except Exception:
        return str(v)


def max_abs_with_time(times, values):
    if not times or not values or len(times) != len(values):
        return (math.nan, math.nan)
    idx = max(range(len(values)), key=lambda i: abs(values[i]))
    return abs(values[idx]), times[idx]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--summary-csv", required=True)
    ap.add_argument("--phase2-csv", required=True)
    ap.add_argument("--zblock-csv", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dir = Path(args.case_dir)
    case_name = case_dir.name

    summary = {r["case_name"]: r for r in read_csv(Path(args.summary_csv))}
    phase2 = {r["case_name"]: r for r in read_csv(Path(args.phase2_csv))}
    zblock_rows = {r["candidate"]: r for r in read_csv(Path(args.zblock_csv))}

    full_hst = case_dir / "harris_full.out1.hst"
    controlled_hst = case_dir / "harris_controlled.out1.hst"
    full = hsm.parse_hst(full_hst)
    controlled = hsm.parse_hst(controlled_hst)

    full_time = full["time"]
    ctl_time = controlled["time"]

    final_full_psi0 = full["m4d_psi0_span"][-1]
    final_ctl_psi0 = controlled["m4d_psi0_span"][-1]
    final_full_proj = full["m4d_psi_proj_span"][-1]
    final_full_w0 = full["m4d_psi_w0_span"][-1]
    final_gap = final_full_w0 - final_full_proj
    max_emf, t_emf = max_abs_with_time(full_time, full["m4d_emf_vw_c_abs"])
    max_mixed_c2, t_c2 = max_abs_with_time(full_time, full["m4d_mixed_c2"])
    max_psi0, t_psi0 = max_abs_with_time(full_time, full["m4d_psi0_span"])
    max_gap, t_gap = max_abs_with_time(full_time, [a - b for a, b in zip(full["m4d_psi_w0_span"], full["m4d_psi_proj_span"])])

    legacy = zblock_rows["legacy_bulk"]
    corrected = zblock_rows["zblock_corrected_bulk"]
    legacy_max = float(legacy["ledger_max_abs_rate"])
    legacy_rms = float(legacy["ledger_rms_abs_rate"])
    corrected_max = float(corrected["ledger_max_abs_rate"])
    corrected_rms = float(corrected["ledger_rms_abs_rate"])

    row = {
        "case_name": case_name,
        "final_controlled_psi0_span": final_ctl_psi0,
        "final_full_psi0_span": final_full_psi0,
        "final_full_psi_proj_span": final_full_proj,
        "final_full_psi_w0_span": final_full_w0,
        "final_full_psiw0_minus_proj": final_gap,
        "max_abs_emf_vw_c_abs": max_emf,
        "time_at_max_abs_emf_vw_c_abs": t_emf,
        "max_abs_mixed_c2": max_mixed_c2,
        "time_at_max_abs_mixed_c2": t_c2,
        "max_abs_psi0_span": max_psi0,
        "time_at_max_abs_psi0_span": t_psi0,
        "max_abs_psiw0_minus_proj": max_gap,
        "time_at_max_abs_psiw0_minus_proj": t_gap,
        "psi0_ratio_full_over_controlled": float(phase2[case_name]["psi0_ratio_full_over_controlled"]),
        "active_transfer_status": phase2[case_name]["active_transfer_status"],
        "predictive_status": phase2[case_name]["predictive_status"],
        "projected_status": phase2[case_name]["projected_status"],
        "local_status": phase2[case_name]["local_status"],
        "transport_status": phase2[case_name]["transport_status"],
        "ledger_status": phase2[case_name]["ledger_status"],
        "top_transport_component": summary[case_name]["top_abs_component"],
        "top_transport_abs_rate": float(summary[case_name]["top_abs_rate"]),
        "legacy_ledger_max_abs_rate": legacy_max,
        "legacy_ledger_rms_abs_rate": legacy_rms,
        "corrected_ledger_max_abs_rate": corrected_max,
        "corrected_ledger_rms_abs_rate": corrected_rms,
        "ledger_max_reduction_factor": legacy_max / corrected_max if corrected_max > 0.0 else math.inf,
        "ledger_rms_reduction_factor": legacy_rms / corrected_rms if corrected_rms > 0.0 else math.inf,
    }

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-206 Post-Patch Physics Reassessment\n\n")
        f.write(f"- case: `{case_name}`\n")
        f.write(f"- final controlled psi0 span: `{fmt(final_ctl_psi0)}`\n")
        f.write(f"- final full psi0 span: `{fmt(final_full_psi0)}`\n")
        f.write(f"- final full psi_proj span: `{fmt(final_full_proj)}`\n")
        f.write(f"- final full psi_w0 span: `{fmt(final_full_w0)}`\n")
        f.write(f"- final psi_w0 - psi_proj: `{fmt(final_gap)}`\n")
        f.write(f"- psi0 ratio full/control: `{fmt(row['psi0_ratio_full_over_controlled'])}`\n")
        f.write(f"- max |emf_vw_c_abs|: `{fmt(max_emf)}` at `t={fmt(t_emf)}`\n")
        f.write(f"- max |mixed_c2|: `{fmt(max_mixed_c2)}` at `t={fmt(t_c2)}`\n")
        f.write(f"- max |psi0 span|: `{fmt(max_psi0)}` at `t={fmt(t_psi0)}`\n")
        f.write(f"- max |psi_w0 - psi_proj|: `{fmt(max_gap)}` at `t={fmt(t_gap)}`\n\n")
        f.write("## Integrity\n\n")
        f.write(f"- projected/local/transport/ledger: `{row['projected_status']}/{row['local_status']}/{row['transport_status']}/{row['ledger_status']}`\n")
        f.write(f"- active transfer: `{row['active_transfer_status']}`\n")
        f.write(f"- predictive status: `{row['predictive_status']}`\n")
        f.write(f"- top transport: `{row['top_transport_component']}={fmt(row['top_transport_abs_rate'])}`\n\n")
        f.write("## Ledger Improvement\n\n")
        f.write(f"- legacy ledger max abs rate: `{fmt(legacy_max)}`\n")
        f.write(f"- corrected ledger max abs rate: `{fmt(corrected_max)}`\n")
        f.write(f"- legacy ledger rms abs rate: `{fmt(legacy_rms)}`\n")
        f.write(f"- corrected ledger rms abs rate: `{fmt(corrected_rms)}`\n")
        f.write(f"- max reduction factor: `{fmt(row['ledger_max_reduction_factor'])}`\n")
        f.write(f"- rms reduction factor: `{fmt(row['ledger_rms_reduction_factor'])}`\n")


if __name__ == "__main__":
    main()
