#!/usr/bin/env python3
"""Summarize Step-213 near-threshold sub-amp010 forcing sweep."""

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


def max_abs_with_time(times, values, tmin=0.0):
    pairs = [(abs(v), t) for t, v in zip(times, values) if t >= tmin and math.isfinite(v)]
    if not pairs:
        return (math.nan, math.nan)
    val, t = max(pairs, key=lambda p: p[0])
    return val, t


def rms_after(times, values, tmin=0.0):
    vals = [v for t, v in zip(times, values) if t >= tmin and math.isfinite(v)]
    if not vals:
        return math.nan
    return math.sqrt(sum(v * v for v in vals) / len(vals))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", required=True)
    ap.add_argument("--summary-csv", required=True)
    ap.add_argument("--phase2-csv", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    summary = {r["case_name"]: r for r in read_csv(Path(args.summary_csv))}
    phase2 = {r["case_name"]: r for r in read_csv(Path(args.phase2_csv))}

    rows = []
    for case_dir in [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]:
        case = case_dir.name
        full_hst = case_dir / "harris_full.out1.hst"
        if not full_hst.exists():
            continue
        full = hsm.parse_hst(full_hst)
        times = full["time"]
        gap = [a - b for a, b in zip(full["m4d_psi_w0_span"], full["m4d_psi_proj_span"])]
        late_gap, t_late_gap = max_abs_with_time(times, gap, tmin=0.5)
        row = {
            "case_name": case,
            "psi0_ratio_full_over_controlled": float(phase2[case]["psi0_ratio_full_over_controlled"]),
            "final_gap": gap[-1],
            "late_max_abs_gap": late_gap,
            "time_at_late_max_abs_gap": t_late_gap,
            "late_rms_gap": rms_after(times, gap, tmin=0.5),
            "projected_status": phase2[case]["projected_status"],
            "local_status": phase2[case]["local_status"],
            "transport_status": phase2[case]["transport_status"],
            "ledger_status": phase2[case]["ledger_status"],
            "top_transport_component": summary[case]["top_abs_component"],
            "top_transport_abs_rate": float(summary[case]["top_abs_rate"]),
        }
        rows.append(row)

    if not rows:
        raise SystemExit("No Step-213 rows computed.")

    rows.sort(key=lambda r: r["case_name"])
    pass_rows = [r for r in rows if r["transport_status"] == "PASS"]
    best_pass = max(pass_rows, key=lambda r: r["late_rms_gap"]) if pass_rows else None
    best_fail = max(rows, key=lambda r: r["late_rms_gap"])

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-213 Post-Patch Low-Amplitude Threshold Sweep\n\n")
        f.write("| case | psi0 ratio | final gap | late max gap | late gap rms | transport | top transport | ledger |\n")
        f.write("|:---|---:|---:|---:|---:|:---|:---|:---|\n")
        for r in rows:
            f.write(
                f"| {r['case_name']} | {fmt(r['psi0_ratio_full_over_controlled'])} | {fmt(r['final_gap'])} | "
                f"{fmt(r['late_max_abs_gap'])} | {fmt(r['late_rms_gap'])} | {r['transport_status']} | "
                f"{r['top_transport_component']}={fmt(r['top_transport_abs_rate'])} | {r['ledger_status']} |\n"
            )
        f.write("\n")
        if best_pass is not None:
            f.write("## Best Transport-Clean Case\n\n")
            f.write(
                f"- case: `{best_pass['case_name']}`\n"
                f"- late rms `psi_w0 - psi_proj`: `{fmt(best_pass['late_rms_gap'])}`\n"
                f"- final `psi_w0 - psi_proj`: `{fmt(best_pass['final_gap'])}`\n\n"
            )
        else:
            f.write("## Best Transport-Clean Case\n\n- none\n\n")
        f.write("## Strongest Forced Case\n\n")
        f.write(
            f"- case: `{best_fail['case_name']}`\n"
            f"- late rms `psi_w0 - psi_proj`: `{fmt(best_fail['late_rms_gap'])}`\n"
            f"- transport: `{best_fail['transport_status']}`\n"
            f"- top transport: `{best_fail['top_transport_component']}={fmt(best_fail['top_transport_abs_rate'])}`\n"
        )
