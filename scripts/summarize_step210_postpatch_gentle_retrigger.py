#!/usr/bin/env python3
"""Summarize gentle post-patch omega retrigger sweep."""

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


def count_local_peaks(times, values, tmin=0.2, min_abs=1.0e-2):
    count = 0
    for i in range(1, len(values) - 1):
        if times[i] < tmin:
            continue
        v = abs(values[i])
        if v < min_abs:
            continue
        if v >= abs(values[i - 1]) and v >= abs(values[i + 1]):
            count += 1
    return count


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
        ctl_hst = case_dir / "harris_controlled.out1.hst"
        if not full_hst.exists() or not ctl_hst.exists():
            continue
        full = hsm.parse_hst(full_hst)
        ctl = hsm.parse_hst(ctl_hst)
        times = full["time"]
        gap = [a - b for a, b in zip(full["m4d_psi_w0_span"], full["m4d_psi_proj_span"])]
        max_gap, t_gap = max_abs_with_time(times, gap)
        late_gap, t_late_gap = max_abs_with_time(times, gap, tmin=0.5)
        late_c2, t_late_c2 = max_abs_with_time(times, full["m4d_mixed_c2"], tmin=0.5)
        peaks = count_local_peaks(times, gap)
        corrected_max, corrected_rms, corrected_final = hsm.em_bulk_ledger_abs_rate_metrics(
            times, full["m4d_em_u_bulk"], full["m4d_int_ja_ea"], full["m4d_int_jw_ew"]
        )
        legacy_max, legacy_rms, legacy_final = hsm.em_bulk_ledger_abs_rate_metrics(
            times, full["m4d_em_u_bulk_legacy"], full["m4d_int_ja_ea"], full["m4d_int_jw_ew"]
        )
        rows.append({
            "case_name": case,
            "psi0_ratio_full_over_controlled": float(phase2[case]["psi0_ratio_full_over_controlled"]),
            "final_gap": gap[-1],
            "max_abs_gap": max_gap,
            "time_at_max_abs_gap": t_gap,
            "late_max_abs_gap": late_gap,
            "time_at_late_max_abs_gap": t_late_gap,
            "late_rms_gap": rms_after(times, gap, tmin=0.5),
            "postburst_peak_count": peaks,
            "late_max_abs_mixed_c2": late_c2,
            "time_at_late_max_abs_mixed_c2": t_late_c2,
            "projected_status": phase2[case]["projected_status"],
            "local_status": phase2[case]["local_status"],
            "transport_status": phase2[case]["transport_status"],
            "ledger_status": phase2[case]["ledger_status"],
            "active_transfer_status": phase2[case]["active_transfer_status"],
            "predictive_status": phase2[case]["predictive_status"],
            "top_transport_component": summary[case]["top_abs_component"],
            "top_transport_abs_rate": float(summary[case]["top_abs_rate"]),
            "legacy_ledger_max_abs_rate": legacy_max,
            "legacy_ledger_rms_abs_rate": legacy_rms,
            "legacy_ledger_final_rate": legacy_final,
            "corrected_ledger_max_abs_rate": corrected_max,
            "corrected_ledger_rms_abs_rate": corrected_rms,
            "corrected_ledger_final_rate": corrected_final,
        })

    if not rows:
        raise SystemExit("No Step-210 rows computed.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-210 Gentle Post-Patch Omega Retrigger Sweep\n\n")
        f.write("| case | psi0 ratio | final gap | late max gap | late gap rms | postburst peaks | late max C2 | proj/local/transport/ledger |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|:---|\n")
        for r in rows:
            f.write(
                f"| {r['case_name']} | {fmt(r['psi0_ratio_full_over_controlled'])} | {fmt(r['final_gap'])} | "
                f"{fmt(r['late_max_abs_gap'])} | {fmt(r['late_rms_gap'])} | {int(r['postburst_peak_count'])} | "
                f"{fmt(r['late_max_abs_mixed_c2'])} | "
                f"{r['projected_status']}/{r['local_status']}/{r['transport_status']}/{r['ledger_status']} |\n"
            )
        f.write("\n## Details\n\n")
        for r in rows:
            f.write(f"### {r['case_name']}\n\n")
            f.write(f"- psi0 ratio full/control: `{fmt(r['psi0_ratio_full_over_controlled'])}`\n")
            f.write(f"- final `psi_w0 - psi_proj`: `{fmt(r['final_gap'])}`\n")
            f.write(f"- peak `|psi_w0 - psi_proj|`: `{fmt(r['max_abs_gap'])}` at `t={fmt(r['time_at_max_abs_gap'])}`\n")
            f.write(f"- late peak `|psi_w0 - psi_proj|`: `{fmt(r['late_max_abs_gap'])}` at `t={fmt(r['time_at_late_max_abs_gap'])}`\n")
            f.write(f"- late rms `psi_w0 - psi_proj`: `{fmt(r['late_rms_gap'])}`\n")
            f.write(f"- postburst peak count: `{r['postburst_peak_count']}`\n")
            f.write(f"- late max `|mixed_c2|`: `{fmt(r['late_max_abs_mixed_c2'])}` at `t={fmt(r['time_at_late_max_abs_mixed_c2'])}`\n")
            f.write(f"- top transport: `{r['top_transport_component']}={fmt(r['top_transport_abs_rate'])}`\n")
            f.write(f"- corrected ledger rms/max: `{fmt(r['corrected_ledger_rms_abs_rate'])}` / `{fmt(r['corrected_ledger_max_abs_rate'])}`\n\n")


if __name__ == "__main__":
    main()
