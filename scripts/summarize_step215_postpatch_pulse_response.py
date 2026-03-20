#!/usr/bin/env python3
"""Summarize structured pulse forcing on the patched branch."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


CASE_META = {
    "patched_static_pulse_ref": {"t_on": None, "t_off": None},
    "patched_pulse_amp0025_short": {"t_on": 0.35, "t_off": 0.40},
    "patched_pulse_amp0025_tau": {"t_on": 0.35, "t_off": 0.65},
    "patched_pulse_amp0025_pair": {"t_on": 0.35, "t_off": 0.75},
    "patched_pulse_amp0025_train": {"t_on": 0.35, "t_off": 0.80},
}


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


def max_abs_with_time(times, values, tmin=0.0, tmax=None):
    pairs = []
    for t, v in zip(times, values):
        if t < tmin:
            continue
        if tmax is not None and t > tmax:
            continue
        if math.isfinite(v):
            pairs.append((abs(v), t))
    if not pairs:
        return (math.nan, math.nan)
    val, t = max(pairs, key=lambda p: p[0])
    return val, t


def rms_between(times, values, tmin=0.0, tmax=None):
    vals = []
    for t, v in zip(times, values):
        if t < tmin:
            continue
        if tmax is not None and t > tmax:
            continue
        if math.isfinite(v):
            vals.append(v)
    if not vals:
        return math.nan
    return math.sqrt(sum(v * v for v in vals) / len(vals))


def count_peaks(times, values, tmin=0.0):
    count = 0
    for i in range(1, len(values) - 1):
        if times[i] < tmin:
            continue
        v = abs(values[i])
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
        if not full_hst.exists() or case not in CASE_META:
            continue
        full = hsm.parse_hst(full_hst)
        times = full["time"]
        gap = [a - b for a, b in zip(full["m4d_psi_w0_span"], full["m4d_psi_proj_span"])]
        meta = CASE_META[case]
        t_on = meta["t_on"]
        t_off = meta["t_off"]
        if t_on is None:
            during_max = math.nan
            post_max, t_post_max = max_abs_with_time(times, gap, tmin=0.5)
            post_rms = rms_between(times, gap, tmin=0.5)
            memory_ratio = math.nan
            post_peaks = count_peaks(times, gap, tmin=0.5)
        else:
            during_max, _ = max_abs_with_time(times, gap, tmin=t_on, tmax=t_off)
            post_max, t_post_max = max_abs_with_time(times, gap, tmin=t_off)
            post_rms = rms_between(times, gap, tmin=t_off)
            memory_ratio = (post_max / during_max) if (during_max and math.isfinite(during_max) and during_max > 0.0) else math.nan
            post_peaks = count_peaks(times, gap, tmin=t_off)
        rows.append({
            "case_name": case,
            "psi0_ratio_full_over_controlled": float(phase2[case]["psi0_ratio_full_over_controlled"]),
            "final_gap": gap[-1],
            "post_force_max_abs_gap": post_max,
            "time_at_post_force_max_abs_gap": t_post_max,
            "post_force_rms_gap": post_rms,
            "memory_ratio": memory_ratio,
            "post_force_peak_count": post_peaks,
            "projected_status": phase2[case]["projected_status"],
            "local_status": phase2[case]["local_status"],
            "transport_status": phase2[case]["transport_status"],
            "ledger_status": phase2[case]["ledger_status"],
            "active_transfer_status": phase2[case]["active_transfer_status"],
            "predictive_status": phase2[case]["predictive_status"],
            "top_transport_component": summary[case]["top_abs_component"],
            "top_transport_abs_rate": float(summary[case]["top_abs_rate"]),
        })

    if not rows:
        raise SystemExit("No Step-215 rows computed.")

    rows.sort(key=lambda r: r["case_name"])
    forced_rows = [r for r in rows if r["case_name"] != "patched_static_pulse_ref"]
    best_clean = [r for r in forced_rows if r["transport_status"] == "PASS"]
    best_clean = max(best_clean, key=lambda r: r["post_force_rms_gap"]) if best_clean else None
    best_forced = max(forced_rows, key=lambda r: r["post_force_rms_gap"])

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-215 Structured Pulse Forcing Response\n\n")
        f.write("| case | psi0 ratio | final gap | post-force max gap | post-force rms gap | memory ratio | post-force peaks | transport | top transport | ledger |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|:---|:---|:---|\n")
        for r in rows:
            f.write(
                f"| {r['case_name']} | {fmt(r['psi0_ratio_full_over_controlled'])} | {fmt(r['final_gap'])} | "
                f"{fmt(r['post_force_max_abs_gap'])} | {fmt(r['post_force_rms_gap'])} | {fmt(r['memory_ratio'])} | "
                f"{int(r['post_force_peak_count'])} | {r['transport_status']} | "
                f"{r['top_transport_component']}={fmt(r['top_transport_abs_rate'])} | {r['ledger_status']} |\n"
            )
        f.write("\n## Best Transport-Clean Forced Case\n\n")
        if best_clean is None:
            f.write("- none\n\n")
        else:
            f.write(
                f"- case: `{best_clean['case_name']}`\n"
                f"- post-force rms `psi_w0 - psi_proj`: `{fmt(best_clean['post_force_rms_gap'])}`\n"
                f"- final `psi_w0 - psi_proj`: `{fmt(best_clean['final_gap'])}`\n\n"
            )
        f.write("## Strongest Forced Case\n\n")
        f.write(
            f"- case: `{best_forced['case_name']}`\n"
            f"- post-force rms `psi_w0 - psi_proj`: `{fmt(best_forced['post_force_rms_gap'])}`\n"
            f"- memory ratio: `{fmt(best_forced['memory_ratio'])}`\n"
            f"- transport: `{best_forced['transport_status']}`\n"
            f"- top transport: `{best_forced['top_transport_component']}={fmt(best_forced['top_transport_abs_rate'])}`\n"
        )


if __name__ == "__main__":
    main()
