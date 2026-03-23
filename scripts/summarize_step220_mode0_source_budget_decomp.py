#!/usr/bin/env python3
"""Summarize Step-220 mode-0 source-budget decomposition."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path

REF_CASE = "patched_env_ref"
FORCED_CASE = "patched_env_train_quad_amp0010"

AGGREGATE_PROXIES = {
    "Src_mode0_total_abs": lambda cols: cols["m4d_int_src_mode0_total_abs"],
    "Src_mode0_plasma_abs": lambda cols: cols["m4d_int_src_mode0_plasma_abs"],
    "Src_mode0_em_abs": lambda cols: cols["m4d_int_src_mode0_em_abs"],
    "Src_mode0_timelike_abs": lambda cols: cols["m4d_int_src_mode0_timelike_abs"],
}

SIGNED_COMPONENTS = {
    "src_momx": "m4d_int_srcmomx_mode0",
    "src_momy": "m4d_int_srcmomy_mode0",
    "src_momz": "m4d_int_srcmomz_mode0",
    "src_momw": "m4d_int_srcmomw_mode0",
    "src_energy": "m4d_int_srcenergy_mode0",
    "src_pi0": "m4d_int_srcpi0_mode0",
    "src_pix": "m4d_int_srcpix_mode0",
    "src_piy": "m4d_int_srcpiy_mode0",
    "src_piz": "m4d_int_srcpiz_mode0",
    "src_piw": "m4d_int_srcpiw_mode0",
    "src_pi0_rhs": "m4d_int_srcpi0_mode0_rhs",
    "src_pi0_mix": "m4d_int_srcpi0_mode0_mix",
    "src_pi0_mix_w0_dynamic": "m4d_int_srcpi0_mode0_mix_w0_dynamic",
    "src_pi0_mix_w0_local_gradient": "m4d_int_srcpi0_mode0_mix_w0_local_gradient",
    "src_pi0_mix_lambda": "m4d_int_srcpi0_mode0_mix_lambda",
    "src_piy_rhs": "m4d_int_srcpiy_mode0_rhs",
    "src_piy_mix": "m4d_int_srcpiy_mode0_mix",
    "src_piy_mix_w0_dynamic": "m4d_int_srcpiy_mode0_mix_w0_dynamic",
    "src_piy_mix_w0_local_gradient": "m4d_int_srcpiy_mode0_mix_w0_local_gradient",
    "src_piy_mix_lambda": "m4d_int_srcpiy_mode0_mix_lambda",
    "src_piw_rhs": "m4d_int_srcpiw_mode0_rhs",
    "src_piw_mix": "m4d_int_srcpiw_mode0_mix",
    "src_piw_mix_w0_dynamic": "m4d_int_srcpiw_mode0_mix_w0_dynamic",
    "src_piw_mix_w0_local_gradient": "m4d_int_srcpiw_mode0_mix_w0_local_gradient",
    "src_piw_mix_lambda": "m4d_int_srcpiw_mode0_mix_lambda",
    "src_em_gauge": "m4d_int_src_em_gauge_mode0",
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


def window_values(times, values, t0, t1=None):
    return [v for t, v in zip(times, values) if t >= t0 and (t1 is None or t < t1) and math.isfinite(v)]


def max_abs(vals):
    if not vals:
        return math.nan
    return max(abs(v) for v in vals)


def rms(vals):
    if not vals:
        return math.nan
    return math.sqrt(sum(v * v for v in vals) / len(vals))


def first_cross_time(times, values, threshold, t0):
    for t, v in zip(times, values):
        if t >= t0 and math.isfinite(v) and abs(v) > threshold:
            return t
    return math.nan


def sample_at_or_after(times, values, t0):
    for t, v in zip(times, values):
        if t >= t0 and math.isfinite(v):
            return v
    return math.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", required=True)
    ap.add_argument("--summary-csv", required=True)
    ap.add_argument("--phase2-csv", required=True)
    ap.add_argument("--gate-on", type=float, required=True)
    ap.add_argument("--gate-off", type=float, required=True)
    ap.add_argument("--eps", type=float, default=1.0e-30)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    summary = {r["case_name"]: r for r in read_csv(Path(args.summary_csv))}
    phase2 = {r["case_name"]: r for r in read_csv(Path(args.phase2_csv))}

    cases = {}
    for case_dir in [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]:
        case = case_dir.name
        full_hst = case_dir / "harris_full.out1.hst"
        if full_hst.exists() and case in (REF_CASE, FORCED_CASE):
            cols = hsm.parse_hst(full_hst)
            cases[case] = {"times": cols["time"], "cols": cols}

    missing = [case for case in (REF_CASE, FORCED_CASE) if case not in cases]
    if missing:
        raise SystemExit(f"Missing required cases: {missing}")

    ref = cases[REF_CASE]
    forced = cases[FORCED_CASE]
    gap = [a - b for a, b in zip(forced["cols"]["m4d_psi_w0_span"], forced["cols"]["m4d_psi_proj_span"])]
    t_release = first_cross_time(forced["times"], gap, 0.0, args.gate_on)
    total_at_release = sample_at_or_after(forced["times"], forced["cols"]["m4d_int_src_mode0_total_abs"], t_release)

    aggregate_rows = []
    for name, getter in AGGREGATE_PROXIES.items():
        ref_vals = getter(ref["cols"])
        forced_vals = getter(forced["cols"])
        ref_env = max_abs(window_values(ref["times"], ref_vals, args.gate_on, None))
        t_cross = first_cross_time(forced["times"], forced_vals, ref_env, args.gate_on)
        aggregate_rows.append({
            "kind": "aggregate",
            "name": name,
            "key": name,
            "ref_env_abs": ref_env,
            "forced_load_max_abs": max_abs(window_values(forced["times"], forced_vals, args.gate_on, args.gate_off)),
            "forced_post_rms_abs": rms(window_values(forced["times"], forced_vals, args.gate_off, None)),
            "t_cross": t_cross,
            "lead": (t_release - t_cross) if math.isfinite(t_release) and math.isfinite(t_cross) else math.nan,
            "value_at_release": sample_at_or_after(forced["times"], forced_vals, t_release),
            "abs_at_release": math.nan,
            "share_of_total_at_release": math.nan,
        })

    component_rows = []
    for name, key in SIGNED_COMPONENTS.items():
        ref_vals = ref["cols"][key]
        forced_vals = forced["cols"][key]
        ref_abs = [abs(v) for v in ref_vals]
        forced_abs = [abs(v) for v in forced_vals]
        ref_env = max_abs(window_values(ref["times"], ref_abs, args.gate_on, None))
        t_cross = first_cross_time(forced["times"], forced_abs, ref_env, args.gate_on)
        value_at_release = sample_at_or_after(forced["times"], forced_vals, t_release)
        abs_at_release = abs(value_at_release) if math.isfinite(value_at_release) else math.nan
        share = abs_at_release / total_at_release if math.isfinite(abs_at_release) and math.isfinite(total_at_release) and total_at_release != 0.0 else math.nan
        component_rows.append({
            "kind": "component",
            "name": name,
            "key": key,
            "ref_env_abs": ref_env,
            "forced_load_max_abs": max_abs(window_values(forced["times"], forced_abs, args.gate_on, args.gate_off)),
            "forced_post_rms_abs": rms(window_values(forced["times"], forced_abs, args.gate_off, None)),
            "t_cross": t_cross,
            "lead": (t_release - t_cross) if math.isfinite(t_release) and math.isfinite(t_cross) else math.nan,
            "value_at_release": value_at_release,
            "abs_at_release": abs_at_release,
            "share_of_total_at_release": share,
        })

    component_rows.sort(
        key=lambda row: (
            math.isnan(row["share_of_total_at_release"]),
            -(row["share_of_total_at_release"] if math.isfinite(row["share_of_total_at_release"]) else -1.0),
        )
    )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "kind",
                "name",
                "key",
                "ref_env_abs",
                "forced_load_max_abs",
                "forced_post_rms_abs",
                "t_cross",
                "lead",
                "value_at_release",
                "abs_at_release",
                "share_of_total_at_release",
            ],
        )
        writer.writeheader()
        for row in aggregate_rows:
            writer.writerow(row)
        for row in component_rows:
            writer.writerow(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-220 Mode-0 Source-Budget Decomposition\n\n")
        f.write("This decomposes the first Step-219 onset clue, `m4d_int_src_mode0_total_abs`, into aggregate buckets and per-channel mode-0 source terms around release onset.\n\n")
        f.write("## Case Summary\n\n")
        f.write("| case | psi0 ratio | projected | local | transport | ledger | top transport |\n")
        f.write("|:---|---:|:---|:---|:---|:---|:---|\n")
        for case_name in (REF_CASE, FORCED_CASE):
            f.write(
                f"| {case_name} | {fmt(phase2[case_name]['psi0_ratio_full_over_controlled'])} | {phase2[case_name]['projected_status']} | {phase2[case_name]['local_status']} | {phase2[case_name]['transport_status']} | {phase2[case_name]['ledger_status']} | {summary[case_name]['top_abs_component']}={fmt(summary[case_name]['top_abs_rate'])} |\n"
            )
        f.write("\n## Release Onset\n\n")
        f.write(f"- `t_R = {fmt(t_release)}` from `R_gap`\n")
        f.write(f"- `Src_mode0_total_abs(t_R) = {fmt(total_at_release)}`\n")
        f.write("\n## Aggregate Buckets\n\n")
        f.write("| bucket | ref envelope | forced load max | forced post rms | crossing time | lead | value at release |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|\n")
        for row in aggregate_rows:
            f.write(
                f"| {row['name']} | {fmt(row['ref_env_abs'])} | {fmt(row['forced_load_max_abs'])} | {fmt(row['forced_post_rms_abs'])} | {fmt(row['t_cross'])} | {fmt(row['lead'])} | {fmt(row['value_at_release'])} |\n"
            )
        f.write("\n## Component Ranking At Release\n\n")
        f.write("| component | abs(value at release) | share of total | crossing time | lead | signed value at release |\n")
        f.write("|:---|---:|---:|---:|---:|---:|\n")
        for row in component_rows:
            f.write(
                f"| {row['name']} | {fmt(row['abs_at_release'])} | {fmt(row['share_of_total_at_release'])} | {fmt(row['t_cross'])} | {fmt(row['lead'])} | {fmt(row['value_at_release'])} |\n"
            )
        leaders = [row for row in component_rows if math.isfinite(row['share_of_total_at_release'])][:8]
        f.write("\n## Leading Components\n\n")
        for row in leaders:
            f.write(
                f"- `{row['name']}` (`{row['key']}`): share `{fmt(row['share_of_total_at_release'])}`, value `{fmt(row['value_at_release'])}`, crossing `{fmt(row['t_cross'])}`\n"
            )
        early = [row for row in component_rows if math.isfinite(row['lead']) and row['lead'] > 0.0]
        coincident = [row for row in component_rows if math.isfinite(row['lead']) and abs(row['lead']) < 1.0e-12]
        if early:
            f.write("\n## Pre-Release Crossings\n\n")
            for row in early:
                f.write(f"- `{row['name']}` leads onset by `{fmt(row['lead'])}`\n")
        else:
            f.write("\n## Pre-Release Crossings\n\n- none\n")
        if coincident:
            f.write("\n## Coincident Crossings\n\n")
            for row in coincident:
                f.write(f"- `{row['name']}` turns on at onset (`t = {fmt(row['t_cross'])}`)\n")


if __name__ == "__main__":
    main()
