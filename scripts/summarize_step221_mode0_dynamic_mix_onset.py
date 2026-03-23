#!/usr/bin/env python3
"""Summarize Step-221 mode-0 dynamic-mix onset decomposition."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path

REF_CASE = "patched_env_ref"
FORCED_CASE = "patched_env_train_quad_amp0010"

SOURCE_GROUPS = {
    "piw": {
        "total": "m4d_int_srcpiw_mode0",
        "rhs": "m4d_int_srcpiw_mode0_rhs",
        "mix": "m4d_int_srcpiw_mode0_mix",
        "mix_w0_dynamic": "m4d_int_srcpiw_mode0_mix_w0_dynamic",
        "mix_w0_local_gradient": "m4d_int_srcpiw_mode0_mix_w0_local_gradient",
        "mix_lambda": "m4d_int_srcpiw_mode0_mix_lambda",
    },
    "piy": {
        "total": "m4d_int_srcpiy_mode0",
        "rhs": "m4d_int_srcpiy_mode0_rhs",
        "mix": "m4d_int_srcpiy_mode0_mix",
        "mix_w0_dynamic": "m4d_int_srcpiy_mode0_mix_w0_dynamic",
        "mix_w0_local_gradient": "m4d_int_srcpiy_mode0_mix_w0_local_gradient",
        "mix_lambda": "m4d_int_srcpiy_mode0_mix_lambda",
    },
}

UPSTREAM_PROXIES = {
    "response_w0_dynamic_mix_a_abs": lambda cols: cols["m4d_int_response_w0_dynamic_mix_a_abs"],
    "response_w0_dynamic_mix_pi_abs": lambda cols: cols["m4d_int_response_w0_dynamic_mix_pi_abs"],
    "response_w0_local_gradient_mix_a_abs": lambda cols: cols["m4d_int_response_w0_local_gradient_mix_a_abs"],
    "response_w0_local_gradient_mix_pi_abs": lambda cols: cols["m4d_int_response_w0_local_gradient_mix_pi_abs"],
    "response_lambda_dynamic_mix_a_abs": lambda cols: cols["m4d_int_response_lambda_dynamic_mix_a_abs"],
    "response_lambda_dynamic_mix_pi_abs": lambda cols: cols["m4d_int_response_lambda_dynamic_mix_pi_abs"],
    "aw_mode1_gradz_power_mix_w0_dynamic_abs": lambda cols: [abs(v) for v in cols["m4d_int_aw_mode1_gradz_power_mix_w0_dynamic"]],
    "aw_mode2_gradz_power_mix_w0_dynamic_abs": lambda cols: [abs(v) for v in cols["m4d_int_aw_mode2_gradz_power_mix_w0_dynamic"]],
    "aw_gradz_power_mix_w0_dynamic_sum_abs": lambda cols: [
        abs(a) + abs(b)
        for a, b in zip(
            cols["m4d_int_aw_mode1_gradz_power_mix_w0_dynamic"],
            cols["m4d_int_aw_mode2_gradz_power_mix_w0_dynamic"],
        )
    ],
    "aw_mode1_gradz_power_mix_w0_local_gradient_abs": lambda cols: [abs(v) for v in cols["m4d_int_aw_mode1_gradz_power_mix_w0_local_gradient"]],
    "aw_mode2_gradz_power_mix_w0_local_gradient_abs": lambda cols: [abs(v) for v in cols["m4d_int_aw_mode2_gradz_power_mix_w0_local_gradient"]],
    "aw_gradz_power_mix_w0_local_gradient_sum_abs": lambda cols: [
        abs(a) + abs(b)
        for a, b in zip(
            cols["m4d_int_aw_mode1_gradz_power_mix_w0_local_gradient"],
            cols["m4d_int_aw_mode2_gradz_power_mix_w0_local_gradient"],
        )
    ],
    "piw_mode0_abs": lambda cols: [abs(v) for v in cols["m4d_piw_mode_0"]],
    "piy_mode0_abs": lambda cols: [abs(v) for v in cols["m4d_piy_mode_0"]],
    "divpiw_mode0_abs": lambda cols: [abs(v) for v in cols["m4d_int_divpiw_mode0"]],
    "divpiy_mode0_abs": lambda cols: [abs(v) for v in cols["m4d_int_divpiy_mode0"]],
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


def max_abs(vals):
    vals = [abs(v) for v in vals if math.isfinite(v)]
    if not vals:
        return math.nan
    return max(vals)


def rms(vals):
    vals = [v for v in vals if math.isfinite(v)]
    if not vals:
        return math.nan
    return math.sqrt(sum(v * v for v in vals) / len(vals))


def window(times, vals, t0, t1=None):
    return [v for t, v in zip(times, vals) if t >= t0 and (t1 is None or t < t1)]


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


def mean_dt(times):
    if len(times) < 2:
        return math.nan
    diffs = [b - a for a, b in zip(times[:-1], times[1:]) if math.isfinite(a) and math.isfinite(b)]
    return sum(diffs) / len(diffs) if diffs else math.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", required=True)
    ap.add_argument("--summary-csv", required=True)
    ap.add_argument("--phase2-csv", required=True)
    ap.add_argument("--gate-on", type=float, required=True)
    ap.add_argument("--gate-off", type=float, required=True)
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
    times = forced["times"]
    dt = mean_dt(times)
    gap = [a - b for a, b in zip(forced["cols"]["m4d_psi_w0_span"], forced["cols"]["m4d_psi_proj_span"])]
    t_release = first_cross_time(times, gap, 0.0, args.gate_on)
    onset_lo = t_release - 2.0 * dt if math.isfinite(t_release) and math.isfinite(dt) else args.gate_on
    onset_hi = t_release + 2.0 * dt if math.isfinite(t_release) and math.isfinite(dt) else args.gate_off

    source_rows = []
    for axis, pieces in SOURCE_GROUPS.items():
        total_key = pieces["total"]
        total_at_release = abs(sample_at_or_after(times, forced["cols"][total_key], t_release))
        for label, key in pieces.items():
            ref_abs = [abs(v) for v in ref["cols"][key]]
            forced_signed = forced["cols"][key]
            forced_abs = [abs(v) for v in forced_signed]
            ref_env = max_abs(window(ref["times"], ref_abs, args.gate_on, None))
            t_cross = first_cross_time(times, forced_abs, ref_env, args.gate_on)
            value_at_release = sample_at_or_after(times, forced_signed, t_release)
            abs_at_release = abs(value_at_release) if math.isfinite(value_at_release) else math.nan
            source_rows.append({
                "group": axis,
                "label": label,
                "key": key,
                "ref_env_abs": ref_env,
                "forced_load_max_abs": max_abs(window(times, forced_abs, args.gate_on, args.gate_off)),
                "forced_onset_rms_abs": rms(window(times, forced_abs, onset_lo, onset_hi)),
                "forced_post_rms_abs": rms(window(times, forced_abs, args.gate_off, None)),
                "t_cross": t_cross,
                "lead": (t_release - t_cross) if math.isfinite(t_release) and math.isfinite(t_cross) else math.nan,
                "value_at_release": value_at_release,
                "abs_at_release": abs_at_release,
                "share_of_axis_total": abs_at_release / total_at_release if math.isfinite(abs_at_release) and math.isfinite(total_at_release) and total_at_release != 0.0 else math.nan,
            })

    source_rows.sort(
        key=lambda r: (
            r["group"],
            math.isnan(r["share_of_axis_total"]),
            -(r["share_of_axis_total"] if math.isfinite(r["share_of_axis_total"]) else -1.0),
        )
    )

    upstream_rows = []
    for name, getter in UPSTREAM_PROXIES.items():
        ref_vals = getter(ref["cols"])
        forced_vals = getter(forced["cols"])
        ref_env = max_abs(window(ref["times"], ref_vals, args.gate_on, None))
        t_cross = first_cross_time(times, forced_vals, ref_env, args.gate_on)
        upstream_rows.append({
            "name": name,
            "ref_env_abs": ref_env,
            "forced_load_max_abs": max_abs(window(times, forced_vals, args.gate_on, args.gate_off)),
            "forced_onset_rms_abs": rms(window(times, forced_vals, onset_lo, onset_hi)),
            "forced_post_rms_abs": rms(window(times, forced_vals, args.gate_off, None)),
            "t_cross": t_cross,
            "lead": (t_release - t_cross) if math.isfinite(t_release) and math.isfinite(t_cross) else math.nan,
            "value_at_release": sample_at_or_after(times, forced_vals, t_release),
        })

    upstream_rows.sort(
        key=lambda r: (
            math.isnan(r["t_cross"]),
            r["t_cross"] if math.isfinite(r["t_cross"]) else math.inf,
            -r["forced_load_max_abs"],
        )
    )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "section",
                "group",
                "label",
                "name",
                "key",
                "ref_env_abs",
                "forced_load_max_abs",
                "forced_onset_rms_abs",
                "forced_post_rms_abs",
                "t_cross",
                "lead",
                "value_at_release",
                "abs_at_release",
                "share_of_axis_total",
            ],
        )
        writer.writeheader()
        for row in source_rows:
            writer.writerow({"section": "source", "name": "", **row})
        for row in upstream_rows:
            writer.writerow(
                {
                    "section": "upstream",
                    "group": "",
                    "label": "",
                    "key": "",
                    "abs_at_release": "",
                    "share_of_axis_total": "",
                    **row,
                }
            )

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-221 Mode-0 Dynamic-Mix Onset Decomposition\n\n")
        f.write(
            "This decomposes the Step-220 onset clue further, focusing on the `piw` and `piy` mode-0 dynamic-mixing source paths and nearby upstream response proxies.\n\n"
        )
        f.write("## Case Summary\n\n")
        f.write("| case | psi0 ratio | projected | local | transport | ledger | top transport |\n")
        f.write("|:---|---:|:---|:---|:---|:---|:---|\n")
        for case_name in (REF_CASE, FORCED_CASE):
            f.write(
                f"| {case_name} | {fmt(phase2[case_name]['psi0_ratio_full_over_controlled'])} | {phase2[case_name]['projected_status']} | {phase2[case_name]['local_status']} | {phase2[case_name]['transport_status']} | {phase2[case_name]['ledger_status']} | {summary[case_name]['top_abs_component']}={fmt(summary[case_name]['top_abs_rate'])} |\n"
            )
        f.write("\n## Release Onset\n\n")
        f.write(f"- `t_R = {fmt(t_release)}` from `R_gap`\n")
        f.write(f"- onset window `[t_R-2dt, t_R+2dt] = [{fmt(onset_lo)}, {fmt(onset_hi)}]` with `dt = {fmt(dt)}`\n")
        for axis in ("piw", "piy"):
            total_key = SOURCE_GROUPS[axis]["total"]
            total_val = sample_at_or_after(times, forced["cols"][total_key], t_release)
            f.write(f"- `{axis}` total at release: `{fmt(total_val)}`\n")

        for axis in ("piw", "piy"):
            rows = [r for r in source_rows if r["group"] == axis]
            f.write(f"\n## {axis.upper()} Source Split\n\n")
            f.write("| piece | share of axis total | value at release | ref envelope | crossing time | lead | onset rms | post rms |\n")
            f.write("|:---|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in rows:
                f.write(
                    f"| {row['label']} | {fmt(row['share_of_axis_total'])} | {fmt(row['value_at_release'])} | {fmt(row['ref_env_abs'])} | {fmt(row['t_cross'])} | {fmt(row['lead'])} | {fmt(row['forced_onset_rms_abs'])} | {fmt(row['forced_post_rms_abs'])} |\n"
                )

        f.write("\n## Upstream Proxy Timing\n\n")
        f.write("| proxy | ref envelope | value at release | crossing time | lead | onset rms | post rms |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|\n")
        for row in upstream_rows:
            f.write(
                f"| {row['name']} | {fmt(row['ref_env_abs'])} | {fmt(row['value_at_release'])} | {fmt(row['t_cross'])} | {fmt(row['lead'])} | {fmt(row['forced_onset_rms_abs'])} | {fmt(row['forced_post_rms_abs'])} |\n"
            )

        leading = [r for r in upstream_rows if math.isfinite(r["lead"]) and r["lead"] > 0.0]
        coincident = [r for r in upstream_rows if math.isfinite(r["lead"]) and abs(r["lead"]) < 1.0e-12]
        f.write("\n## Timing Interpretation\n\n")
        if leading:
            f.write("Pre-release upstream crossings:\n")
            for row in leading:
                f.write(f"- `{row['name']}` leads onset by `{fmt(row['lead'])}`\n")
        else:
            f.write("- no upstream proxy crosses before release onset\n")
        if coincident:
            f.write("Coincident upstream crossings:\n")
            for row in coincident[:10]:
                f.write(f"- `{row['name']}` turns on at onset (`t = {fmt(row['t_cross'])}`)\n")


if __name__ == "__main__":
    main()
