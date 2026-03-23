#!/usr/bin/env python3
"""Summarize Step-222 dynamic response-state onset probe."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path

REF_CASE = "patched_env_ref"
FORCED_CASE = "patched_env_train_quad_amp0010"

PROXIES = {
    "R_gap": lambda cols: [a - b for a, b in zip(cols["m4d_psi_w0_span"], cols["m4d_psi_proj_span"])],
    "response_w0_local": lambda cols: cols["m4d_response_w0_local"],
    "response_w0_local_dot": lambda cols: cols["m4d_response_w0_local_dot"],
    "response_w0_local_dz": lambda cols: cols["m4d_response_w0_local_dz"],
    "response_w0_local_gate": lambda cols: cols["m4d_response_w0_local_gate"],
    "response_w0_dot_local_total": lambda cols: cols["m4d_response_w0_dot_local_total"],
    "response_w0_mix_rate_dot": lambda cols: cols["m4d_response_w0_mix_rate_dot"],
    "response_w0_mix_rate_grad": lambda cols: cols["m4d_response_w0_mix_rate_grad"],
    "response_lambda_mix_rate": lambda cols: cols["m4d_response_lambda_mix_rate"],
    "int_response_w0_dynamic_mix_a": lambda cols: cols["m4d_int_response_w0_dynamic_mix_a"],
    "int_response_w0_dynamic_mix_pi": lambda cols: cols["m4d_int_response_w0_dynamic_mix_pi"],
    "int_response_w0_local_gradient_mix_a": lambda cols: cols["m4d_int_response_w0_local_gradient_mix_a"],
    "int_response_w0_local_gradient_mix_pi": lambda cols: cols["m4d_int_response_w0_local_gradient_mix_pi"],
    "srcpiw_mode0_mix_w0_dynamic": lambda cols: cols["m4d_int_srcpiw_mode0_mix_w0_dynamic"],
    "srcpiy_mode0_mix_w0_dynamic": lambda cols: cols["m4d_int_srcpiy_mode0_mix_w0_dynamic"],
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
    return [v for t, v in zip(times, vals) if t >= t0 and (t1 is None or t < t1) and math.isfinite(v)]


def first_cross_time(times, values, threshold, t0):
    for t, v in zip(times, values):
        if t >= t0 and math.isfinite(v) and abs(v) > threshold:
            return t
    return math.nan


def first_nonzero_time(times, values, t0, atol=1.0e-30):
    for t, v in zip(times, values):
        if t >= t0 and math.isfinite(v) and abs(v) > atol:
            return t
    return math.nan


def sample_at_or_after(times, values, t0):
    for t, v in zip(times, values):
        if t >= t0 and math.isfinite(v):
            return v
    return math.nan


def mean_dt(times):
    diffs = [b - a for a, b in zip(times[:-1], times[1:]) if math.isfinite(a) and math.isfinite(b)]
    return (sum(diffs) / len(diffs)) if diffs else math.nan


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
    release = PROXIES["R_gap"](forced["cols"])
    t_release = first_cross_time(times, release, 0.0, args.gate_on)
    onset_lo = t_release - 2.0 * dt if math.isfinite(t_release) and math.isfinite(dt) else args.gate_on
    onset_hi = t_release + 2.0 * dt if math.isfinite(t_release) and math.isfinite(dt) else args.gate_off

    rows = []
    ordered = []
    for name, getter in PROXIES.items():
        ref_vals = getter(ref["cols"])
        forced_vals = getter(forced["cols"])
        ref_env = max_abs(window(ref["times"], ref_vals, args.gate_on, None))
        t_cross = first_cross_time(times, forced_vals, ref_env, args.gate_on)
        t_nonzero = first_nonzero_time(times, forced_vals, args.gate_on)
        lead = (t_release - t_cross) if math.isfinite(t_release) and math.isfinite(t_cross) else math.nan
        rows.append({
            "proxy": name,
            "ref_env_abs": ref_env,
            "forced_load_max_abs": max_abs(window(times, forced_vals, args.gate_on, args.gate_off)),
            "forced_onset_rms": rms(window(times, forced_vals, onset_lo, onset_hi)),
            "forced_post_rms": rms(window(times, forced_vals, args.gate_off, None)),
            "t_nonzero": t_nonzero,
            "t_cross": t_cross,
            "lead": lead,
            "value_at_release": sample_at_or_after(times, forced_vals, t_release),
        })
        if math.isfinite(t_cross):
            ordered.append((name, t_cross, lead))

    ordered.sort(key=lambda item: item[1])
    pre_release = [row for row in ordered if row[0] != "R_gap" and math.isfinite(row[2]) and row[2] > 0.0]
    coincident = [row for row in ordered if row[0] != "R_gap" and math.isfinite(row[2]) and abs(row[2]) <= (dt if math.isfinite(dt) else 0.0)]

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "proxy",
                "ref_env_abs",
                "forced_load_max_abs",
                "forced_onset_rms",
                "forced_post_rms",
                "t_nonzero",
                "t_cross",
                "lead",
                "value_at_release",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-222 Dynamic Response-State Probe\n\n")
        f.write("This rerun instruments the signed local dynamic `w0` response state and signed dynamic/local-gradient mixing budgets to see whether any upstream state leads the mode-0 EM source onset.\n\n")
        f.write("## Case Summary\n\n")
        f.write("| case | psi0 ratio | projected | local | transport | ledger | top transport |\n")
        f.write("|:---|---:|:---|:---|:---|:---|:---|\n")
        for case_name in (REF_CASE, FORCED_CASE):
            f.write(
                f"| {case_name} | {fmt(phase2[case_name]['psi0_ratio_full_over_controlled'])} | {phase2[case_name]['projected_status']} | {phase2[case_name]['local_status']} | {phase2[case_name]['transport_status']} | {phase2[case_name]['ledger_status']} | {summary[case_name]['top_abs_component']}={fmt(summary[case_name]['top_abs_rate'])} |\n"
            )
        f.write("\n## Timing\n\n")
        f.write(f"- gate-on `t_on = {fmt(args.gate_on)}`\n")
        f.write(f"- gate-off `t_off = {fmt(args.gate_off)}`\n")
        f.write(f"- mean history dt `= {fmt(dt)}`\n")
        f.write(f"- release onset `t_R = {fmt(t_release)}` from `R_gap`\n")
        f.write("\n## Proxy Windows\n\n")
        f.write("| proxy | ref envelope abs | forced load max abs | forced onset rms | forced post rms | t_nonzero | t_cross | lead | value at release |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            f.write(
                f"| {row['proxy']} | {fmt(row['ref_env_abs'])} | {fmt(row['forced_load_max_abs'])} | {fmt(row['forced_onset_rms'])} | {fmt(row['forced_post_rms'])} | {fmt(row['t_nonzero'])} | {fmt(row['t_cross'])} | {fmt(row['lead'])} | {fmt(row['value_at_release'])} |\n"
            )
        f.write("\n## Crossing Order\n\n")
        if ordered:
            for name, t_cross, lead in ordered:
                f.write(f"- `{name}` crosses at `t = {fmt(t_cross)}` with lead `t_R - t = {fmt(lead)}`\n")
        else:
            f.write("- no forced proxy exceeded its clean-reference envelope\n")
        f.write("\n## Assessment\n\n")
        if pre_release:
            f.write("Leading proxies before release onset:\n")
            for name, t_cross, lead in pre_release:
                f.write(f"- `{name}` leads by `{fmt(lead)}`\n")
        else:
            f.write("- no signed dynamic-response proxy leads release onset\n")
        if coincident:
            f.write("\nCoincident proxies at onset:\n")
            for name, t_cross, lead in coincident:
                f.write(f"- `{name}` is coincident within one history step (`lead = {fmt(lead)}`)\n")
        f.write("\n## Interpretation Guide\n\n")
        f.write("- `response_w0_local*`: signed local forcing state and slope\n")
        f.write("- `response_w0_mix_rate_*`: signed dynamic vs local-gradient mixing rates before source conversion\n")
        f.write("- `int_response_w0_*`: signed integrated dynamic/local-gradient mode-mixing budgets\n")
        f.write("- `srcpiw/srcpiy_mode0_mix_w0_dynamic`: downstream EM source sinks that already turned on at onset in Step-221\n")


if __name__ == "__main__":
    main()
