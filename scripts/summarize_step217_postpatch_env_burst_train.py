#!/usr/bin/env python3
"""Summarize the post-patch environment-style burst-train sweep."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


CASE_META = {
    'patched_env_ref': {'t_off': None},
    'patched_env_train_triplet_amp0005': {'t_off': 0.90},
    'patched_env_train_triplet_amp0010': {'t_off': 0.90},
    'patched_env_train_quad_amp0005': {'t_off': 0.90},
    'patched_env_train_quad_amp0010': {'t_off': 0.90},
}


def load_hsm(script_dir: Path):
    module_path = script_dir / 'harris_scan_matrix.py'
    spec = importlib.util.spec_from_file_location('harris_scan_matrix', module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Could not load module spec from {module_path}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return 'nan'
    try:
        return f'{float(v):.6e}'
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
    ap.add_argument('--glob', required=True)
    ap.add_argument('--summary-csv', required=True)
    ap.add_argument('--phase2-csv', required=True)
    ap.add_argument('--out-csv', required=True)
    ap.add_argument('--out-md', required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    summary = {r['case_name']: r for r in read_csv(Path(args.summary_csv))}
    phase2 = {r['case_name']: r for r in read_csv(Path(args.phase2_csv))}

    rows = []
    for case_dir in [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]:
        case = case_dir.name
        full_hst = case_dir / 'harris_full.out1.hst'
        if not full_hst.exists() or case not in CASE_META:
            continue
        full = hsm.parse_hst(full_hst)
        times = full['time']
        gap = [a - b for a, b in zip(full['m4d_psi_w0_span'], full['m4d_psi_proj_span'])]
        t_off = CASE_META[case]['t_off']
        if t_off is None:
            post_max = 0.0
            t_post = math.nan
            post_rms = 0.0
            memory_ratio = math.nan
        else:
            during_max, _ = max_abs_with_time(times, gap, tmin=0.30)
            post_max, t_post = max_abs_with_time(times, gap, tmin=t_off)
            post_rms = rms_after(times, gap, tmin=t_off)
            memory_ratio = (post_max / during_max) if during_max > 0.0 else math.nan
        rows.append({
            'case_name': case,
            'psi0_ratio_full_over_controlled': float(phase2[case]['psi0_ratio_full_over_controlled']),
            'final_gap': gap[-1],
            'post_force_max_abs_gap': post_max,
            'time_at_post_force_max_abs_gap': t_post,
            'post_force_rms_gap': post_rms,
            'memory_ratio': memory_ratio,
            'projected_status': phase2[case]['projected_status'],
            'local_status': phase2[case]['local_status'],
            'transport_status': phase2[case]['transport_status'],
            'ledger_status': phase2[case]['ledger_status'],
            'top_transport_component': summary[case]['top_abs_component'],
            'top_transport_abs_rate': float(summary[case]['top_abs_rate']),
        })

    if not rows:
        raise SystemExit('No Step-217 rows computed.')

    rows.sort(key=lambda r: r['case_name'])
    forced_rows = [r for r in rows if r['case_name'] != 'patched_env_ref']
    best_clean = [r for r in forced_rows if r['transport_status'] == 'PASS']
    best_clean = max(best_clean, key=lambda r: r['post_force_rms_gap']) if best_clean else None
    best_forced = max(forced_rows, key=lambda r: r['post_force_rms_gap'])

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open('w', encoding='utf-8') as f:
        f.write('# Step-217 Post-Patch Environment Burst-Train Sweep\n\n')
        f.write('| case | psi0 ratio | final gap | post-force max gap | post-force rms gap | memory ratio | transport | top transport | ledger |\n')
        f.write('|:---|---:|---:|---:|---:|---:|:---|:---|:---|\n')
        for r in rows:
            f.write(
                f"| {r['case_name']} | {fmt(r['psi0_ratio_full_over_controlled'])} | {fmt(r['final_gap'])} | "
                f"{fmt(r['post_force_max_abs_gap'])} | {fmt(r['post_force_rms_gap'])} | {fmt(r['memory_ratio'])} | "
                f"{r['transport_status']} | {r['top_transport_component']}={fmt(r['top_transport_abs_rate'])} | {r['ledger_status']} |\n"
            )
        f.write('\n## Best Transport-Clean Forced Case\n\n')
        if best_clean is None:
            f.write('- none\n\n')
        else:
            f.write(
                f"- case: `{best_clean['case_name']}`\n"
                f"- post-force rms `psi_w0 - psi_proj`: `{fmt(best_clean['post_force_rms_gap'])}`\n"
                f"- final `psi_w0 - psi_proj`: `{fmt(best_clean['final_gap'])}`\n\n"
            )
        f.write('## Strongest Forced Case\n\n')
        f.write(
            f"- case: `{best_forced['case_name']}`\n"
            f"- post-force rms `psi_w0 - psi_proj`: `{fmt(best_forced['post_force_rms_gap'])}`\n"
            f"- memory ratio: `{fmt(best_forced['memory_ratio'])}`\n"
            f"- transport: `{best_forced['transport_status']}`\n"
            f"- top transport: `{best_forced['top_transport_component']}={fmt(best_forced['top_transport_abs_rate'])}`\n"
        )


if __name__ == '__main__':
    main()
