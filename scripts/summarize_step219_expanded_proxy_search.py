#!/usr/bin/env python3
"""Summarize an expanded proxy search on the Step-218 A/B outputs."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path

REF_CASE = 'patched_env_ref'
FORCED_CASE = 'patched_env_train_quad_amp0010'
PROXIES = {
    'R_gap': lambda cols, eps, sigma: [a - b for a, b in zip(cols['m4d_psi_w0_span'], cols['m4d_psi_proj_span'])],
    'Lz_mixed_cz2': lambda cols, eps, sigma: cols['m4d_mixed_cz2'],
    'Ltot_mixed_c2': lambda cols, eps, sigma: cols['m4d_mixed_c2'],
    'Lew_mixed_ew2': lambda cols, eps, sigma: cols['m4d_mixed_ew2'],
    'Bz_zblock_correction': lambda cols, eps, sigma: cols['m4d_em_u_bulk_zblock_correction'],
    'Qz_partition': lambda cols, eps, sigma: [bz / (lz + eps) for bz, lz in zip(cols['m4d_em_u_bulk_zblock_correction'], cols['m4d_mixed_cz2'])],
    'Chi_sigma_star': lambda cols, eps, sigma: [bz / sigma for bz in cols['m4d_em_u_bulk_zblock_correction']],
    'JwEw_abs': lambda cols, eps, sigma: [abs(v) for v in cols['m4d_int_jw_ew']],
    'S_leak_abs': lambda cols, eps, sigma: cols['m4d_int_s_leak_abs'],
    'EMF_vw_c_abs': lambda cols, eps, sigma: cols['m4d_emf_vw_c_abs'],
    'EMF_cov_vxb_abs': lambda cols, eps, sigma: cols['m4d_emf_cov_vxb_abs'],
    'Src_mode0_total_abs': lambda cols, eps, sigma: cols['m4d_int_src_mode0_total_abs'],
    'Bridge_power_mode0_abs': lambda cols, eps, sigma: cols['m4d_int_bridge_power_mode0_abs'],
    'Geometry_shift_power_mode0_abs': lambda cols, eps, sigma: cols['m4d_int_geometry_shift_power_mode0_abs'],
    'Psi0_span': lambda cols, eps, sigma: cols['m4d_psi0_span'],
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


def window_values(times, values, t0, t1=None):
    return [v for t, v in zip(times, values) if t >= t0 and (t1 is None or t < t1) and math.isfinite(v)]


def mean_abs(vals):
    if not vals:
        return math.nan
    return sum(abs(v) for v in vals) / len(vals)


def rms(vals):
    if not vals:
        return math.nan
    return math.sqrt(sum(v * v for v in vals) / len(vals))


def max_abs(vals):
    if not vals:
        return math.nan
    return max(abs(v) for v in vals)


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
    ap.add_argument('--glob', required=True)
    ap.add_argument('--summary-csv', required=True)
    ap.add_argument('--phase2-csv', required=True)
    ap.add_argument('--gate-on', type=float, required=True)
    ap.add_argument('--gate-off', type=float, required=True)
    ap.add_argument('--sigma-star', type=float, default=2.07614e-01)
    ap.add_argument('--eps', type=float, default=1.0e-30)
    ap.add_argument('--out-csv', required=True)
    ap.add_argument('--out-md', required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    summary = {r['case_name']: r for r in read_csv(Path(args.summary_csv))}
    phase2 = {r['case_name']: r for r in read_csv(Path(args.phase2_csv))}

    cases = {}
    for case_dir in [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]:
        case = case_dir.name
        full_hst = case_dir / 'harris_full.out1.hst'
        if full_hst.exists() and case in (REF_CASE, FORCED_CASE):
            cols = hsm.parse_hst(full_hst)
            times = cols['time']
            proxies = {name: fn(cols, args.eps, args.sigma_star) for name, fn in PROXIES.items()}
            cases[case] = {'times': times, 'proxies': proxies}

    missing = [case for case in (REF_CASE, FORCED_CASE) if case not in cases]
    if missing:
        raise SystemExit(f'Missing required cases: {missing}')

    ref = cases[REF_CASE]
    forced = cases[FORCED_CASE]

    rows = []
    for case_name, case_data in cases.items():
        times = case_data['times']
        for proxy_name, values in case_data['proxies'].items():
            rows.append({
                'case_name': case_name,
                'proxy': proxy_name,
                'pre_mean_abs': mean_abs(window_values(times, values, 0.0, args.gate_on)),
                'load_max_abs': max_abs(window_values(times, values, args.gate_on, args.gate_off)),
                'post_max_abs': max_abs(window_values(times, values, args.gate_off, None)),
                'post_rms': rms(window_values(times, values, args.gate_off, None)),
            })

    ref_envelopes = {
        proxy_name: max_abs(window_values(ref['times'], values, args.gate_on, None))
        for proxy_name, values in ref['proxies'].items()
    }
    forced_cross = {
        proxy_name: first_cross_time(forced['times'], values, ref_envelopes[proxy_name], args.gate_on)
        for proxy_name, values in forced['proxies'].items()
    }
    t_release = forced_cross['R_gap']
    snapshot = {
        proxy_name: sample_at_or_after(forced['times'], values, t_release)
        for proxy_name, values in forced['proxies'].items()
    }

    ranked = []
    for proxy_name, t_cross in forced_cross.items():
        if not math.isfinite(t_cross):
            continue
        lead = t_release - t_cross if math.isfinite(t_release) else math.nan
        ranked.append((proxy_name, t_cross, lead))
    ranked.sort(key=lambda item: item[1])

    candidate_rows = [row for row in ranked if row[0] != 'R_gap']
    pre_release = [row for row in candidate_rows if math.isfinite(row[2]) and row[2] >= 0.0]

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['case_name', 'proxy', 'pre_mean_abs', 'load_max_abs', 'post_max_abs', 'post_rms'])
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open('w', encoding='utf-8') as f:
        f.write('# Step-219 Expanded Proxy Search\n\n')
        f.write('This extends the Step-218 A/B with leakage, topology, and low-order mode-0 proxy channels already present in the HST output.\n\n')
        f.write('## Case Summary\n\n')
        f.write('| case | psi0 ratio | projected | local | transport | ledger | top transport |\n')
        f.write('|:---|---:|:---|:---|:---|:---|:---|\n')
        for case_name in (REF_CASE, FORCED_CASE):
            f.write(
                f"| {case_name} | {fmt(phase2[case_name]['psi0_ratio_full_over_controlled'])} | {phase2[case_name]['projected_status']} | {phase2[case_name]['local_status']} | {phase2[case_name]['transport_status']} | {phase2[case_name]['ledger_status']} | {summary[case_name]['top_abs_component']}={fmt(summary[case_name]['top_abs_rate'])} |\n"
            )
        f.write('\n## Proxy Windows\n\n')
        f.write('| case | proxy | pre mean abs | load max abs | post max abs | post rms |\n')
        f.write('|:---|:---|---:|---:|---:|---:|\n')
        for row in rows:
            f.write(
                f"| {row['case_name']} | {row['proxy']} | {fmt(row['pre_mean_abs'])} | {fmt(row['load_max_abs'])} | {fmt(row['post_max_abs'])} | {fmt(row['post_rms'])} |\n"
            )
        f.write('\n## Forced-Case Crossing Times\n\n')
        for proxy_name, t_cross, lead in ranked:
            f.write(f'- `{proxy_name}` crosses at `t = {fmt(t_cross)}` with lead `t_R - t = {fmt(lead)}`\n')
        missing_cross = [name for name, t in forced_cross.items() if not math.isfinite(t)]
        if missing_cross:
            f.write('\nNo crossing detected for: ' + ', '.join(f'`{name}`' for name in missing_cross) + '\n')
        f.write('\n## Snapshot At Release Onset\n\n')
        f.write(f'- release onset `t_R = {fmt(t_release)}` from `R_gap`\n')
        for proxy_name in (
            'Lz_mixed_cz2', 'Ltot_mixed_c2', 'Lew_mixed_ew2', 'Bz_zblock_correction',
            'JwEw_abs', 'S_leak_abs', 'EMF_vw_c_abs', 'EMF_cov_vxb_abs',
            'Src_mode0_total_abs', 'Bridge_power_mode0_abs', 'Geometry_shift_power_mode0_abs',
            'Psi0_span', 'Qz_partition', 'Chi_sigma_star',
        ):
            f.write(f'- `{proxy_name}(t_R)` = `{fmt(snapshot[proxy_name])}`\n')
        f.write('\n## Candidate Assessment\n\n')
        if pre_release:
            for proxy_name, t_cross, lead in pre_release:
                f.write(f'- `{proxy_name}` crosses no later than release onset with lead `{fmt(lead)}`\n')
        else:
            f.write('- No non-`R_gap` proxy crosses before or at release onset.\n')
        f.write('\n## Proxy Families\n\n')
        f.write('- support loading: `Lz_mixed_cz2`, `Ltot_mixed_c2`, `Lew_mixed_ew2`\n')
        f.write('- geometry/partition: `Bz_zblock_correction`, `Qz_partition`, `Chi_sigma_star`\n')
        f.write('- leakage/work: `JwEw_abs`, `S_leak_abs`\n')
        f.write('- topology EMF: `EMF_vw_c_abs`, `EMF_cov_vxb_abs`\n')
        f.write('- low-order mode-0: `Src_mode0_total_abs`, `Bridge_power_mode0_abs`, `Geometry_shift_power_mode0_abs`, `Psi0_span`\n')


if __name__ == '__main__':
    main()
