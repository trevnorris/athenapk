#!/usr/bin/env python3
"""Summarize A_w midpoint-pi A/B results."""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path

def read_csv(path: Path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def fmt(v):
    try:
        return f"{float(v):.6e}"
    except Exception:
        return str(v)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summary-csv', required=True)
    ap.add_argument('--transport-csv', required=True)
    ap.add_argument('--out-csv', required=True)
    ap.add_argument('--out-md', required=True)
    args=ap.parse_args()

    summary = {r['case']: r for r in read_csv(Path(args.summary_csv))}
    transport = {r['case']: r for r in read_csv(Path(args.transport_csv))}
    cases=['baseline_quadrature','aw_midpoint_pi']
    rows=[]
    for case in cases:
        s=summary.get(case, {})
        t=transport.get(case, {})
        rows.append({
            'case': case,
            'psi0_ratio': s.get('ratio_full_over_controlled_psi0','nan'),
            'psiw0_minus_proj': s.get('full_final_psi_w0_minus_psi_proj_span','nan'),
            'ledger_max_abs_rate': s.get('full_em_bulk_ledger_max_abs_rate','nan'),
            'ledger_status': s.get('full_em_bulk_ledger_status','nan'),
            'projected_status': s.get('full_closure_status','nan'),
            'local_status': s.get('full_closure_local_mode0_status','nan'),
            'transport_status': s.get('full_transport_status','nan'),
            'piw_transport_abs_rate': t.get('top_abs_rate','nan'),
        })
    out_csv=Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open('w', newline='', encoding='utf-8') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)
    out_md=Path(args.out_md)
    with out_md.open('w', encoding='utf-8') as f:
        f.write('# Step-195 A_w Midpoint-Pi A/B\n\n')
        f.write('| case | psi0 ratio | psi_w0-psi_proj | ledger max abs rate | ledger | projected/local/transport | top transport abs rate |\n')
        f.write('|:---|---:|---:|---:|:---|:---|---:|\n')
        for r in rows:
            f.write(f"| {r['case']} | {fmt(r['psi0_ratio'])} | {fmt(r['psiw0_minus_proj'])} | {fmt(r['ledger_max_abs_rate'])} | {r['ledger_status']} | {r['projected_status']}/{r['local_status']}/{r['transport_status']} | {fmt(r['piw_transport_abs_rate'])} |\n")

if __name__=='__main__':
    main()
