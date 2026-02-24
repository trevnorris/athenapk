#!/usr/bin/env python3
"""Summarize Step-114 overlap/phase-lock sweep results.

Ranks cases by max |JwEw| and highlights whether overlap ingredients are active.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _f(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def _fmt(x: float) -> str:
    if x != x:
        return "nan"
    return f"{x:.6e}"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--summary-csv", required=True, help="Input summary CSV path")
    p.add_argument("--out-csv", required=True, help="Output ranked CSV path")
    p.add_argument("--out-md", required=True, help="Output ranked markdown path")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    inp = Path(args.summary_csv)
    rows = []
    with inp.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(
                {
                    "case_name": r["case_name"],
                    "max_jw_ew": _f(r.get("full_max_abs_jw_ew", "nan")),
                    "max_jw_l1": _f(r.get("full_max_abs_jw_l1", "nan")),
                    "max_ew_l1": _f(r.get("full_max_abs_ew_l1", "nan")),
                    "max_c2": _f(r.get("full_max_abs_mixed_c2", "nan")),
                    "max_semw": _f(r.get("full_max_abs_em_leak_w", "nan")),
                    "max_emf": _f(r.get("full_max_abs_emf_vw_c_abs", "nan")),
                    "psi_w0_minus_proj": _f(r.get("full_final_psi_w0_minus_psi_proj_span", "nan")),
                    "psi0_ratio": _f(r.get("ratio_full_over_controlled_psi0", "nan")),
                    "gate_jw_ew_active": r.get("gate_jw_ew_active", "N/A"),
                    "gate_jw_active": r.get("gate_jw_active", "N/A"),
                    "gate_ew_active": r.get("gate_ew_active", "N/A"),
                    "gate_c_active": r.get("gate_c_active", "N/A"),
                    "transport": r.get("full_transport_closure_status", "N/A"),
                    "ledger": r.get("full_em_bulk_ledger_status", "N/A"),
                }
            )

    rows.sort(key=lambda r: (r["max_jw_ew"] if r["max_jw_ew"] == r["max_jw_ew"] else -1.0), reverse=True)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case_name",
                "max_jw_ew",
                "max_jw_l1",
                "max_ew_l1",
                "max_c2",
                "max_semw",
                "max_emf",
                "psi_w0_minus_proj",
                "psi0_ratio",
                "gate_jw_ew_active",
                "gate_jw_active",
                "gate_ew_active",
                "gate_c_active",
                "transport",
                "ledger",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    out_md = Path(args.out_md)
    with out_md.open("w") as f:
        f.write("# Step-114 Overlap / Phase-Lock Sweep Summary\n\n")
        f.write(f"Input summary: `{inp}`\n\n")
        f.write("| case | max|JwEw| | max|Jw| | max|Ew| | max|C2| | max|S_EMw| | max|emf_vw_c| | psi_w0-psi_proj | psi0 ratio | gate JwEw/Jw/Ew/C | transport/ledger |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|\n")
        for r in rows:
            f.write(
                "| "
                + f"{r['case_name']} | {_fmt(r['max_jw_ew'])} | {_fmt(r['max_jw_l1'])} | {_fmt(r['max_ew_l1'])} | "
                + f"{_fmt(r['max_c2'])} | {_fmt(r['max_semw'])} | {_fmt(r['max_emf'])} | "
                + f"{_fmt(r['psi_w0_minus_proj'])} | {_fmt(r['psi0_ratio'])} | "
                + f"{r['gate_jw_ew_active']}/{r['gate_jw_active']}/{r['gate_ew_active']}/{r['gate_c_active']} | "
                + f"{r['transport']}/{r['ledger']} |\n"
            )

        if rows:
            best = rows[0]
            f.write("\n## Best Case By max|JwEw|\n\n")
            f.write(f"- case: `{best['case_name']}`\n")
            f.write(f"- max|JwEw|: `{_fmt(best['max_jw_ew'])}`\n")
            f.write(f"- max|Jw|: `{_fmt(best['max_jw_l1'])}`\n")
            f.write(f"- max|Ew|: `{_fmt(best['max_ew_l1'])}`\n")
            f.write(f"- max|C2|: `{_fmt(best['max_c2'])}`\n")


if __name__ == "__main__":
    main()
