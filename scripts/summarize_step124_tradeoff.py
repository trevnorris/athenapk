#!/usr/bin/env python3
"""Combine phase2 + controlled-bounds summaries into a tradeoff view."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List


def _read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _get(row: Dict[str, str], key: str, default: str = "") -> str:
    value = row.get(key, default)
    if value is None:
        return default
    return value


def _status_dual(controlled: str, active: str, projected: str, predictive: str) -> tuple[str, str]:
    relaxed = "PASS" if controlled == "PASS" and active == "PASS" and projected == "PASS" else "FAIL"
    strict = (
        "PASS"
        if controlled == "PASS" and active == "PASS" and projected == "PASS" and predictive == "PASS"
        else "FAIL"
    )
    return relaxed, strict


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase2-csv", required=True)
    parser.add_argument("--controlled-bounds-csv", required=True)
    parser.add_argument("--summary-csv", required=False, default="")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    phase2_rows = _read_csv(Path(args.phase2_csv))
    bounds_rows = _read_csv(Path(args.controlled_bounds_csv))
    summary_rows = _read_csv(Path(args.summary_csv)) if args.summary_csv else []

    phase2_by_case = {_get(r, "case_name"): r for r in phase2_rows}
    bounds_by_case = {_get(r, "case_name"): r for r in bounds_rows}
    summary_by_case = {_get(r, "case_name"): r for r in summary_rows}

    case_names = sorted(set(phase2_by_case) | set(bounds_by_case) | set(summary_by_case))
    out_rows: List[Dict[str, str]] = []
    for case_name in case_names:
        p = phase2_by_case.get(case_name, {})
        b = bounds_by_case.get(case_name, {})
        s = summary_by_case.get(case_name, {})

        controlled_status = _get(b, "status", "UNKNOWN")
        active_transfer = _get(p, "active_transfer_status", "UNKNOWN")
        projected_status = _get(p, "projected_status", "UNKNOWN")
        predictive_status = _get(p, "predictive_status", "UNKNOWN")
        relaxed, strict = _status_dual(
            controlled=controlled_status,
            active=active_transfer,
            projected=projected_status,
            predictive=predictive_status,
        )

        out_rows.append(
            {
                "case_name": case_name,
                "controlled_bounds_status": controlled_status,
                "active_transfer_status": active_transfer,
                "predictive_status": predictive_status,
                "projected_status": projected_status,
                "dual_gate_relaxed_status": relaxed,
                "dual_gate_strict_status": strict,
                "psi0_ratio_full_over_controlled": _get(p, "psi0_ratio_full_over_controlled", _get(s, "ratio_full_over_controlled_psi0")),
                "psi_w0_minus_proj": _get(p, "psi_w0_minus_proj", _get(s, "full_final_psi_w0_minus_psi_proj_span")),
                "max_abs_jw_ew": _get(p, "max_abs_jw_ew", _get(s, "full_max_abs_jw_ew")),
                "max_abs_s_leak_abs": _get(p, "max_abs_s_leak_abs", _get(s, "full_max_abs_s_leak_abs")),
                "max_abs_emf_vw_c_abs": _get(p, "max_abs_emf_vw_c_abs", _get(s, "full_max_abs_emf_vw_c_abs")),
            }
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_name",
        "controlled_bounds_status",
        "active_transfer_status",
        "predictive_status",
        "projected_status",
        "dual_gate_relaxed_status",
        "dual_gate_strict_status",
        "psi0_ratio_full_over_controlled",
        "psi_w0_minus_proj",
        "max_abs_jw_ew",
        "max_abs_s_leak_abs",
        "max_abs_emf_vw_c_abs",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    relaxed_pass = sum(1 for r in out_rows if r["dual_gate_relaxed_status"] == "PASS")
    strict_pass = sum(1 for r in out_rows if r["dual_gate_strict_status"] == "PASS")

    out_md = Path(args.out_md)
    lines: List[str] = []
    lines.append("# Step-124 Tradeoff Summary")
    lines.append("")
    lines.append(f"Cases: `{len(out_rows)}`")
    lines.append(f"- Relaxed dual gate PASS (`controlled + active_transfer + projected`): `{relaxed_pass}/{len(out_rows)}`")
    lines.append(f"- Strict dual gate PASS (`controlled + active_transfer + projected + predictive`): `{strict_pass}/{len(out_rows)}`")
    lines.append("")
    lines.append("| case | controlled | active transfer | predictive | projected | relaxed | strict | psi0 ratio | psi_w0-psi_proj | max |JwEw| | max |S_leak| | max |emf_vw_c| |")
    lines.append("|:---|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|")
    for r in out_rows:
        lines.append(
            f"| {r['case_name']} | {r['controlled_bounds_status']} | {r['active_transfer_status']} | {r['predictive_status']} | {r['projected_status']} | {r['dual_gate_relaxed_status']} | {r['dual_gate_strict_status']} | {r['psi0_ratio_full_over_controlled']} | {r['psi_w0_minus_proj']} | {r['max_abs_jw_ew']} | {r['max_abs_s_leak_abs']} | {r['max_abs_emf_vw_c_abs']} |"
        )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(out_rows)}")
    print(f"relaxed_pass,{relaxed_pass}/{len(out_rows)}")
    print(f"strict_pass,{strict_pass}/{len(out_rows)}")


if __name__ == "__main__":
    main()

