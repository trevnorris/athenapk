#!/usr/bin/env python3
"""Summarize Step-119 transverse-channel pinning audit."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary-csv", required=True, help="Input summary CSV from summarize_transverse_channel_activation.py")
    p.add_argument("--out-csv", required=True, help="Output focused pinning CSV")
    p.add_argument("--out-md", required=True, help="Output focused pinning markdown")
    p.add_argument("--jw-floor", type=float, default=1.0e-2)
    p.add_argument("--ew-floor", type=float, default=2.0e-2)
    p.add_argument("--c-floor", type=float, default=1.0e-1)
    p.add_argument("--jw-ew-floor", type=float, default=1.0e-12)
    p.add_argument("--s-leak-floor", type=float, default=1.0e-12)
    return p.parse_args()


def as_float(row: Dict[str, str], key: str) -> float:
    raw = row.get(key, "")
    try:
        v = float(raw)
    except ValueError:
        return float("nan")
    return v


def gate(value: float, floor: float) -> bool:
    return math.isfinite(value) and abs(value) >= floor


def classify(row: Dict[str, str], args: argparse.Namespace) -> str:
    g_jw = gate(as_float(row, "full_max_abs_jw_l1"), args.jw_floor)
    g_ew = gate(as_float(row, "full_max_abs_ew_l1"), args.ew_floor)
    g_c = gate(as_float(row, "full_max_abs_mixed_c2"), args.c_floor)
    g_work = gate(as_float(row, "full_max_abs_jw_ew"), args.jw_ew_floor)
    g_leak = gate(as_float(row, "full_max_abs_s_leak_abs"), args.s_leak_floor)

    if g_jw and g_ew and g_c and g_work and g_leak:
        return "TRANSVERSE_ACTIVE"
    if not (g_jw or g_ew or g_c):
        return "PINNED_CONTROLLED_LIMIT"
    return "PARTIAL_ACTIVATION"


def fmt(x: float) -> str:
    if not math.isfinite(x):
        return "nan"
    return f"{x:.6e}"


def main() -> None:
    args = parse_args()
    in_path = Path(args.summary_csv)
    rows = list(csv.DictReader(in_path.open()))

    out_rows: List[Dict[str, str]] = []
    for row in rows:
        out_rows.append(
            {
                "case_name": row["case_name"],
                "classification": classify(row, args),
                "max_abs_jw_l1": fmt(as_float(row, "full_max_abs_jw_l1")),
                "max_abs_ew_l1": fmt(as_float(row, "full_max_abs_ew_l1")),
                "max_abs_mixed_c2": fmt(as_float(row, "full_max_abs_mixed_c2")),
                "max_abs_jw_ew": fmt(as_float(row, "full_max_abs_jw_ew")),
                "max_abs_s_leak_abs": fmt(as_float(row, "full_max_abs_s_leak_abs")),
                "max_abs_em_leak_w": fmt(as_float(row, "full_max_abs_em_leak_w")),
                "max_abs_emf_vw_c_abs": fmt(as_float(row, "full_max_abs_emf_vw_c_abs")),
                "psi0_ratio": fmt(as_float(row, "ratio_full_over_controlled_psi0")),
                "psi_w0_minus_proj": fmt(as_float(row, "full_final_psi_w0_minus_psi_proj_span")),
                "projected_status": row.get("full_closure_projected_status", "N/A"),
                "local_status": row.get("full_closure_local_mode0_status", "N/A"),
                "transport_status": row.get("full_transport_closure_status", "N/A"),
                "ledger_status": row.get("full_em_bulk_ledger_status", "N/A"),
            }
        )

    class_order = {"TRANSVERSE_ACTIVE": 0, "PARTIAL_ACTIVATION": 1, "PINNED_CONTROLLED_LIMIT": 2}
    out_rows.sort(
        key=lambda r: (
            class_order.get(r["classification"], 99),
            -abs(float(r["max_abs_jw_ew"])) if r["max_abs_jw_ew"] != "nan" else float("inf"),
        )
    )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()) if out_rows else [])
        if out_rows:
            writer.writeheader()
            writer.writerows(out_rows)

    out_md = Path(args.out_md)
    with out_md.open("w") as f:
        f.write("# Step-119 Transverse Pinning Summary\n\n")
        f.write(f"Parsed `{len(out_rows)}` cases from `{in_path}`.\n\n")
        f.write("Activation floors:\n")
        f.write(f"- `|Jw_l1|max >= {args.jw_floor:.3e}`\n")
        f.write(f"- `|Ew_l1|max >= {args.ew_floor:.3e}`\n")
        f.write(f"- `|mixed_c2|max >= {args.c_floor:.3e}`\n")
        f.write(f"- `|JwEw|max >= {args.jw_ew_floor:.3e}`\n")
        f.write(f"- `|S_leak_abs|max >= {args.s_leak_floor:.3e}`\n\n")
        if not out_rows:
            f.write("No rows parsed.\n")
            return
        headers = list(out_rows[0].keys())
        f.write("| " + " | ".join(h.replace("_", " ") for h in headers) + " |\n")
        f.write("|" + "|".join([":---"] * len(headers)) + "|\n")
        for r in out_rows:
            f.write("| " + " | ".join(r[h] for h in headers) + " |\n")


if __name__ == "__main__":
    main()
