#!/usr/bin/env python3
"""Summarize Step-121 Phase-2 predictive transfer checks."""

from __future__ import annotations

import argparse
import glob
import importlib.util
import math
from pathlib import Path
from typing import Dict, List


def load_hsm(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def series_max_abs(series):
    if series is None or len(series) == 0:
        return math.nan
    try:
        return max(abs(float(v)) for v in series)
    except (TypeError, ValueError):
        return math.nan


def safe_ratio(numer, denom):
    if (
        numer is None
        or denom is None
        or not math.isfinite(numer)
        or not math.isfinite(denom)
        or abs(denom) == 0.0
    ):
        return math.nan
    return numer / denom


def fmt(value):
    if value is None:
        return "nan"
    if isinstance(value, str):
        return value
    if isinstance(value, float) and not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def corr_pass(v: float, threshold: float) -> bool:
    return math.isfinite(v) and abs(v) >= threshold


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Glob for case dirs")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-md", required=True)
    p.add_argument("--jw-ew-floor", type=float, default=1.0e-12)
    p.add_argument("--s-leak-floor", type=float, default=1.0e-12)
    p.add_argument("--corr-threshold", type=float, default=5.0e-1)
    args = p.parse_args()

    dirs = [Path(pth) for pth in sorted(glob.glob(args.glob))]
    if not dirs:
        raise SystemExit(f"No directories matched: {args.glob}")

    hsm = load_hsm(Path(__file__).resolve().parent)

    rows: List[Dict[str, object]] = []
    for d in dirs:
        controlled_hst = d / "harris_controlled.out1.hst"
        full_hst = d / "harris_full.out1.hst"
        if not controlled_hst.exists() or not full_hst.exists():
            continue

        controlled_cols = hsm.parse_hst(controlled_hst)
        full_cols = hsm.parse_hst(full_hst)

        controlled = hsm.analyze_case(
            "controlled",
            controlled_cols,
            5.0e-2,
            1.0e-8,
            1.0e-8,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )
        full = hsm.analyze_case(
            "full",
            full_cols,
            5.0e-2,
            1.0e-8,
            1.0e-8,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )

        max_abs_jw_ew = series_max_abs(full_cols.get("m4d_int_jw_ew"))
        max_abs_s_leak_abs = series_max_abs(full_cols.get("m4d_int_s_leak_abs"))
        max_abs_emf_vw_c_abs = series_max_abs(full_cols.get("m4d_emf_vw_c_abs"))
        max_abs_brane_epar2 = series_max_abs(full_cols.get("m4d_brane_epar2"))

        corr_names = [
            "corr_psi0_jw_ew",
            "corr_psi0_s_leak_abs",
            "corr_psi0_emf_vw_c_abs",
            "corr_psiw0_jw_ew",
            "corr_psiw0_s_leak_abs",
            "corr_psiw0_emf_vw_c_abs",
        ]
        corr_vals = [full.get(name, math.nan) for name in corr_names]
        corr_hits = sum(corr_pass(v, args.corr_threshold) for v in corr_vals)
        finite_corrs = sum(math.isfinite(v) for v in corr_vals)

        active_transfer = (
            math.isfinite(max_abs_jw_ew)
            and math.isfinite(max_abs_s_leak_abs)
            and max_abs_jw_ew >= args.jw_ew_floor
            and max_abs_s_leak_abs >= args.s_leak_floor
        )
        predictive_status = (
            "PASS"
            if active_transfer and corr_hits >= 2
            else ("PARTIAL" if active_transfer else "FAIL")
        )

        rows.append(
            {
                "case_name": d.name,
                "psi0_ratio_full_over_controlled": safe_ratio(
                    full.get("final_psi0_span", math.nan),
                    controlled.get("final_psi0_span", math.nan),
                ),
                "psi_w0_minus_proj": full.get("final_psi_w0_minus_psi_proj_span", math.nan),
                "max_abs_jw_ew": max_abs_jw_ew,
                "max_abs_s_leak_abs": max_abs_s_leak_abs,
                "max_abs_emf_vw_c_abs": max_abs_emf_vw_c_abs,
                "max_abs_brane_epar2": max_abs_brane_epar2,
                "corr_psi0_jw_ew": full.get("corr_psi0_jw_ew", math.nan),
                "corr_psi0_s_leak_abs": full.get("corr_psi0_s_leak_abs", math.nan),
                "corr_psi0_emf_vw_c_abs": full.get("corr_psi0_emf_vw_c_abs", math.nan),
                "corr_psiw0_jw_ew": full.get("corr_psiw0_jw_ew", math.nan),
                "corr_psiw0_s_leak_abs": full.get("corr_psiw0_s_leak_abs", math.nan),
                "corr_psiw0_emf_vw_c_abs": full.get("corr_psiw0_emf_vw_c_abs", math.nan),
                "corr_hits_ge_threshold": corr_hits,
                "corr_finite_count": finite_corrs,
                "active_transfer_status": "PASS" if active_transfer else "FAIL",
                "predictive_status": predictive_status,
                "projected_status": full.get("closure_projected_status", "N/A"),
                "local_status": full.get("closure_local_mode0_status", "N/A"),
                "transport_status": full.get("transport_closure_status", "N/A"),
                "ledger_status": full.get("em_bulk_ledger_status", "N/A"),
            }
        )

    rows.sort(
        key=lambda r: (
            {"PASS": 0, "PARTIAL": 1, "FAIL": 2}.get(str(r["predictive_status"]), 9),
            -abs(float(r["max_abs_jw_ew"])) if math.isfinite(float(r["max_abs_jw_ew"])) else -1.0,
        )
    )

    out_csv = Path(args.out_csv)
    out_md = Path(args.out_md)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    import csv

    with out_csv.open("w", newline="") as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-121 Phase-2 Predictive Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        f.write("Gates:\n")
        f.write(f"- active transfer: `max|JwEw| >= {args.jw_ew_floor:.3e}` and `max|S_leak_abs| >= {args.s_leak_floor:.3e}`\n")
        f.write(f"- predictive correlation threshold: `|corr| >= {args.corr_threshold:.3f}` (need >=2 hits)\n\n")
        if not rows:
            f.write("No rows parsed.\n")
            return
        headers = list(rows[0].keys())
        f.write("| " + " | ".join(h.replace("_", " ") for h in headers) + " |\n")
        f.write("|" + "|".join([":---"] * len(headers)) + "|\n")
        for r in rows:
            f.write("| " + " | ".join(fmt(r[h]) for h in headers) + " |\n")


if __name__ == "__main__":
    main()
