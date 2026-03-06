#!/usr/bin/env python3
"""Decompose transport-closure residuals by mode-0 conserved component."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


COMPONENTS = [
    ("momx", "m4d_momx_mode_0", "m4d_int_divmomx_mode0", "m4d_int_srcmomx_mode0"),
    ("momy", "m4d_momy_mode_0", "m4d_int_divmomy_mode0", "m4d_int_srcmomy_mode0"),
    ("momz", "m4d_momz_mode_0", "m4d_int_divmomz_mode0", "m4d_int_srcmomz_mode0"),
    ("momw", "m4d_momw_mode_0", "m4d_int_divmomw_mode0", "m4d_int_srcmomw_mode0"),
    ("energy", "m4d_energy_mode_0", "m4d_int_divenergy_mode0", "m4d_int_srcenergy_mode0"),
    ("pi0", "m4d_pi0_mode_0", "m4d_int_divpi0_mode0", "m4d_int_srcpi0_mode0"),
    ("pix", "m4d_pix_mode_0", "m4d_int_divpix_mode0", "m4d_int_srcpix_mode0"),
    ("piy", "m4d_piy_mode_0", "m4d_int_divpiy_mode0", "m4d_int_srcpiy_mode0"),
    ("piz", "m4d_piz_mode_0", "m4d_int_divpiz_mode0", "m4d_int_srcpiz_mode0"),
    ("piw", "m4d_piw_mode_0", "m4d_int_divpiw_mode0", "m4d_int_srcpiw_mode0"),
]


def load_hsm(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def maybe_col(cols, key):
    values = cols.get(key)
    if values is None:
        return None
    return values


def fmt(value: float) -> str:
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.6e}"


def status_from_metrics(max_norm: float, max_abs: float, norm_tol: float, abs_tol: float) -> str:
    if not math.isfinite(max_norm):
        return "N/A"
    return "PASS" if (max_norm <= norm_tol or max_abs <= abs_tol) else "FAIL"


def format_top(entries, key: str, n: int = 3) -> str:
    ranked = [row for row in entries if math.isfinite(float(row[key]))]
    if not ranked:
        return "none"
    ranked.sort(key=lambda row: float(row[key]), reverse=True)
    parts = []
    for row in ranked[:n]:
        parts.append(f"{row['component']}={fmt(row[key])}")
    return ", ".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Absolute or relative case-dir glob")
    parser.add_argument("--out-csv", required=True, help="Per-component CSV output path")
    parser.add_argument("--out-md", required=True, help="Markdown summary output path")
    parser.add_argument("--summary-csv", required=True, help="Per-case summary CSV output path")
    parser.add_argument("--baseline-case", default="baseline_decomp")
    parser.add_argument("--closure-norm-tol", type=float, default=5.0e-2)
    parser.add_argument("--closure-abs-rate-tol", type=float, default=1.0e-8)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]

    component_rows: list[dict[str, object]] = []
    case_summaries: list[dict[str, object]] = []

    for case_dir in case_dirs:
        hst = case_dir / "harris_full.out1.hst"
        if not hst.exists():
            continue
        cols = hsm.parse_hst(hst)
        times = cols.get("time", [])
        case_entries: list[dict[str, object]] = []
        for component, quantity_key, div_key, src_key in COMPONENTS:
            quantity = maybe_col(cols, quantity_key)
            div = maybe_col(cols, div_key)
            src = maybe_col(cols, src_key)
            max_norm, rms_norm, max_abs_rate, final_rate = hsm.transport_balance_metrics(
                times,
                quantity,
                div,
                src,
                args.closure_abs_rate_tol,
            )
            row = {
                "case_name": case_dir.name,
                "component": component,
                "max_norm": max_norm,
                "rms_norm": rms_norm,
                "max_abs_rate": max_abs_rate,
                "final_rate": final_rate,
                "status": status_from_metrics(
                    max_norm, max_abs_rate, args.closure_norm_tol, args.closure_abs_rate_tol
                ),
            }
            component_rows.append(row)
            case_entries.append(row)

        top_abs = max(
            case_entries,
            key=lambda row: float(row["max_abs_rate"]) if math.isfinite(float(row["max_abs_rate"])) else -1.0,
        )
        top_norm = max(
            case_entries,
            key=lambda row: float(row["max_norm"]) if math.isfinite(float(row["max_norm"])) else -1.0,
        )
        case_summaries.append(
            {
                "case_name": case_dir.name,
                "top_abs_component": top_abs["component"],
                "top_abs_rate": top_abs["max_abs_rate"],
                "top_abs_status": top_abs["status"],
                "top_norm_component": top_norm["component"],
                "top_norm_value": top_norm["max_norm"],
                "top_norm_status": top_norm["status"],
                "top3_abs": format_top(case_entries, "max_abs_rate"),
                "top3_norm": format_top(case_entries, "max_norm"),
            }
        )

    baseline_top_abs = math.nan
    baseline_top_norm = math.nan
    for row in case_summaries:
        if row["case_name"] == args.baseline_case:
            baseline_top_abs = float(row["top_abs_rate"])
            baseline_top_norm = float(row["top_norm_value"])
            break

    for row in case_summaries:
        top_abs_rate = float(row["top_abs_rate"])
        top_norm_value = float(row["top_norm_value"])
        row["top_abs_reduction_vs_baseline"] = (
            baseline_top_abs / top_abs_rate
            if math.isfinite(baseline_top_abs)
            and math.isfinite(top_abs_rate)
            and top_abs_rate != 0.0
            else math.nan
        )
        row["top_norm_reduction_vs_baseline"] = (
            baseline_top_norm / top_norm_value
            if math.isfinite(baseline_top_norm)
            and math.isfinite(top_norm_value)
            and top_norm_value != 0.0
            else math.nan
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "case_name",
            "component",
            "max_norm",
            "rms_norm",
            "max_abs_rate",
            "final_rate",
            "status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in component_rows:
            writer.writerow(row)

    summary_csv = Path(args.summary_csv)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "case_name",
            "top_abs_component",
            "top_abs_rate",
            "top_abs_status",
            "top_norm_component",
            "top_norm_value",
            "top_norm_status",
            "top_abs_reduction_vs_baseline",
            "top_norm_reduction_vs_baseline",
            "top3_abs",
            "top3_norm",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in case_summaries:
            writer.writerow(row)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as handle:
        handle.write("# Step-149 Transport Component Decomposition\n\n")
        handle.write(f"Parsed `{len(case_summaries)}` cases from `{args.glob}`.\n\n")
        handle.write(
            "Transport component metrics are recomputed directly from the full-case `.hst` series "
            "using `transport_balance_metrics` on each mode-0 conserved variable.\n\n"
        )
        if case_summaries:
            handle.write(
                "| case | top abs-rate component | top abs rate | abs status | top norm component | top norm | norm status | abs reduction vs baseline | top-3 abs | top-3 norm |\n"
            )
            handle.write(
                "|:---|:---|---:|:---|:---|---:|:---|---:|:---|:---|\n"
            )
            for row in case_summaries:
                handle.write(
                    f"| {row['case_name']} | {row['top_abs_component']} | {fmt(row['top_abs_rate'])} | "
                    f"{row['top_abs_status']} | {row['top_norm_component']} | {fmt(row['top_norm_value'])} | "
                    f"{row['top_norm_status']} | {fmt(row['top_abs_reduction_vs_baseline'])} | "
                    f"{row['top3_abs']} | {row['top3_norm']} |\n"
                )

            baseline_rows = [row for row in component_rows if row["case_name"] == args.baseline_case]
            if baseline_rows:
                ranked_baseline = sorted(
                    baseline_rows,
                    key=lambda row: float(row["max_abs_rate"])
                    if math.isfinite(float(row["max_abs_rate"]))
                    else -1.0,
                    reverse=True,
                )
                handle.write("\n## Baseline Ranking\n\n")
                handle.write("| component | max norm | rms norm | max abs rate | final rate | status |\n")
                handle.write("|:---|---:|---:|---:|---:|:---|\n")
                for row in ranked_baseline:
                    handle.write(
                        f"| {row['component']} | {fmt(row['max_norm'])} | {fmt(row['rms_norm'])} | "
                        f"{fmt(row['max_abs_rate'])} | {fmt(row['final_rate'])} | {row['status']} |\n"
                    )


if __name__ == "__main__":
    main()
