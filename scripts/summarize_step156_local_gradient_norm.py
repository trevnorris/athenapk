#!/usr/bin/env python3
"""Summarize Step-156 local-gradient norm A/B probe."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_harris_scan_matrix(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_case_rows(run_log: Path) -> tuple[dict[str, str], dict[str, str]]:
    lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()
    header = None
    controlled = None
    full = None
    for i, line in enumerate(lines):
        if not line.startswith("case,final_psi0_span,"):
            continue
        header = line
        for j in range(i + 1, min(i + 24, len(lines))):
            if lines[j].startswith("controlled,"):
                controlled = lines[j]
            elif lines[j].startswith("full,"):
                full = lines[j]
        if header and controlled and full:
            break
    if header is None or controlled is None or full is None:
        raise RuntimeError(f"Could not find controlled/full rows in {run_log}")
    keys = next(csv.reader([header]))
    controlled_vals = next(csv.reader([controlled]))
    full_vals = next(csv.reader([full]))
    return dict(zip(keys, controlled_vals)), dict(zip(keys, full_vals))


def max_val(values):
    if not values:
        return math.nan
    return max(values)


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def safe_ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True)
    parser.add_argument("--baseline-case", default="legacy_combined")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        hst = case_dir / "harris_full.out1.hst"
        if not run_log.exists() or not hst.exists():
            continue
        _, full = parse_case_rows(run_log)
        cols = hsm.parse_hst(hst)

        pi0_transport = safe_float(full.get("pi0_transport_max_abs_rate", "nan"))
        piy_transport = safe_float(full.get("piy_transport_max_abs_rate", "nan"))
        piw_transport = safe_float(full.get("piw_transport_max_abs_rate", "nan"))
        transport_values = {"pi0": pi0_transport, "piy": piy_transport, "piw": piw_transport}
        top_channel = max(transport_values, key=lambda key: transport_values[key])

        row = {
            "case_name": case_dir.name,
            "transport_status": full.get("transport_closure_status", ""),
            "ledger_status": full.get("em_bulk_ledger_status", ""),
            "top_pi_channel": top_channel,
            "top_pi_transport_max_abs_rate": transport_values[top_channel],
            "pi0_transport_max_abs_rate": pi0_transport,
            "piy_transport_max_abs_rate": piy_transport,
            "piw_transport_max_abs_rate": piw_transport,
            "ledger_max_abs_rate": safe_float(full.get("em_bulk_ledger_max_abs_rate", "nan")),
            "selected_grad_abs_max": max_val(cols.get("m4d_response_w0_local_grad_abs", [])),
            "quadrature_grad_abs_max": max_val(cols.get("m4d_response_w0_local_grad_quadrature_abs", [])),
            "overlap_signed_abs_max": max(abs(v) for v in cols.get("m4d_response_w0_local_grad_overlap", [math.nan]) if math.isfinite(v)) if cols.get("m4d_response_w0_local_grad_overlap") else math.nan,
            "overlap_abs_max": max_val(cols.get("m4d_response_w0_local_grad_overlap_abs", [])),
            "mix_rate_abs_max": max_val(cols.get("m4d_int_response_w0_local_gradient_mix_rate_abs", [])),
            "mix_pi_abs_max": max_val(cols.get("m4d_int_response_w0_local_gradient_mix_pi_abs", [])),
        }
        rows.append(row)

    baseline = next((row for row in rows if row["case_name"] == args.baseline_case), None)
    for row in rows:
        if baseline is None:
            row["top_pi_transport_over_baseline"] = math.nan
            row["ledger_over_baseline"] = math.nan
            row["mix_rate_over_baseline"] = math.nan
        else:
            row["top_pi_transport_over_baseline"] = safe_ratio(
                float(row["top_pi_transport_max_abs_rate"]),
                float(baseline["top_pi_transport_max_abs_rate"]),
            )
            row["ledger_over_baseline"] = safe_ratio(
                float(row["ledger_max_abs_rate"]),
                float(baseline["ledger_max_abs_rate"]),
            )
            row["mix_rate_over_baseline"] = safe_ratio(
                float(row["mix_rate_abs_max"]),
                float(baseline["mix_rate_abs_max"]),
            )

    fieldnames = [
        "case_name",
        "transport_status",
        "ledger_status",
        "top_pi_channel",
        "top_pi_transport_max_abs_rate",
        "top_pi_transport_over_baseline",
        "pi0_transport_max_abs_rate",
        "piy_transport_max_abs_rate",
        "piw_transport_max_abs_rate",
        "ledger_max_abs_rate",
        "ledger_over_baseline",
        "selected_grad_abs_max",
        "quadrature_grad_abs_max",
        "overlap_signed_abs_max",
        "overlap_abs_max",
        "mix_rate_abs_max",
        "mix_rate_over_baseline",
        "mix_pi_abs_max",
    ]

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    lines = [
        "# Step-156 Local-Gradient Norm A/B Summary",
        "",
        "This compares alternative local-gradient norm forms for the `w0_local_gradient` operator.",
        "",
        "## Case Summary",
        "",
        "| Case | Transport | Ledger | Top pi | Top pi max abs rate | Top pi / baseline | Ledger / baseline | Selected grad max | Quadrature grad max | Overlap abs max | Mix-rate max | Mix-rate / baseline |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {case_name} | {transport_status} | {ledger_status} | {top_pi_channel} | {top_pi_transport_max_abs_rate} | {top_pi_transport_over_baseline} | {ledger_over_baseline} | {selected_grad_abs_max} | {quadrature_grad_abs_max} | {overlap_abs_max} | {mix_rate_abs_max} | {mix_rate_over_baseline} |".format(
                case_name=row["case_name"],
                transport_status=row["transport_status"],
                ledger_status=row["ledger_status"],
                top_pi_channel=row["top_pi_channel"],
                top_pi_transport_max_abs_rate=fmt(float(row["top_pi_transport_max_abs_rate"])),
                top_pi_transport_over_baseline=fmt(float(row["top_pi_transport_over_baseline"])),
                ledger_over_baseline=fmt(float(row["ledger_over_baseline"])),
                selected_grad_abs_max=fmt(float(row["selected_grad_abs_max"])),
                quadrature_grad_abs_max=fmt(float(row["quadrature_grad_abs_max"])),
                overlap_abs_max=fmt(float(row["overlap_abs_max"])),
                mix_rate_abs_max=fmt(float(row["mix_rate_abs_max"])),
                mix_rate_over_baseline=fmt(float(row["mix_rate_over_baseline"])),
            )
        )

    lines.append("")
    lines.append("## Channel Detail")
    lines.append("")
    lines.append("| Case | pi0 transport | piy transport | piw transport | overlap |")
    lines.append("| --- | ---: | ---: | ---: | ---: |")
    for row in rows:
        lines.append(
            "| {case_name} | {pi0_transport_max_abs_rate} | {piy_transport_max_abs_rate} | {piw_transport_max_abs_rate} | {overlap_signed_abs_max} |".format(
                case_name=row["case_name"],
                pi0_transport_max_abs_rate=fmt(float(row["pi0_transport_max_abs_rate"])),
                piy_transport_max_abs_rate=fmt(float(row["piy_transport_max_abs_rate"])),
                piw_transport_max_abs_rate=fmt(float(row["piw_transport_max_abs_rate"])),
                overlap_signed_abs_max=fmt(float(row["overlap_signed_abs_max"])),
            )
        )

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
