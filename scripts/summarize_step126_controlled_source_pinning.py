#!/usr/bin/env python3
"""Summarize controlled-source pinning audit from per-case run logs."""

from __future__ import annotations

import argparse
import csv
import glob
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _to_float(value: Optional[str]) -> float:
    if value is None:
        return float("nan")
    text = value.strip()
    if text == "":
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def _fmt(value: float) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def _safe_ratio(numer: float, denom: float) -> float:
    if not math.isfinite(numer) or not math.isfinite(denom) or denom == 0.0:
        return float("nan")
    return numer / denom


def _load_case_rows(run_log: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    text = run_log.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    header = None
    controlled = None
    full = None
    for idx, line in enumerate(lines):
        if not line.startswith("case,final_psi0_span,"):
            continue
        header = line
        for j in range(idx + 1, min(idx + 12, len(lines))):
            if lines[j].startswith("controlled,"):
                controlled = lines[j]
            elif lines[j].startswith("full,"):
                full = lines[j]
        if header and controlled and full:
            break

    if header is None or controlled is None or full is None:
        raise RuntimeError(f"Could not find controlled/full CSV rows in {run_log}")

    keys = next(csv.reader([header]))
    c_vals = next(csv.reader([controlled]))
    f_vals = next(csv.reader([full]))
    return dict(zip(keys, c_vals)), dict(zip(keys, f_vals))


def _read_bounds(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {row.get("case_name", ""): row for row in rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Case directory glob (e.g., /path/out/*)")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--out-md", required=True, help="Output Markdown path")
    parser.add_argument(
        "--controlled-bounds-csv",
        default="",
        help="Optional Step-123 bounds CSV to attach status labels",
    )
    parser.add_argument(
        "--baseline-case",
        default="ctrl_baseline",
        help="Case name used as baseline for reduction factors",
    )
    args = parser.parse_args()

    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
    bounds_by_case = _read_bounds(Path(args.controlled_bounds_csv)) if args.controlled_bounds_csv else {}

    rows: List[Dict[str, str]] = []
    for case_dir in case_dirs:
        if not case_dir.is_dir():
            continue
        case_name = case_dir.name
        run_log = case_dir / "run.log"
        exit_code_file = case_dir / "exit_code.txt"
        if not run_log.exists():
            continue

        try:
            controlled, full = _load_case_rows(run_log)
        except Exception as exc:
            rows.append(
                {
                    "case_name": case_name,
                    "status": "FAIL_PARSE",
                    "exit_code": "",
                    "error": str(exc),
                }
            )
            continue

        exit_code = ""
        if exit_code_file.exists():
            exit_code = exit_code_file.read_text(encoding="utf-8", errors="replace").strip()

        status = bounds_by_case.get(case_name, {}).get("status", "UNKNOWN")

        row = {
            "case_name": case_name,
            "status": status,
            "exit_code": exit_code,
            "controlled_final_psi0_span": _fmt(_to_float(controlled.get("final_psi0_span"))),
            "controlled_final_psi_proj_span": _fmt(_to_float(controlled.get("final_psi_proj_span"))),
            "controlled_final_psi_w0_span": _fmt(_to_float(controlled.get("final_psi_w0_span"))),
            "controlled_final_int_src_mode0_total_abs": _fmt(
                _to_float(controlled.get("final_int_src_mode0_total_abs"))
            ),
            "controlled_final_int_src_em_current_abs": _fmt(
                _to_float(controlled.get("final_int_src_em_current_abs"))
            ),
            "controlled_final_int_src_em_timelike_abs": _fmt(
                _to_float(controlled.get("final_int_src_em_timelike_abs"))
            ),
            "controlled_final_int_srcmomw_mode0": _fmt(_to_float(controlled.get("final_int_srcmomw_mode0"))),
            "controlled_final_int_srcenergy_mode0": _fmt(_to_float(controlled.get("final_int_srcenergy_mode0"))),
            "full_final_psi0_span": _fmt(_to_float(full.get("final_psi0_span"))),
            "full_final_jw_ew": _fmt(_to_float(full.get("final_jw_ew"))),
            "full_final_s_leak_abs": _fmt(_to_float(full.get("final_s_leak_abs"))),
            "error": "",
        }
        rows.append(row)

    if not rows:
        raise SystemExit(f"No case directories with run.log matched: {args.glob}")

    # Choose baseline row.
    baseline = None
    for row in rows:
        if row["case_name"] == args.baseline_case:
            baseline = row
            break
    if baseline is None:
        baseline = rows[0]

    b_psi0 = _to_float(baseline.get("controlled_final_psi0_span"))
    b_psip = _to_float(baseline.get("controlled_final_psi_proj_span"))
    b_psiw = _to_float(baseline.get("controlled_final_psi_w0_span"))
    b_src_total = _to_float(baseline.get("controlled_final_int_src_mode0_total_abs"))

    for row in rows:
        c_psi0 = _to_float(row.get("controlled_final_psi0_span"))
        c_psip = _to_float(row.get("controlled_final_psi_proj_span"))
        c_psiw = _to_float(row.get("controlled_final_psi_w0_span"))
        c_src_total = _to_float(row.get("controlled_final_int_src_mode0_total_abs"))

        row["psi0_reduction_vs_baseline"] = _fmt(_safe_ratio(b_psi0, c_psi0))
        row["psiproj_reduction_vs_baseline"] = _fmt(_safe_ratio(b_psip, c_psip))
        row["psiw0_reduction_vs_baseline"] = _fmt(_safe_ratio(b_psiw, c_psiw))
        row["src_mode0_total_reduction_vs_baseline"] = _fmt(_safe_ratio(b_src_total, c_src_total))
        row["psi0_ratio_full_over_controlled"] = _fmt(
            _safe_ratio(_to_float(row.get("full_final_psi0_span")), c_psi0)
        )

    def _psi0_sort_key(r: Dict[str, str]) -> float:
        v = _to_float(r.get("controlled_final_psi0_span"))
        if math.isnan(v):
            return float("inf")
        return abs(v)

    rows.sort(key=_psi0_sort_key)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_name",
        "status",
        "exit_code",
        "controlled_final_psi0_span",
        "controlled_final_psi_proj_span",
        "controlled_final_psi_w0_span",
        "full_final_psi0_span",
        "psi0_ratio_full_over_controlled",
        "psi0_reduction_vs_baseline",
        "psiproj_reduction_vs_baseline",
        "psiw0_reduction_vs_baseline",
        "controlled_final_int_src_mode0_total_abs",
        "src_mode0_total_reduction_vs_baseline",
        "controlled_final_int_src_em_current_abs",
        "controlled_final_int_src_em_timelike_abs",
        "controlled_final_int_srcmomw_mode0",
        "controlled_final_int_srcenergy_mode0",
        "full_final_jw_ew",
        "full_final_s_leak_abs",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    best = rows[0]
    best_case = best["case_name"]
    best_psi0 = best["controlled_final_psi0_span"]

    out_md = Path(args.out_md)
    lines: List[str] = []
    lines.append("# Step-126 Controlled Source Pinning Audit")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.")
    lines.append(f"Baseline case for reductions: `{baseline['case_name']}`.")
    lines.append("")
    lines.append(f"Best controlled-|psi0| case: `{best_case}` (`{best_psi0}`).")
    lines.append("")
    lines.append(
        "| case | status | ctrl psi0 | psi0 reduction vs base | ctrl psi_proj | ctrl psi_w0 | src_mode0_total | src reduction vs base | src_em_current | src_em_timelike | src_momw | src_energy | full/ctrl psi0 | full JwEw | full S_leak_abs |"
    )
    lines.append("|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            "| {case_name} | {status} | {controlled_final_psi0_span} | {psi0_reduction_vs_baseline} | {controlled_final_psi_proj_span} | {controlled_final_psi_w0_span} | {controlled_final_int_src_mode0_total_abs} | {src_mode0_total_reduction_vs_baseline} | {controlled_final_int_src_em_current_abs} | {controlled_final_int_src_em_timelike_abs} | {controlled_final_int_srcmomw_mode0} | {controlled_final_int_srcenergy_mode0} | {psi0_ratio_full_over_controlled} | {full_final_jw_ew} | {full_final_s_leak_abs} |".format(
                **row
            )
        )
    lines.append("")
    out_md.write_text("\n".join(lines), encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    print(f"baseline_case,{baseline['case_name']}")
    print(f"best_case,{best_case}")
    print(f"best_controlled_psi0,{best_psi0}")


if __name__ == "__main__":
    main()
