#!/usr/bin/env python3
"""Summarize controlled non-source pinning audit from run logs."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import re
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


def _last_cmdline_value(run_log: Path, key: str) -> str:
    pattern = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    try:
        text = run_log.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    # Prefer the explicit "command,..." line emitted by run scripts.
    cmd_line = ""
    for line in text.splitlines():
        if line.startswith("command,"):
            cmd_line = line
            break
    if not cmd_line:
        return ""
    matches = pattern.findall(cmd_line)
    if not matches:
        return ""
    return matches[-1]


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
        help="Case name used as baseline for reductions",
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

        c_psi0 = _to_float(controlled.get("final_psi0_span"))
        c_psip = _to_float(controlled.get("final_psi_proj_span"))
        c_psiw = _to_float(controlled.get("final_psi_w0_span"))
        f_psi0 = _to_float(full.get("final_psi0_span"))

        row = {
            "case_name": case_name,
            "status": status,
            "exit_code": exit_code,
            "controlled_final_psi0_span": _fmt(c_psi0),
            "controlled_final_psi_proj_span": _fmt(c_psip),
            "controlled_final_psi_w0_span": _fmt(c_psiw),
            "full_final_psi0_span": _fmt(f_psi0),
            "psi0_ratio_full_over_controlled": _fmt(_safe_ratio(f_psi0, c_psi0)),
            "controlled_final_int_div_mode0_total_abs": _fmt(
                _to_float(controlled.get("final_int_div_mode0_total_abs"))
            ),
            "controlled_final_int_div_mode0_plasma_abs": _fmt(
                _to_float(controlled.get("final_int_div_mode0_plasma_abs"))
            ),
            "controlled_final_int_div_mode0_em_abs": _fmt(
                _to_float(controlled.get("final_int_div_mode0_em_abs"))
            ),
            "controlled_closure_local_mode0_final_max_abs_rate": _fmt(
                _to_float(controlled.get("closure_local_mode0_final_max_abs_rate"))
            ),
            "controlled_em_bulk_ledger_max_abs_rate": _fmt(
                _to_float(controlled.get("em_bulk_ledger_max_abs_rate"))
            ),
            "full_final_jw_ew": _fmt(_to_float(full.get("final_jw_ew"))),
            "full_final_s_leak_abs": _fmt(_to_float(full.get("final_s_leak_abs"))),
            "diffusion_resistivity": _last_cmdline_value(run_log, "diffusion/resistivity"),
            "plasma_rho_divj_gain": _last_cmdline_value(run_log, "modes4d/plasma_rho_divj_gain"),
            "plasma_pressure_rusanov_gain": _last_cmdline_value(
                run_log, "modes4d/plasma_pressure_rusanov_gain"
            ),
            "plasma_pressure_signal_speed_cap": _last_cmdline_value(
                run_log, "modes4d/plasma_pressure_signal_speed_cap"
            ),
            "plasma_pressure_flux_relative_cap": _last_cmdline_value(
                run_log, "modes4d/plasma_pressure_flux_relative_cap"
            ),
            "error": "",
        }
        rows.append(row)

    if not rows:
        raise SystemExit(f"No case directories with run.log matched: {args.glob}")

    baseline = None
    for row in rows:
        if row["case_name"] == args.baseline_case:
            baseline = row
            break
    if baseline is None:
        baseline = rows[0]

    b_psi0 = _to_float(baseline.get("controlled_final_psi0_span"))
    b_div_total = _to_float(baseline.get("controlled_final_int_div_mode0_total_abs"))
    b_ledger = _to_float(baseline.get("controlled_em_bulk_ledger_max_abs_rate"))

    for row in rows:
        c_psi0 = _to_float(row.get("controlled_final_psi0_span"))
        c_div_total = _to_float(row.get("controlled_final_int_div_mode0_total_abs"))
        c_ledger = _to_float(row.get("controlled_em_bulk_ledger_max_abs_rate"))
        row["psi0_reduction_vs_baseline"] = _fmt(_safe_ratio(b_psi0, c_psi0))
        row["div_mode0_total_reduction_vs_baseline"] = _fmt(_safe_ratio(b_div_total, c_div_total))
        row["ledger_reduction_vs_baseline"] = _fmt(_safe_ratio(b_ledger, c_ledger))

    def _psi0_sort_key(r: Dict[str, str]) -> float:
        v = _to_float(r.get("controlled_final_psi0_span"))
        if math.isnan(v):
            return float("inf")
        return abs(v)

    rows.sort(key=_psi0_sort_key)
    best = rows[0]

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
        "controlled_final_int_div_mode0_total_abs",
        "controlled_final_int_div_mode0_plasma_abs",
        "controlled_final_int_div_mode0_em_abs",
        "div_mode0_total_reduction_vs_baseline",
        "controlled_closure_local_mode0_final_max_abs_rate",
        "controlled_em_bulk_ledger_max_abs_rate",
        "ledger_reduction_vs_baseline",
        "full_final_jw_ew",
        "full_final_s_leak_abs",
        "diffusion_resistivity",
        "plasma_rho_divj_gain",
        "plasma_pressure_rusanov_gain",
        "plasma_pressure_signal_speed_cap",
        "plasma_pressure_flux_relative_cap",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.out_md)
    lines: List[str] = []
    lines.append("# Step-127 Non-Source Controlled Pinning Audit")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.")
    lines.append(f"Baseline case for reductions: `{baseline['case_name']}`.")
    lines.append("")
    lines.append(
        f"Best controlled-|psi0| case: `{best['case_name']}` (`{best['controlled_final_psi0_span']}`)."
    )
    lines.append("")
    lines.append(
        "| case | status | ctrl psi0 | psi0 reduction | div_mode0_total | div reduction | ledger max rate | ledger reduction | rho_divj | rusanov | speed cap | flux cap | resistivity | full/ctrl psi0 | full JwEw | full S_leak_abs |"
    )
    lines.append("|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|")
    for row in rows:
        lines.append(
            "| {case_name} | {status} | {controlled_final_psi0_span} | {psi0_reduction_vs_baseline} | {controlled_final_int_div_mode0_total_abs} | {div_mode0_total_reduction_vs_baseline} | {controlled_em_bulk_ledger_max_abs_rate} | {ledger_reduction_vs_baseline} | {plasma_rho_divj_gain} | {plasma_pressure_rusanov_gain} | {plasma_pressure_signal_speed_cap} | {plasma_pressure_flux_relative_cap} | {diffusion_resistivity} | {psi0_ratio_full_over_controlled} | {full_final_jw_ew} | {full_final_s_leak_abs} |".format(
                **row
            )
        )
    lines.append("")
    out_md.write_text("\n".join(lines), encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    print(f"baseline_case,{baseline['case_name']}")
    print(f"best_case,{best['case_name']}")
    print(f"best_controlled_psi0,{best['controlled_final_psi0_span']}")


if __name__ == "__main__":
    main()
