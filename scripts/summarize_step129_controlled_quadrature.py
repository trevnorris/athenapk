#!/usr/bin/env python3
"""Summarize Step-129 controlled quadrature-collapse probe outputs."""

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
    lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()

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


def _last_cmdline_value(run_log: Path, key: str) -> str:
    pattern = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    text = run_log.read_text(encoding="utf-8", errors="replace")
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
    parser.add_argument("--glob", required=True, help="Case directory glob")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--out-md", required=True, help="Output Markdown path")
    parser.add_argument(
        "--baseline-case",
        default="ctrl_baseline",
        help="Case name used as baseline for reductions",
    )
    args = parser.parse_args()

    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
    rows: List[Dict[str, str]] = []

    for case_dir in case_dirs:
        if not case_dir.is_dir():
            continue
        case_name = case_dir.name
        run_log = case_dir / "run.log"
        if not run_log.exists():
            continue
        exit_code_file = case_dir / "exit_code.txt"
        exit_code = exit_code_file.read_text(encoding="utf-8", errors="replace").strip() if exit_code_file.exists() else ""
        try:
            controlled, full = _load_case_rows(run_log)
        except Exception as exc:
            rows.append(
                {
                    "case_name": case_name,
                    "exit_code": exit_code,
                    "status": "FAIL_PARSE",
                    "error": str(exc),
                }
            )
            continue

        c_psi0 = _to_float(controlled.get("final_psi0_span"))
        c_psip = _to_float(controlled.get("final_psi_proj_span"))
        c_psiw = _to_float(controlled.get("final_psi_w0_span"))
        f_psi0 = _to_float(full.get("final_psi0_span"))
        row = {
            "case_name": case_name,
            "exit_code": exit_code,
            "status": "PASS" if exit_code == "0" else "FAIL",
            "controlled_final_psi0_span": _fmt(c_psi0),
            "controlled_final_psi_proj_span": _fmt(c_psip),
            "controlled_final_psi_w0_span": _fmt(c_psiw),
            "full_final_psi0_span": _fmt(f_psi0),
            "psi0_ratio_full_over_controlled": _fmt(_safe_ratio(f_psi0, c_psi0)),
            "controlled_closure_local_mode0_final_max_abs_rate": _fmt(
                _to_float(controlled.get("closure_local_mode0_final_max_abs_rate"))
            ),
            "controlled_em_bulk_ledger_max_abs_rate": _fmt(
                _to_float(controlled.get("em_bulk_ledger_max_abs_rate"))
            ),
            "full_final_jw_ew": _fmt(_to_float(full.get("final_jw_ew"))),
            "full_final_s_leak_abs": _fmt(_to_float(full.get("final_s_leak_abs"))),
            "controlled_hard_limit_enable": _last_cmdline_value(
                run_log, "modes4d/hard_controlled_limit_enable"
            ),
            "controlled_n_quadrature": _last_cmdline_value(run_log, "modes4d/n_quadrature"),
            "diffusion_resistivity": _last_cmdline_value(run_log, "diffusion/resistivity"),
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
    b_ledger = _to_float(baseline.get("controlled_em_bulk_ledger_max_abs_rate"))

    for row in rows:
        row["psi0_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_psi0, _to_float(row.get("controlled_final_psi0_span")))
        )
        row["ledger_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_ledger, _to_float(row.get("controlled_em_bulk_ledger_max_abs_rate")))
        )

    def _sort_key(r: Dict[str, str]) -> float:
        v = _to_float(r.get("controlled_final_psi0_span"))
        if math.isnan(v):
            return float("inf")
        return abs(v)

    rows.sort(key=_sort_key)
    best = rows[0]

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_name",
        "exit_code",
        "status",
        "controlled_final_psi0_span",
        "controlled_final_psi_proj_span",
        "controlled_final_psi_w0_span",
        "full_final_psi0_span",
        "psi0_ratio_full_over_controlled",
        "psi0_reduction_vs_baseline",
        "controlled_closure_local_mode0_final_max_abs_rate",
        "controlled_em_bulk_ledger_max_abs_rate",
        "ledger_reduction_vs_baseline",
        "full_final_jw_ew",
        "full_final_s_leak_abs",
        "controlled_hard_limit_enable",
        "controlled_n_quadrature",
        "diffusion_resistivity",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.out_md)
    lines: List[str] = []
    lines.append("# Step-129 Controlled Quadrature-Collapse Probe")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from `{args.glob}`.")
    lines.append(f"Baseline case: `{baseline['case_name']}`.")
    lines.append(
        "Best (lowest controlled psi0) case: "
        f"`{best['case_name']}` with controlled psi0=`{best.get('controlled_final_psi0_span', 'nan')}`."
    )
    lines.append("")
    lines.append("| case | status | controlled psi0 | psi0 reduction vs baseline | ledger reduction vs baseline | hard-limit | n_quadrature (controlled) | resistivity |")
    lines.append("|---|---|---:|---:|---:|---|---:|---|")
    for row in rows:
        lines.append(
            "| "
            f"{row['case_name']} | {row['status']} | {row['controlled_final_psi0_span']} | "
            f"{row['psi0_reduction_vs_baseline']} | {row['ledger_reduction_vs_baseline']} | "
            f"{row['controlled_hard_limit_enable'] or 'false'} | "
            f"{row['controlled_n_quadrature'] or 'default'} | {row['diffusion_resistivity'] or 'default'} |"
        )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    print(f"best_case,{best['case_name']}")
    print(f"best_controlled_psi0,{best.get('controlled_final_psi0_span', 'nan')}")


if __name__ == "__main__":
    main()
