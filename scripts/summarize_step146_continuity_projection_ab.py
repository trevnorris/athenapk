#!/usr/bin/env python3
"""Summarize Step-146 continuity-consistent current projection A/B matrix."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import shlex
from pathlib import Path


def _to_float(value: str | None) -> float:
    if value is None:
        return math.nan
    s = value.strip()
    if s == "":
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def _fmt(x: float) -> str:
    if not math.isfinite(x):
        return "nan"
    return f"{x:.6e}"


def _safe_ratio(a: float, b: float) -> float:
    if not math.isfinite(a) or not math.isfinite(b) or b == 0.0:
        return math.nan
    return a / b


def _parse_case_rows(run_log: Path):
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
    cvals = next(csv.reader([controlled]))
    fvals = next(csv.reader([full]))
    return dict(zip(keys, cvals)), dict(zip(keys, fvals))


def _parse_projection_flags(run_log: Path) -> tuple[str, str]:
    text = run_log.read_text(encoding="utf-8", errors="replace")
    cmd_line = ""
    for line in text.splitlines():
        if line.startswith("command,"):
            cmd_line = line.split(",", 1)[1]
            break
    if not cmd_line:
        return "unset(default=false)", "unset(default=false)"

    tokens = shlex.split(cmd_line)
    controlled = "unset(default=false)"
    full = "unset(default=false)"
    key = "modes4d/continuity_consistent_current_projection_enable="

    for i, tok in enumerate(tokens[:-1]):
        val = tokens[i + 1]
        if tok == "--controlled-arg" and val.startswith(key):
            controlled = val.split("=", 1)[1]
        elif tok == "--full-arg" and val.startswith(key):
            full = val.split("=", 1)[1]
    return controlled, full


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Case glob (e.g. /path/out/*)")
    p.add_argument("--baseline-case", default="baseline_projection_off")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-md", required=True)
    args = p.parse_args()

    case_dirs = [Path(x) for x in sorted(glob.glob(args.glob)) if Path(x).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        if not run_log.exists():
            continue
        case_name = case_dir.name
        exit_code = ""
        exit_path = case_dir / "exit_code.txt"
        if exit_path.exists():
            exit_code = exit_path.read_text(encoding="utf-8", errors="replace").strip()

        controlled_flag, full_flag = _parse_projection_flags(run_log)
        try:
            controlled, full = _parse_case_rows(run_log)
        except Exception as exc:
            rows.append(
                {
                    "case_name": case_name,
                    "status": "FAIL_PARSE",
                    "exit_code": exit_code,
                    "projection_controlled": controlled_flag,
                    "projection_full": full_flag,
                    "error": str(exc),
                }
            )
            continue

        c_psi0 = _to_float(controlled.get("final_psi0_span"))
        c_div0 = _to_float(controlled.get("final_int_div_mode0_total_abs"))
        c_local_rate = _to_float(controlled.get("closure_local_mode0_final_max_abs_rate"))
        c_ledger_rate = _to_float(controlled.get("em_bulk_ledger_max_abs_rate"))

        f_psi0 = _to_float(full.get("final_psi0_span"))
        f_div0 = _to_float(full.get("final_int_div_mode0_total_abs"))
        f_local_rate = _to_float(full.get("closure_local_mode0_final_max_abs_rate"))
        f_ledger_rate = _to_float(full.get("em_bulk_ledger_max_abs_rate"))
        f_jwew = abs(_to_float(full.get("final_jw_ew")))
        f_sleak = _to_float(full.get("final_s_leak_abs"))

        rows.append(
            {
                "case_name": case_name,
                "status": "PASS" if exit_code == "0" else "FAIL_RUN",
                "exit_code": exit_code,
                "projection_controlled": controlled_flag,
                "projection_full": full_flag,
                "controlled_final_psi0_span": _fmt(c_psi0),
                "controlled_final_int_div_mode0_total_abs": _fmt(c_div0),
                "controlled_closure_local_mode0_final_max_abs_rate": _fmt(c_local_rate),
                "controlled_em_bulk_ledger_max_abs_rate": _fmt(c_ledger_rate),
                "full_final_psi0_span": _fmt(f_psi0),
                "full_final_int_div_mode0_total_abs": _fmt(f_div0),
                "full_closure_local_mode0_final_max_abs_rate": _fmt(f_local_rate),
                "full_em_bulk_ledger_max_abs_rate": _fmt(f_ledger_rate),
                "full_final_jw_ew_abs": _fmt(f_jwew),
                "full_final_s_leak_abs": _fmt(f_sleak),
                "full_closure_local_mode0_status": full.get("closure_local_mode0_status", ""),
                "full_transport_closure_status": full.get("transport_closure_status", ""),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", ""),
                "error": "",
            }
        )

    baseline = next((r for r in rows if r.get("case_name") == args.baseline_case), None)
    if baseline is None and rows:
        baseline = rows[0]

    b_c_div0 = _to_float(baseline.get("controlled_final_int_div_mode0_total_abs")) if baseline else math.nan
    b_c_ledger = _to_float(baseline.get("controlled_em_bulk_ledger_max_abs_rate")) if baseline else math.nan
    b_f_div0 = _to_float(baseline.get("full_final_int_div_mode0_total_abs")) if baseline else math.nan
    b_f_local = _to_float(baseline.get("full_closure_local_mode0_final_max_abs_rate")) if baseline else math.nan
    b_f_ledger = _to_float(baseline.get("full_em_bulk_ledger_max_abs_rate")) if baseline else math.nan

    for r in rows:
        r["controlled_div_mode0_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_c_div0, _to_float(r.get("controlled_final_int_div_mode0_total_abs")))
        )
        r["controlled_ledger_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_c_ledger, _to_float(r.get("controlled_em_bulk_ledger_max_abs_rate")))
        )
        r["full_div_mode0_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_f_div0, _to_float(r.get("full_final_int_div_mode0_total_abs")))
        )
        r["full_local_rate_reduction_vs_baseline"] = _fmt(
            _safe_ratio(
                b_f_local, _to_float(r.get("full_closure_local_mode0_final_max_abs_rate"))
            )
        )
        r["full_ledger_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_f_ledger, _to_float(r.get("full_em_bulk_ledger_max_abs_rate")))
        )

    rows.sort(key=lambda r: r.get("case_name", ""))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_name",
        "status",
        "exit_code",
        "projection_controlled",
        "projection_full",
        "controlled_final_psi0_span",
        "controlled_final_int_div_mode0_total_abs",
        "controlled_closure_local_mode0_final_max_abs_rate",
        "controlled_em_bulk_ledger_max_abs_rate",
        "controlled_div_mode0_reduction_vs_baseline",
        "controlled_ledger_reduction_vs_baseline",
        "full_final_psi0_span",
        "full_final_int_div_mode0_total_abs",
        "full_closure_local_mode0_final_max_abs_rate",
        "full_em_bulk_ledger_max_abs_rate",
        "full_div_mode0_reduction_vs_baseline",
        "full_local_rate_reduction_vs_baseline",
        "full_ledger_reduction_vs_baseline",
        "full_final_jw_ew_abs",
        "full_final_s_leak_abs",
        "full_closure_local_mode0_status",
        "full_transport_closure_status",
        "full_em_bulk_ledger_status",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.out_md)
    lines = []
    lines.append("# Step-146 Continuity-Consistent Current Projection A/B")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from `{args.glob}`.")
    if baseline:
        lines.append(f"Baseline case: `{baseline['case_name']}`.")
    lines.append("")
    lines.append(
        "| case | status | proj(ctrl/full) | ctrl div0 | ctrl ledger | full div0 | full local rate | full ledger | full JwEw | full Sleak | full local/transport/ledger |"
    )
    lines.append("|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|:---|")
    for r in rows:
        lines.append(
            "| {case_name} | {status} | {projection_controlled}/{projection_full} | "
            "{controlled_final_int_div_mode0_total_abs} | {controlled_em_bulk_ledger_max_abs_rate} | "
            "{full_final_int_div_mode0_total_abs} | {full_closure_local_mode0_final_max_abs_rate} | "
            "{full_em_bulk_ledger_max_abs_rate} | {full_final_jw_ew_abs} | {full_final_s_leak_abs} | "
            "{full_closure_local_mode0_status}/{full_transport_closure_status}/{full_em_bulk_ledger_status} |".format(
                **r
            )
        )
    lines.append("")
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    if baseline:
        print(f"baseline_case,{baseline['case_name']}")


if __name__ == "__main__":
    main()
