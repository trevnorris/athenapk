#!/usr/bin/env python3
"""Summarize Step-145 controlled plasma continuity-coupling audit."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import re
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
        for j in range(i + 1, min(i + 16, len(lines))):
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


def _last_cmdline_value(run_log: Path, key: str) -> str:
    pat = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    text = run_log.read_text(encoding="utf-8", errors="replace")
    cmd = ""
    for line in text.splitlines():
        if line.startswith("command,"):
            cmd = line
            break
    if not cmd:
        return ""
    vals = pat.findall(cmd)
    return vals[-1] if vals else ""


def _read_bounds(path: Path):
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {r.get("case_name", ""): r for r in rows}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Case glob (e.g. /path/out/*)")
    p.add_argument("--baseline-case", default="c3_baseline_strict")
    p.add_argument("--controlled-bounds-csv", default="")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-md", required=True)
    args = p.parse_args()

    case_dirs = [Path(x) for x in sorted(glob.glob(args.glob)) if Path(x).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    bounds = _read_bounds(Path(args.controlled_bounds_csv)) if args.controlled_bounds_csv else {}

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

        try:
            controlled, full = _parse_case_rows(run_log)
        except Exception as exc:
            rows.append(
                {
                    "case_name": case_name,
                    "status": "FAIL_PARSE",
                    "exit_code": exit_code,
                    "error": str(exc),
                }
            )
            continue

        c_psi0 = _to_float(controlled.get("final_psi0_span"))
        c_psip = _to_float(controlled.get("final_psi_proj_span"))
        c_psiw = _to_float(controlled.get("final_psi_w0_span"))
        c_div = _to_float(controlled.get("final_int_div_mode0_total_abs"))
        c_cont = _to_float(controlled.get("closure_local_mode0_final_max_abs_rate"))
        c_ledger = _to_float(controlled.get("em_bulk_ledger_max_abs_rate"))

        f_psi0 = _to_float(full.get("final_psi0_span"))
        f_jwew = _to_float(full.get("final_jw_ew"))
        f_sleak = _to_float(full.get("final_s_leak_abs"))
        f_cont_status = full.get("closure_local_mode0_status", "")
        f_transport_status = full.get("transport_closure_status", "")
        f_ledger_status = full.get("em_bulk_ledger_status", "")

        rows.append(
            {
                "case_name": case_name,
                "status": bounds.get(case_name, {}).get("status", "UNKNOWN"),
                "exit_code": exit_code,
                "controlled_final_psi0_span": _fmt(c_psi0),
                "controlled_final_psi_proj_span": _fmt(c_psip),
                "controlled_final_psi_w0_span": _fmt(c_psiw),
                "controlled_final_int_div_mode0_total_abs": _fmt(c_div),
                "controlled_closure_local_mode0_final_max_abs_rate": _fmt(c_cont),
                "controlled_em_bulk_ledger_max_abs_rate": _fmt(c_ledger),
                "full_final_psi0_span": _fmt(f_psi0),
                "psi0_ratio_full_over_controlled": _fmt(_safe_ratio(f_psi0, c_psi0)),
                "full_final_jw_ew": _fmt(f_jwew),
                "full_final_s_leak_abs": _fmt(f_sleak),
                "full_closure_local_mode0_status": f_cont_status,
                "full_transport_closure_status": f_transport_status,
                "full_em_bulk_ledger_status": f_ledger_status,
                "modes4d_em_source_current_gain": _last_cmdline_value(
                    run_log, "modes4d/em_source_current_gain"
                ),
                "modes4d_plasma_rho_divj_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_rho_divj_gain"
                ),
                "modes4d_plasma_pressure_rusanov_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_pressure_rusanov_gain"
                ),
                "modes4d_plasma_pressure_signal_speed_cap": _last_cmdline_value(
                    run_log, "modes4d/plasma_pressure_signal_speed_cap"
                ),
                "modes4d_plasma_pressure_flux_relative_cap": _last_cmdline_value(
                    run_log, "modes4d/plasma_pressure_flux_relative_cap"
                ),
                "hard_controlled_limit_strict_solver_path": _last_cmdline_value(
                    run_log, "modes4d/hard_controlled_limit_strict_solver_path"
                ),
                "error": "",
            }
        )

    # Baseline-reduction columns.
    baseline = next((r for r in rows if r.get("case_name") == args.baseline_case), None)
    if baseline is None and rows:
        baseline = rows[0]
    b_psi0 = _to_float(baseline.get("controlled_final_psi0_span")) if baseline else math.nan
    b_div = _to_float(baseline.get("controlled_final_int_div_mode0_total_abs")) if baseline else math.nan
    b_cont = _to_float(baseline.get("controlled_closure_local_mode0_final_max_abs_rate")) if baseline else math.nan
    b_ledger = _to_float(baseline.get("controlled_em_bulk_ledger_max_abs_rate")) if baseline else math.nan
    for r in rows:
        c_psi0 = _to_float(r.get("controlled_final_psi0_span"))
        c_div = _to_float(r.get("controlled_final_int_div_mode0_total_abs"))
        c_cont = _to_float(r.get("controlled_closure_local_mode0_final_max_abs_rate"))
        c_ledger = _to_float(r.get("controlled_em_bulk_ledger_max_abs_rate"))
        r["psi0_reduction_vs_baseline"] = _fmt(_safe_ratio(b_psi0, c_psi0))
        r["div_mode0_total_reduction_vs_baseline"] = _fmt(_safe_ratio(b_div, c_div))
        r["local_mode0_rate_reduction_vs_baseline"] = _fmt(_safe_ratio(b_cont, c_cont))
        r["ledger_reduction_vs_baseline"] = _fmt(_safe_ratio(b_ledger, c_ledger))

    def _sort_key(r):
        x = _to_float(r.get("controlled_em_bulk_ledger_max_abs_rate"))
        return abs(x) if math.isfinite(x) else float("inf")

    rows.sort(key=_sort_key)
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
        "controlled_final_int_div_mode0_total_abs",
        "controlled_closure_local_mode0_final_max_abs_rate",
        "controlled_em_bulk_ledger_max_abs_rate",
        "psi0_reduction_vs_baseline",
        "div_mode0_total_reduction_vs_baseline",
        "local_mode0_rate_reduction_vs_baseline",
        "ledger_reduction_vs_baseline",
        "full_final_psi0_span",
        "psi0_ratio_full_over_controlled",
        "full_final_jw_ew",
        "full_final_s_leak_abs",
        "full_closure_local_mode0_status",
        "full_transport_closure_status",
        "full_em_bulk_ledger_status",
        "modes4d_em_source_current_gain",
        "modes4d_plasma_rho_divj_gain",
        "modes4d_plasma_pressure_rusanov_gain",
        "modes4d_plasma_pressure_signal_speed_cap",
        "modes4d_plasma_pressure_flux_relative_cap",
        "hard_controlled_limit_strict_solver_path",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.out_md)
    lines = []
    lines.append("# Step-145 Controlled Plasma Continuity-Coupling Audit")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from `{args.glob}`.")
    if baseline:
        lines.append(f"Baseline case: `{baseline['case_name']}`.")
    lines.append(
        f"Best controlled ledger case: `{best['case_name']}` "
        f"(`{best['controlled_em_bulk_ledger_max_abs_rate']}`)."
    )
    lines.append("")
    lines.append(
        "| case | status | ctrl psi0 | psi0 red | ctrl div_mode0 | div red | ctrl local rate | local red | ctrl ledger rate | ledger red | current | rho_divj | rusanov | full JwEw | full S_leak_abs | full local/transport/ledger |"
    )
    lines.append("|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|")
    for r in rows:
        lines.append(
            "| {case_name} | {status} | {controlled_final_psi0_span} | {psi0_reduction_vs_baseline} | "
            "{controlled_final_int_div_mode0_total_abs} | {div_mode0_total_reduction_vs_baseline} | "
            "{controlled_closure_local_mode0_final_max_abs_rate} | {local_mode0_rate_reduction_vs_baseline} | "
            "{controlled_em_bulk_ledger_max_abs_rate} | {ledger_reduction_vs_baseline} | "
            "{modes4d_em_source_current_gain} | {modes4d_plasma_rho_divj_gain} | "
            "{modes4d_plasma_pressure_rusanov_gain} | {full_final_jw_ew} | {full_final_s_leak_abs} | "
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
    print(f"best_case,{best['case_name']}")
    print(f"best_controlled_ledger_rate,{best['controlled_em_bulk_ledger_max_abs_rate']}")


if __name__ == "__main__":
    main()
