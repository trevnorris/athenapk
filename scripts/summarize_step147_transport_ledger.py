#!/usr/bin/env python3
"""Summarize Step-147 transport/ledger isolation cases."""

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


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Case glob (e.g. /path/out/*)")
    p.add_argument("--baseline-case", default="baseline_proj_full")
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

        try:
            _, full = _parse_case_rows(run_log)
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

        transport_norm_sum = 0.0
        transport_keys = [
            "momx_transport_max_norm",
            "momy_transport_max_norm",
            "momz_transport_max_norm",
            "momw_transport_max_norm",
            "energy_transport_max_norm",
            "pi0_transport_max_norm",
            "pix_transport_max_norm",
            "piy_transport_max_norm",
            "piz_transport_max_norm",
            "piw_transport_max_norm",
        ]
        for k in transport_keys:
            v = _to_float(full.get(k))
            if math.isfinite(v):
                transport_norm_sum += abs(v)

        rows.append(
            {
                "case_name": case_name,
                "status": "PASS" if exit_code == "0" else "FAIL_RUN",
                "exit_code": exit_code,
                "full_closure_local_mode0_status": full.get("closure_local_mode0_status", ""),
                "full_transport_closure_status": full.get("transport_closure_status", ""),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", ""),
                "full_closure_local_mode0_final_max_abs_rate": _fmt(
                    _to_float(full.get("closure_local_mode0_final_max_abs_rate"))
                ),
                "full_transport_norm_sum": _fmt(transport_norm_sum),
                "full_em_bulk_ledger_max_abs_rate": _fmt(
                    _to_float(full.get("em_bulk_ledger_max_abs_rate"))
                ),
                "full_final_int_src_mode0_total_abs": _fmt(
                    _to_float(full.get("final_int_src_mode0_total_abs"))
                ),
                "full_final_int_div_mode0_total_abs": _fmt(
                    _to_float(full.get("final_int_div_mode0_total_abs"))
                ),
                "full_final_jw_ew_abs": _fmt(abs(_to_float(full.get("final_jw_ew")))),
                "full_final_s_leak_abs": _fmt(_to_float(full.get("final_s_leak_abs"))),
                "modes4d_cont_projection": _last_cmdline_value(
                    run_log, "modes4d/continuity_consistent_current_projection_enable"
                ),
                "modes4d_em_source_current_gain": _last_cmdline_value(
                    run_log, "modes4d/em_source_current_gain"
                ),
                "modes4d_em_source_timelike_gain": _last_cmdline_value(
                    run_log, "modes4d/em_source_timelike_gain"
                ),
                "modes4d_em_source_damping_gain": _last_cmdline_value(
                    run_log, "modes4d/em_source_damping_gain"
                ),
                "modes4d_plasma_momw_source_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_momw_source_gain"
                ),
                "modes4d_plasma_momw_pressure_source_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_momw_pressure_source_gain"
                ),
                "modes4d_plasma_rho_divj_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_rho_divj_gain"
                ),
                "modes4d_plasma_pressure_transport_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_pressure_transport_gain"
                ),
                "modes4d_plasma_pressure_rusanov_gain": _last_cmdline_value(
                    run_log, "modes4d/plasma_pressure_rusanov_gain"
                ),
                "error": "",
            }
        )

    baseline = next((r for r in rows if r.get("case_name") == args.baseline_case), None)
    if baseline is None and rows:
        baseline = rows[0]

    b_local = _to_float(baseline.get("full_closure_local_mode0_final_max_abs_rate")) if baseline else math.nan
    b_transport = _to_float(baseline.get("full_transport_norm_sum")) if baseline else math.nan
    b_ledger = _to_float(baseline.get("full_em_bulk_ledger_max_abs_rate")) if baseline else math.nan

    for r in rows:
        r["local_rate_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_local, _to_float(r.get("full_closure_local_mode0_final_max_abs_rate")))
        )
        r["transport_sum_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_transport, _to_float(r.get("full_transport_norm_sum")))
        )
        r["ledger_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_ledger, _to_float(r.get("full_em_bulk_ledger_max_abs_rate")))
        )

    rows.sort(key=lambda r: r.get("case_name", ""))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_name",
        "status",
        "exit_code",
        "full_closure_local_mode0_status",
        "full_transport_closure_status",
        "full_em_bulk_ledger_status",
        "full_closure_local_mode0_final_max_abs_rate",
        "full_transport_norm_sum",
        "full_em_bulk_ledger_max_abs_rate",
        "full_final_int_src_mode0_total_abs",
        "full_final_int_div_mode0_total_abs",
        "full_final_jw_ew_abs",
        "full_final_s_leak_abs",
        "local_rate_reduction_vs_baseline",
        "transport_sum_reduction_vs_baseline",
        "ledger_reduction_vs_baseline",
        "modes4d_cont_projection",
        "modes4d_em_source_current_gain",
        "modes4d_em_source_timelike_gain",
        "modes4d_em_source_damping_gain",
        "modes4d_plasma_momw_source_gain",
        "modes4d_plasma_momw_pressure_source_gain",
        "modes4d_plasma_rho_divj_gain",
        "modes4d_plasma_pressure_transport_gain",
        "modes4d_plasma_pressure_rusanov_gain",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.out_md)
    lines = []
    lines.append("# Step-147 Transport/Ledger Isolation Summary")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from `{args.glob}`.")
    if baseline:
        lines.append(f"Baseline case: `{baseline['case_name']}`.")
    lines.append("")
    lines.append(
        "| case | status | local | transport | ledger | local rate | transport sum | ledger max rate | JwEw | Sleak | local/transport/ledger reduction vs baseline |"
    )
    lines.append("|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|:---|")
    for r in rows:
        lines.append(
            "| {case_name} | {status} | {full_closure_local_mode0_status} | "
            "{full_transport_closure_status} | {full_em_bulk_ledger_status} | "
            "{full_closure_local_mode0_final_max_abs_rate} | {full_transport_norm_sum} | "
            "{full_em_bulk_ledger_max_abs_rate} | {full_final_jw_ew_abs} | "
            "{full_final_s_leak_abs} | {local_rate_reduction_vs_baseline} / "
            "{transport_sum_reduction_vs_baseline} / {ledger_reduction_vs_baseline} |".format(
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
