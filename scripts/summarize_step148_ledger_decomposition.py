#!/usr/bin/env python3
"""Summarize Step-148 ledger decomposition matrix."""

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


def _transport_norm_sum(row: dict[str, str]) -> float:
    keys = [
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
    total = 0.0
    for k in keys:
        v = _to_float(row.get(k))
        if math.isfinite(v):
            total += abs(v)
    return total


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Case glob (e.g. /path/out/*)")
    p.add_argument("--baseline-case", default="baseline_decomp")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-md", required=True)
    args = p.parse_args()

    case_dirs = [Path(x) for x in sorted(glob.glob(args.glob)) if Path(x).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows: list[dict[str, str]] = []
    for case_dir in case_dirs:
        case_name = case_dir.name
        run_log = case_dir / "run.log"
        exit_path = case_dir / "exit_code.txt"
        exit_code = exit_path.read_text(encoding="utf-8", errors="replace").strip() if exit_path.exists() else ""

        if not run_log.exists():
            rows.append(
                {
                    "case_name": case_name,
                    "status": "FAIL_PARSE",
                    "exit_code": exit_code,
                    "error": "missing run.log",
                }
            )
            continue

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

        row = {
            "case_name": case_name,
            "status": "PASS" if exit_code == "0" else "FAIL_RUN",
            "exit_code": exit_code,
            "full_local_status": full.get("closure_local_mode0_status", ""),
            "full_transport_status": full.get("transport_closure_status", ""),
            "full_ledger_status": full.get("em_bulk_ledger_status", ""),
            "full_local_rate": _fmt(_to_float(full.get("closure_local_mode0_final_max_abs_rate"))),
            "full_transport_norm_sum": _fmt(_transport_norm_sum(full)),
            "full_ledger_max_abs_rate": _fmt(_to_float(full.get("em_bulk_ledger_max_abs_rate"))),
            "full_bridge_best_label": full.get("em_bulk_ledger_bridge_best_label", ""),
            "full_bridge_best_max_abs_rate": _fmt(_to_float(full.get("em_bulk_ledger_bridge_best_max_abs_rate"))),
            "full_div_mode0_total_abs": _fmt(_to_float(full.get("final_int_div_mode0_total_abs"))),
            "full_div_mode0_plasma_abs": _fmt(_to_float(full.get("final_int_div_mode0_plasma_abs"))),
            "full_div_mode0_em_abs": _fmt(_to_float(full.get("final_int_div_mode0_em_abs"))),
            "full_src_mode0_total_abs": _fmt(_to_float(full.get("final_int_src_mode0_total_abs"))),
            "full_src_mode0_plasma_abs": _fmt(_to_float(full.get("final_int_src_mode0_plasma_abs"))),
            "full_src_mode0_em_abs": _fmt(_to_float(full.get("final_int_src_mode0_em_abs"))),
            "full_src_mode0_timelike_abs": _fmt(_to_float(full.get("final_int_src_mode0_timelike_abs"))),
            "full_src_em_laplacian_abs": _fmt(_to_float(full.get("final_int_src_em_laplacian_abs"))),
            "full_src_em_mass_abs": _fmt(_to_float(full.get("final_int_src_em_mass_abs"))),
            "full_src_em_current_abs": _fmt(_to_float(full.get("final_int_src_em_current_abs"))),
            "full_src_em_damping_abs": _fmt(_to_float(full.get("final_int_src_em_damping_abs"))),
            "full_src_em_timelike_abs": _fmt(_to_float(full.get("final_int_src_em_timelike_abs"))),
            "full_jwew_abs": _fmt(abs(_to_float(full.get("final_jw_ew")))),
            "full_sleak_abs": _fmt(_to_float(full.get("final_s_leak_abs"))),
            "arg_em_conservative_transport": _last_cmdline_value(run_log, "modes4d/em_conservative_transport"),
            "arg_cont_projection": _last_cmdline_value(
                run_log, "modes4d/continuity_consistent_current_projection_enable"
            ),
            "arg_em_source_current_gain": _last_cmdline_value(run_log, "modes4d/em_source_current_gain"),
            "arg_em_source_timelike_gain": _last_cmdline_value(run_log, "modes4d/em_source_timelike_gain"),
            "arg_em_source_damping_gain": _last_cmdline_value(run_log, "modes4d/em_source_damping_gain"),
            "arg_plasma_force_source_gain": _last_cmdline_value(run_log, "modes4d/plasma_force_source_gain"),
            "arg_plasma_momw_source_gain": _last_cmdline_value(run_log, "modes4d/plasma_momw_source_gain"),
            "arg_plasma_momw_pressure_source_gain": _last_cmdline_value(
                run_log, "modes4d/plasma_momw_pressure_source_gain"
            ),
            "arg_plasma_energy_source_gain": _last_cmdline_value(run_log, "modes4d/plasma_energy_source_gain"),
            "arg_plasma_w_flux_source_gain": _last_cmdline_value(run_log, "modes4d/plasma_w_flux_source_gain"),
            "arg_plasma_rho_divj_gain": _last_cmdline_value(run_log, "modes4d/plasma_rho_divj_gain"),
            "arg_plasma_pressure_transport_gain": _last_cmdline_value(
                run_log, "modes4d/plasma_pressure_transport_gain"
            ),
            "arg_plasma_pressure_rusanov_gain": _last_cmdline_value(
                run_log, "modes4d/plasma_pressure_rusanov_gain"
            ),
            "error": "",
        }
        rows.append(row)

    baseline = next((r for r in rows if r.get("case_name") == args.baseline_case), None)
    if baseline is None and rows:
        baseline = rows[0]

    b_local = _to_float(baseline.get("full_local_rate")) if baseline else math.nan
    b_transport = _to_float(baseline.get("full_transport_norm_sum")) if baseline else math.nan
    b_ledger = _to_float(baseline.get("full_ledger_max_abs_rate")) if baseline else math.nan
    for r in rows:
        r["local_reduction_vs_baseline"] = _fmt(_safe_ratio(b_local, _to_float(r.get("full_local_rate"))))
        r["transport_reduction_vs_baseline"] = _fmt(
            _safe_ratio(b_transport, _to_float(r.get("full_transport_norm_sum")))
        )
        r["ledger_reduction_vs_baseline"] = _fmt(_safe_ratio(b_ledger, _to_float(r.get("full_ledger_max_abs_rate"))))

    rows.sort(key=lambda x: x.get("case_name", ""))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_name",
        "status",
        "exit_code",
        "full_local_status",
        "full_transport_status",
        "full_ledger_status",
        "full_local_rate",
        "full_transport_norm_sum",
        "full_ledger_max_abs_rate",
        "local_reduction_vs_baseline",
        "transport_reduction_vs_baseline",
        "ledger_reduction_vs_baseline",
        "full_bridge_best_label",
        "full_bridge_best_max_abs_rate",
        "full_div_mode0_total_abs",
        "full_div_mode0_plasma_abs",
        "full_div_mode0_em_abs",
        "full_src_mode0_total_abs",
        "full_src_mode0_plasma_abs",
        "full_src_mode0_em_abs",
        "full_src_mode0_timelike_abs",
        "full_src_em_laplacian_abs",
        "full_src_em_mass_abs",
        "full_src_em_current_abs",
        "full_src_em_damping_abs",
        "full_src_em_timelike_abs",
        "full_jwew_abs",
        "full_sleak_abs",
        "arg_em_conservative_transport",
        "arg_cont_projection",
        "arg_em_source_current_gain",
        "arg_em_source_timelike_gain",
        "arg_em_source_damping_gain",
        "arg_plasma_force_source_gain",
        "arg_plasma_momw_source_gain",
        "arg_plasma_momw_pressure_source_gain",
        "arg_plasma_energy_source_gain",
        "arg_plasma_w_flux_source_gain",
        "arg_plasma_rho_divj_gain",
        "arg_plasma_pressure_transport_gain",
        "arg_plasma_pressure_rusanov_gain",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.out_md)
    lines: list[str] = []
    lines.append("# Step-148 Ledger Decomposition Summary")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from `{args.glob}`.")
    if baseline:
        lines.append(f"Baseline case: `{baseline['case_name']}`.")
    lines.append("")
    lines.append(
        "| case | status | local | transport | ledger | local rate | transport sum | ledger max rate | bridge best | JwEw | Sleak | local/transport/ledger reduction vs baseline |"
    )
    lines.append("|:---|:---|:---|:---|:---|---:|---:|---:|:---|---:|---:|:---|")
    for r in rows:
        lines.append(
            "| {case_name} | {status} | {full_local_status} | {full_transport_status} | "
            "{full_ledger_status} | {full_local_rate} | {full_transport_norm_sum} | "
            "{full_ledger_max_abs_rate} | {full_bridge_best_label} ({full_bridge_best_max_abs_rate}) | "
            "{full_jwew_abs} | {full_sleak_abs} | {local_reduction_vs_baseline} / "
            "{transport_reduction_vs_baseline} / {ledger_reduction_vs_baseline} |".format(**r)
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
