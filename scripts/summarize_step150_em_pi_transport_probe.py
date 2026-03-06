#!/usr/bin/env python3
"""Summarize Step-150 EM pi transport probe cases."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import re
from pathlib import Path


def to_float(value: str | None) -> float:
    if value is None:
        return math.nan
    text = value.strip()
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def safe_ratio(baseline: float, value: float) -> float:
    if not math.isfinite(baseline) or not math.isfinite(value) or value == 0.0:
        return math.nan
    return baseline / value


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


def last_cmdline_value(run_log: Path, key: str) -> str:
    pattern = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    text = run_log.read_text(encoding="utf-8", errors="replace")
    cmdline = ""
    for line in text.splitlines():
        if line.startswith("command,"):
            cmdline = line
            break
    if not cmdline:
        return ""
    matches = pattern.findall(cmdline)
    return matches[-1] if matches else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Case glob (absolute or relative)")
    parser.add_argument("--baseline-case", default="baseline_pi_transport")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows: list[dict[str, str]] = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        if not run_log.exists():
            continue

        exit_code = ""
        exit_path = case_dir / "exit_code.txt"
        if exit_path.exists():
            exit_code = exit_path.read_text(encoding="utf-8", errors="replace").strip()

        try:
            _, full = parse_case_rows(run_log)
        except Exception as exc:
            rows.append(
                {
                    "case_name": case_dir.name,
                    "status": "FAIL_PARSE",
                    "exit_code": exit_code,
                    "error": str(exc),
                }
            )
            continue

        rows.append(
            {
                "case_name": case_dir.name,
                "status": "PASS" if exit_code == "0" else "FAIL_RUN",
                "exit_code": exit_code,
                "full_closure_local_mode0_status": full.get("closure_local_mode0_status", ""),
                "full_transport_closure_status": full.get("transport_closure_status", ""),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", ""),
                "full_em_bulk_ledger_max_abs_rate": fmt(
                    to_float(full.get("em_bulk_ledger_max_abs_rate"))
                ),
                "full_final_jw_ew_abs": fmt(abs(to_float(full.get("final_jw_ew")))),
                "full_final_s_leak_abs": fmt(to_float(full.get("final_s_leak_abs"))),
                "full_pi0_transport_max_norm": fmt(
                    to_float(full.get("pi0_transport_max_norm"))
                ),
                "full_pi0_transport_max_abs_rate": fmt(
                    to_float(full.get("pi0_transport_max_abs_rate"))
                ),
                "full_pi0_transport_status": full.get("pi0_transport_status", ""),
                "full_piy_transport_max_norm": fmt(
                    to_float(full.get("piy_transport_max_norm"))
                ),
                "full_piy_transport_max_abs_rate": fmt(
                    to_float(full.get("piy_transport_max_abs_rate"))
                ),
                "full_piy_transport_status": full.get("piy_transport_status", ""),
                "full_piw_transport_max_norm": fmt(
                    to_float(full.get("piw_transport_max_norm"))
                ),
                "full_piw_transport_max_abs_rate": fmt(
                    to_float(full.get("piw_transport_max_abs_rate"))
                ),
                "full_piw_transport_status": full.get("piw_transport_status", ""),
                "piw_mode1_amp": last_cmdline_value(run_log, "problem/harris_4d/piw_mode1_amp"),
                "piy_mode1_amp": last_cmdline_value(run_log, "problem/harris_4d/piy_mode1_amp"),
                "a0_mode2_amp": last_cmdline_value(run_log, "problem/harris_4d/a0_mode2_amp"),
                "piw_mode2_amp": last_cmdline_value(run_log, "problem/harris_4d/piw_mode2_amp"),
                "error": "",
            }
        )

    baseline = next((row for row in rows if row.get("case_name") == args.baseline_case), None)
    if baseline is None and rows:
        baseline = rows[0]

    baseline_pi0 = to_float(baseline.get("full_pi0_transport_max_abs_rate")) if baseline else math.nan
    baseline_piy = to_float(baseline.get("full_piy_transport_max_abs_rate")) if baseline else math.nan
    baseline_piw = to_float(baseline.get("full_piw_transport_max_abs_rate")) if baseline else math.nan
    baseline_ledger = (
        to_float(baseline.get("full_em_bulk_ledger_max_abs_rate")) if baseline else math.nan
    )

    for row in rows:
        row["pi0_reduction_vs_baseline"] = fmt(
            safe_ratio(baseline_pi0, to_float(row.get("full_pi0_transport_max_abs_rate")))
        )
        row["piy_reduction_vs_baseline"] = fmt(
            safe_ratio(baseline_piy, to_float(row.get("full_piy_transport_max_abs_rate")))
        )
        row["piw_reduction_vs_baseline"] = fmt(
            safe_ratio(baseline_piw, to_float(row.get("full_piw_transport_max_abs_rate")))
        )
        row["ledger_reduction_vs_baseline"] = fmt(
            safe_ratio(
                baseline_ledger, to_float(row.get("full_em_bulk_ledger_max_abs_rate"))
            )
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_name",
        "status",
        "exit_code",
        "full_closure_local_mode0_status",
        "full_transport_closure_status",
        "full_em_bulk_ledger_status",
        "full_em_bulk_ledger_max_abs_rate",
        "full_final_jw_ew_abs",
        "full_final_s_leak_abs",
        "full_pi0_transport_max_norm",
        "full_pi0_transport_max_abs_rate",
        "full_pi0_transport_status",
        "full_piy_transport_max_norm",
        "full_piy_transport_max_abs_rate",
        "full_piy_transport_status",
        "full_piw_transport_max_norm",
        "full_piw_transport_max_abs_rate",
        "full_piw_transport_status",
        "pi0_reduction_vs_baseline",
        "piy_reduction_vs_baseline",
        "piw_reduction_vs_baseline",
        "ledger_reduction_vs_baseline",
        "piw_mode1_amp",
        "piy_mode1_amp",
        "a0_mode2_amp",
        "piw_mode2_amp",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    lines = [
        "# Step-150 EM Pi Transport Probe Summary",
        "",
        f"Parsed `{len(rows)}` cases from `{args.glob}`.",
    ]
    if baseline:
        lines.append(f"Baseline case: `{baseline['case_name']}`.")
    lines.extend(
        [
            "",
            "| case | status | local | transport | ledger | pi0 abs rate | piy abs rate | piw abs rate | pi0/piy/piw reduction vs baseline | ledger reduction | seeds (piw1, piy1, a0m2, piw2) |",
            "|:---|:---|:---|:---|:---|---:|---:|---:|:---|---:|:---|",
        ]
    )
    for row in rows:
        lines.append(
            "| {case_name} | {status} | {full_closure_local_mode0_status} | "
            "{full_transport_closure_status} | {full_em_bulk_ledger_status} | "
            "{full_pi0_transport_max_abs_rate} | {full_piy_transport_max_abs_rate} | "
            "{full_piw_transport_max_abs_rate} | {pi0_reduction_vs_baseline} / "
            "{piy_reduction_vs_baseline} / {piw_reduction_vs_baseline} | "
            "{ledger_reduction_vs_baseline} | {piw_mode1_amp}, {piy_mode1_amp}, "
            "{a0_mode2_amp}, {piw_mode2_amp} |".format(**row)
        )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    if baseline:
        print(f"baseline_case,{baseline['case_name']}")


if __name__ == "__main__":
    main()
