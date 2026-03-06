#!/usr/bin/env python3
"""Summarize controlled-branch boundedness from harris_scan_matrix run logs."""

from __future__ import annotations

import argparse
import csv
import glob
import math
from pathlib import Path
from typing import Dict, List, Optional


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


def _load_case_rows(run_log: Path) -> tuple[Dict[str, str], Dict[str, str]]:
    text = run_log.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    header = None
    controlled = None
    full = None
    for idx, line in enumerate(lines):
        if not line.startswith("case,final_psi0_span,"):
            continue
        header = line
        # The next lines are expected to contain controlled/full CSV rows.
        for j in range(idx + 1, min(idx + 8, len(lines))):
            if lines[j].startswith("controlled,"):
                controlled = lines[j]
            elif lines[j].startswith("full,"):
                full = lines[j]
        if header and controlled and full:
            break

    if header is None or controlled is None or full is None:
        raise RuntimeError(f"Could not find controlled/full CSV rows in {run_log}")

    keys = next(csv.reader([header]))
    controlled_vals = next(csv.reader([controlled]))
    full_vals = next(csv.reader([full]))
    return dict(zip(keys, controlled_vals)), dict(zip(keys, full_vals))


def _status_from_limits(v0: float, vp: float, vw: float, lim0: float, limp: float, limw: float) -> str:
    if any(math.isnan(v) or math.isinf(v) for v in (v0, vp, vw)):
        return "FAIL_NONFINITE"
    if abs(v0) > lim0 or abs(vp) > limp or abs(vw) > limw:
        return "FAIL_UNBOUNDED"
    return "PASS"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Case directory glob (e.g., /path/out/*)")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--out-md", required=True, help="Output Markdown path")
    parser.add_argument("--psi0-max", type=float, default=1.0e3, help="Max allowed |controlled final_psi0_span|")
    parser.add_argument("--psiproj-max", type=float, default=1.0e3, help="Max allowed |controlled final_psi_proj_span|")
    parser.add_argument("--psiw0-max", type=float, default=1.0e3, help="Max allowed |controlled final_psi_w0_span|")
    args = parser.parse_args()

    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
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
            controlled_row, full_row = _load_case_rows(run_log)
        except Exception as exc:  # pragma: no cover - diagnostics only
            rows.append(
                {
                    "case_name": case_name,
                    "status": "FAIL_PARSE",
                    "exit_code": "",
                    "controlled_final_psi0_span": "nan",
                    "controlled_final_psi_proj_span": "nan",
                    "controlled_final_psi_w0_span": "nan",
                    "full_final_psi0_span": "nan",
                    "ratio_full_over_controlled_psi0_from_log": "nan",
                    "error": str(exc),
                }
            )
            continue

        c_psi0 = _to_float(controlled_row.get("final_psi0_span"))
        c_proj = _to_float(controlled_row.get("final_psi_proj_span"))
        c_w0 = _to_float(controlled_row.get("final_psi_w0_span"))
        f_psi0 = _to_float(full_row.get("final_psi0_span"))
        ratio = float("nan")
        if c_psi0 != 0.0 and math.isfinite(c_psi0) and math.isfinite(f_psi0):
            ratio = f_psi0 / c_psi0

        status = _status_from_limits(c_psi0, c_proj, c_w0, args.psi0_max, args.psiproj_max, args.psiw0_max)
        exit_code = ""
        if exit_code_file.exists():
            exit_code = exit_code_file.read_text(encoding="utf-8", errors="replace").strip()

        rows.append(
            {
                "case_name": case_name,
                "status": status,
                "exit_code": exit_code,
                "controlled_final_psi0_span": f"{c_psi0:.9e}",
                "controlled_final_psi_proj_span": f"{c_proj:.9e}",
                "controlled_final_psi_w0_span": f"{c_w0:.9e}",
                "full_final_psi0_span": f"{f_psi0:.9e}",
                "ratio_full_over_controlled_psi0_from_log": f"{ratio:.9e}",
                "error": "",
            }
        )

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
        "ratio_full_over_controlled_psi0_from_log",
        "error",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    pass_count = sum(1 for row in rows if row["status"] == "PASS")
    out_md = Path(args.out_md)
    lines: List[str] = []
    lines.append("# Step-123 Controlled Bounds Summary")
    lines.append("")
    lines.append(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.")
    lines.append("")
    lines.append("Bounds:")
    lines.append(f"- `|controlled final_psi0_span| <= {args.psi0_max:.3e}`")
    lines.append(f"- `|controlled final_psi_proj_span| <= {args.psiproj_max:.3e}`")
    lines.append(f"- `|controlled final_psi_w0_span| <= {args.psiw0_max:.3e}`")
    lines.append("")
    lines.append(f"Pass count: `{pass_count}/{len(rows)}`")
    lines.append("")
    lines.append("| case | status | exit code | ctrl psi0 | ctrl psi_proj | ctrl psi_w0 | full psi0 | full/ctrl psi0 |")
    lines.append("|:---|:---|:---|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            "| {case_name} | {status} | {exit_code} | {controlled_final_psi0_span} | {controlled_final_psi_proj_span} | {controlled_final_psi_w0_span} | {full_final_psi0_span} | {ratio_full_over_controlled_psi0_from_log} |".format(
                **row
            )
        )
    lines.append("")
    out_md.write_text("\n".join(lines), encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    print(f"pass_count,{pass_count}/{len(rows)}")


if __name__ == "__main__":
    main()
