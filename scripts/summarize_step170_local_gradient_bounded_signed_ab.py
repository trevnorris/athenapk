#!/usr/bin/env python3
"""Summarize Step-170 local-gradient bounded signed operator A/B."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_hsm(script_dir: Path):
    path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("hsm", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def parse_case_rows(run_log: Path) -> tuple[dict[str, str], dict[str, str]]:
    lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()
    header = controlled = full = None
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


def sf(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return math.nan


def ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def fmt(v: float) -> str:
    if not math.isfinite(v):
        return "nan"
    return f"{v:.6e}"


def max_abs(values):
    finite = [abs(v) for v in values if math.isfinite(v)]
    return max(finite) if finite else math.nan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True)
    parser.add_argument("--baseline-case", default="baseline_quadrature")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]

    rows = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        hst = case_dir / "harris_full.out1.hst"
        if not run_log.exists() or not hst.exists():
            continue
        controlled, full = parse_case_rows(run_log)
        cols = hsm.parse_hst(hst)
        transport_vals = {
            "pi0": sf(full.get("pi0_transport_max_abs_rate", "nan")),
            "piy": sf(full.get("piy_transport_max_abs_rate", "nan")),
            "piw": sf(full.get("piw_transport_max_abs_rate", "nan")),
        }
        top_pi = max(transport_vals, key=lambda k: transport_vals[k])
        row = {
            "case_name": case_dir.name,
            "local_status": full.get("closure_local_mode0_status", ""),
            "transport_status": full.get("transport_closure_status", ""),
            "ledger_status": full.get("em_bulk_ledger_status", ""),
            "ledger_max_abs_rate": sf(full.get("em_bulk_ledger_max_abs_rate", "nan")),
            "top_pi_channel": top_pi,
            "top_pi_transport_max_abs_rate": transport_vals[top_pi],
            "psi0_ratio": ratio(sf(full.get("final_psi0_span", "nan")), sf(controlled.get("final_psi0_span", "nan"))),
            "psiw0_ratio": ratio(sf(full.get("final_psi_w0_span", "nan")), sf(controlled.get("final_psi_w0_span", "nan"))),
            "jw_ew": sf(full.get("final_jw_ew", "nan")),
            "s_leak_abs": sf(full.get("final_s_leak_abs", "nan")),
            "selected_grad_abs_max": max_abs(cols.get("m4d_response_w0_local_grad_abs", [])),
            "local_dz_abs_max": max_abs(cols.get("m4d_response_w0_local_dz_abs", [])),
        }
        rows.append(row)

    baseline = next((r for r in rows if r["case_name"] == args.baseline_case), None)
    for row in rows:
        row["ledger_over_baseline"] = ratio(row["ledger_max_abs_rate"], baseline["ledger_max_abs_rate"]) if baseline else math.nan
        row["top_pi_over_baseline"] = ratio(row["top_pi_transport_max_abs_rate"], baseline["top_pi_transport_max_abs_rate"]) if baseline else math.nan

    fieldnames = list(rows[0].keys())
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    lines = [
        "# Step-170 Local-Gradient Bounded Signed Operator A/B",
        "",
        "Compare the corrected quadrature baseline against sign-preserving local-gradient operator variants.",
        "",
        "| Case | Local | Transport | Ledger | Ledger max abs rate | Ledger / baseline | Top pi | Top pi max abs rate | Top pi / baseline | psi0 ratio | psiw0 ratio | JwEw | Sleak_abs | Selected grad abs max | Local dz abs max |",
        "| --- | --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {case_name} | {local_status} | {transport_status} | {ledger_status} | {ledger_max_abs_rate} | {ledger_over_baseline} | {top_pi_channel} | {top_pi_transport_max_abs_rate} | {top_pi_over_baseline} | {psi0_ratio} | {psiw0_ratio} | {jw_ew} | {s_leak_abs} | {selected_grad_abs_max} | {local_dz_abs_max} |".format(
                case_name=row["case_name"],
                local_status=row["local_status"],
                transport_status=row["transport_status"],
                ledger_status=row["ledger_status"],
                ledger_max_abs_rate=fmt(row["ledger_max_abs_rate"]),
                ledger_over_baseline=fmt(row["ledger_over_baseline"]),
                top_pi_channel=row["top_pi_channel"],
                top_pi_transport_max_abs_rate=fmt(row["top_pi_transport_max_abs_rate"]),
                top_pi_over_baseline=fmt(row["top_pi_over_baseline"]),
                psi0_ratio=fmt(row["psi0_ratio"]),
                psiw0_ratio=fmt(row["psiw0_ratio"]),
                jw_ew=fmt(row["jw_ew"]),
                s_leak_abs=fmt(row["s_leak_abs"]),
                selected_grad_abs_max=fmt(row["selected_grad_abs_max"]),
                local_dz_abs_max=fmt(row["local_dz_abs_max"]),
            )
        )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")


if __name__ == "__main__":
    main()
