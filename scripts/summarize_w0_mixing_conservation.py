#!/usr/bin/env python3
"""Summarize Step-110 pure w0-dynamic-mixing conservation A/B runs."""

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_harris_scan_matrix_module(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def parse_arg_from_log(log_path: Path, key: str) -> float:
    if not log_path.exists():
        return math.nan
    token = f"{key}="
    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        for _ in range(32):
            line = f.readline()
            if not line:
                break
            idx = line.find(token)
            if idx < 0:
                continue
            tail = line[idx + len(token) :].strip()
            val = tail.split()[0].split(",")[0]
            try:
                return float(val)
            except ValueError:
                return math.nan
    return math.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Glob for case directories")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--out-md", required=True, help="Output markdown path")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    rows = []
    for case_dir in sorted(Path(p) for p in glob.glob(args.glob)):
        controlled_hst = case_dir / "harris_controlled.out1.hst"
        full_hst = case_dir / "harris_full.out1.hst"
        if not controlled_hst.exists() or not full_hst.exists():
            continue

        controlled_cols = hsm.parse_hst(controlled_hst)
        full_cols = hsm.parse_hst(full_hst)
        controlled = hsm.analyze_case(
            "controlled",
            controlled_cols,
            5.0e-2,
            1.0e-8,
            1.0e-8,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )
        full = hsm.analyze_case(
            "full",
            full_cols,
            5.0e-2,
            1.0e-8,
            1.0e-8,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )
        log_path = case_dir / "run.log"
        psi0_ctrl = controlled.get("final_psi0_span", math.nan)
        psi0_full = full.get("final_psi0_span", math.nan)
        psi0_ratio = (
            psi0_full / psi0_ctrl
            if math.isfinite(psi0_full) and math.isfinite(psi0_ctrl) and abs(psi0_ctrl) > 0.0
            else math.nan
        )
        rows.append(
            {
                "case_name": case_dir.name,
                "response_w0_dot_init": parse_arg_from_log(
                    log_path, "modes4d/response_w0_dot_init"
                ),
                "response_w0_dynamic_mixing_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_dynamic_mixing_gain"
                ),
                "psi0_ratio_full_over_controlled": psi0_ratio,
                "full_final_psi_w0_minus_psi_proj_span": full.get(
                    "final_psi_w0_minus_psi_proj_span", math.nan
                ),
                "full_final_em_u_bulk": full.get("final_em_u_bulk", math.nan),
                "full_final_em_u_resolved": full.get("final_em_u_resolved", math.nan),
                "full_final_em_u_sub": full.get("final_em_u_sub", math.nan),
                "full_final_int_response_w0_dynamic_mix_a_abs": full.get(
                    "final_int_response_w0_dynamic_mix_a_abs", math.nan
                ),
                "full_final_int_response_w0_dynamic_mix_pi_abs": full.get(
                    "final_int_response_w0_dynamic_mix_pi_abs", math.nan
                ),
                "full_final_jw_ew": full.get("final_jw_ew", math.nan),
                "full_final_s_leak_abs": full.get("final_s_leak_abs", math.nan),
                "full_closure_projected_status": full.get("closure_projected_status", "N/A"),
                "full_closure_local_mode0_status": full.get(
                    "closure_local_mode0_status", "N/A"
                ),
                "full_transport_closure_status": full.get("transport_closure_status", "N/A"),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
            }
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        else:
            f.write("")

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# W0 Dynamic Mixing Conservation A/B Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid case outputs found.\n")
        else:
            f.write(
                "| case | w0_dot_init | mix_gain | psi0 ratio | psi_w0-psi_proj | em_u_bulk | em_u_resolved | em_u_sub | mix_a_abs | mix_pi_abs | jw_ew | s_leak_abs | projected/local/transport/ledger |\n"
            )
            f.write(
                "|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|\n"
            )
            for r in rows:
                status = (
                    f"{r['full_closure_projected_status']}/"
                    f"{r['full_closure_local_mode0_status']}/"
                    f"{r['full_transport_closure_status']}/"
                    f"{r['full_em_bulk_ledger_status']}"
                )
                f.write(
                    "| "
                    + " | ".join(
                        [
                            r["case_name"],
                            fmt(r["response_w0_dot_init"]),
                            fmt(r["response_w0_dynamic_mixing_gain"]),
                            fmt(r["psi0_ratio_full_over_controlled"]),
                            fmt(r["full_final_psi_w0_minus_psi_proj_span"]),
                            fmt(r["full_final_em_u_bulk"]),
                            fmt(r["full_final_em_u_resolved"]),
                            fmt(r["full_final_em_u_sub"]),
                            fmt(r["full_final_int_response_w0_dynamic_mix_a_abs"]),
                            fmt(r["full_final_int_response_w0_dynamic_mix_pi_abs"]),
                            fmt(r["full_final_jw_ew"]),
                            fmt(r["full_final_s_leak_abs"]),
                            status,
                        ]
                    )
                    + " |\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
