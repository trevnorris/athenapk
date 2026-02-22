#!/usr/bin/env python3
"""Summarize projected-continuity closure diagnostics for scan case directories."""

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


def analyze_case_dir(case_dir: Path, hsm):
    controlled_hst = case_dir / "harris_controlled.out1.hst"
    full_hst = case_dir / "harris_full.out1.hst"
    if not controlled_hst.exists() or not full_hst.exists():
        return None

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

    return {
        "case_dir": str(case_dir),
        "case_name": case_dir.name,
        "controlled_closure_projected_status": controlled.get(
            "closure_projected_status", "N/A"
        ),
        "controlled_closure_local_mode0_status": controlled.get(
            "closure_local_mode0_status", "N/A"
        ),
        "controlled_closure_projected_max_abs_rate": controlled.get(
            "closure_projected_max_abs_rate", math.nan
        ),
        "controlled_closure_local_mode0_max_abs_rate": controlled.get(
            "closure_local_mode0_max_abs_rate", math.nan
        ),
        "full_closure_projected_status": full.get("closure_projected_status", "N/A"),
        "full_closure_local_mode0_status": full.get("closure_local_mode0_status", "N/A"),
        "full_closure_projected_max_abs_rate": full.get(
            "closure_projected_max_abs_rate", math.nan
        ),
        "full_closure_projected_final_norm": full.get(
            "closure_projected_final_norm", math.nan
        ),
        "full_closure_projected_dcharge_final": full.get(
            "closure_projected_dcharge_final", math.nan
        ),
        "full_closure_projected_dsleak_final": full.get(
            "closure_projected_dsleak_final", math.nan
        ),
        "full_closure_projected_ddivj_final": full.get(
            "closure_projected_ddivj_final", math.nan
        ),
        "full_closure_projected_final_residual": full.get(
            "closure_projected_final_residual", math.nan
        ),
        "full_closure_local_mode0_max_abs_rate": full.get(
            "closure_local_mode0_max_abs_rate", math.nan
        ),
        "full_transport_closure_status": full.get("transport_closure_status", "N/A"),
        "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
        "full_final_psi0_span": full.get("final_psi0_span", math.nan),
        "full_final_psi_proj_span": full.get("final_psi_proj_span", math.nan),
        "full_final_psi_w0_span": full.get("final_psi_w0_span", math.nan),
        "full_final_psi_w0_minus_psi_proj_span": full.get(
            "final_psi_w0_minus_psi_proj_span", math.nan
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        required=True,
        help="Glob for scan case directories (each with harris_controlled/full .hst)",
    )
    parser.add_argument("--out-csv", required=True, help="Path to write summary CSV")
    parser.add_argument("--out-md", required=True, help="Path to write summary markdown")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)
    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob))]

    rows = []
    for case_dir in case_dirs:
        row = analyze_case_dir(case_dir, hsm)
        if row is not None:
            rows.append(row)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Projected Closure Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid cases found.\n")
        else:
            f.write(
                "| case | full projected | full local mode0 | projected max rate | projected final norm | local mode0 max rate | transport | ledger | psi_w0-psi_proj |\n"
            )
            f.write(
                "|:---|:---|:---|---:|---:|---:|:---|:---|---:|\n"
            )
            for r in rows:
                f.write(
                    "| "
                    + " | ".join(
                        [
                            r["case_name"],
                            r["full_closure_projected_status"],
                            r["full_closure_local_mode0_status"],
                            fmt(r["full_closure_projected_max_abs_rate"]),
                            fmt(r["full_closure_projected_final_norm"]),
                            fmt(r["full_closure_local_mode0_max_abs_rate"]),
                            r["full_transport_closure_status"],
                            r["full_em_bulk_ledger_status"],
                            fmt(r["full_final_psi_w0_minus_psi_proj_span"]),
                        ]
                    )
                    + " |\n"
                )
            f.write("\n")
            f.write("## Final Projected Continuity Components (full case)\n\n")
            f.write(
                "| case | dcharge_dt | dsleak_dt | ddivj_dt | residual |\n"
            )
            f.write("|:---|---:|---:|---:|---:|\n")
            for r in rows:
                f.write(
                    "| "
                    + " | ".join(
                        [
                            r["case_name"],
                            fmt(r["full_closure_projected_dcharge_final"]),
                            fmt(r["full_closure_projected_dsleak_final"]),
                            fmt(r["full_closure_projected_ddivj_final"]),
                            fmt(r["full_closure_projected_final_residual"]),
                        ]
                    )
                    + " |\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
