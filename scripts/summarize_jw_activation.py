#!/usr/bin/env python3
"""Summarize Jw-channel activation probe outputs."""

import argparse
import csv
import glob
import importlib.util
import math
import re
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


def safe_ratio(numer, denom):
    if (
        numer is None
        or denom is None
        or not math.isfinite(numer)
        or not math.isfinite(denom)
        or abs(denom) == 0.0
    ):
        return math.nan
    return numer / denom


def series_final(series):
    if series is None or len(series) == 0:
        return math.nan
    return float(series[-1])


def series_max_abs(series):
    if series is None or len(series) == 0:
        return math.nan
    return max(abs(float(v)) for v in series)


def parse_arg_from_log(log_path: Path, key: str) -> float:
    if not log_path.exists():
        return math.nan
    pattern = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    try:
        with log_path.open("r", encoding="utf-8", errors="ignore") as f:
            for _ in range(10):
                line = f.readline()
                if not line:
                    break
                m = pattern.search(line)
                if m is not None:
                    try:
                        return float(m.group(1))
                    except ValueError:
                        return math.nan
    except OSError:
        return math.nan
    return math.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        default="/projects/fluid-engine/out/step104_jw_channel_activation_*/*",
        help="Glob for output case directories",
    )
    parser.add_argument(
        "--out-csv",
        default="/projects/fluid-engine/out/step104_jw_channel_activation_summary.csv",
        help="Path to write summary CSV",
    )
    parser.add_argument(
        "--out-md",
        default="/projects/fluid-engine/out/step104_jw_channel_activation_summary.md",
        help="Path to write summary markdown",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
    if not dirs:
        raise SystemExit(f"No directories matched: {args.glob}")

    rows = []
    for d in dirs:
        controlled_hst = d / "harris_controlled.out1.hst"
        full_hst = d / "harris_full.out1.hst"
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

        rows.append(
            {
                "case_name": d.name,
                "drift_momw_scale": parse_arg_from_log(
                    d / "run.log", "problem/harris_4d/drift_momw_scale"
                ),
                "full_final_jw_l1": full.get("final_jw_l1", math.nan),
                "full_max_abs_jw_l1": series_max_abs(full_cols.get("m4d_jw_l1")),
                "full_final_jw_l2": full.get("final_jw_l2", math.nan),
                "full_final_ew_l1": full.get("final_ew_l1", math.nan),
                "full_final_ew_l2": full.get("final_ew_l2", math.nan),
                "full_final_jw_ew": full.get("final_jw_ew", math.nan),
                "full_max_abs_jw_ew": series_max_abs(full_cols.get("m4d_int_jw_ew")),
                "full_final_s_leak_abs": full.get("final_s_leak_abs", math.nan),
                "full_max_abs_s_leak_abs": series_max_abs(full_cols.get("m4d_int_s_leak_abs")),
                "full_final_emf_vw_c_abs": full.get("final_emf_vw_c_abs", math.nan),
                "full_final_emf_cov_vxb_abs": full.get("final_emf_cov_vxb_abs", math.nan),
                "full_final_momw_mode0": series_final(full_cols.get("m4d_momw_mode_0")),
                "full_max_abs_momw_mode0": series_max_abs(full_cols.get("m4d_momw_mode_0")),
                "full_final_psi0_span": full.get("final_psi0_span", math.nan),
                "full_final_psi_proj_span": full.get("final_psi_proj_span", math.nan),
                "full_final_psi_w0_span": full.get("final_psi_w0_span", math.nan),
                "full_final_psi_w0_minus_psi_proj_span": full.get(
                    "final_psi_w0_minus_psi_proj_span", math.nan
                ),
                "ratio_full_over_controlled_psi0": safe_ratio(
                    full.get("final_psi0_span", math.nan),
                    controlled.get("final_psi0_span", math.nan),
                ),
                "full_final_int_src_mode0_total_abs": full.get(
                    "final_int_src_mode0_total_abs", math.nan
                ),
                "full_final_mixed_ew2": full.get("final_mixed_ew2", math.nan),
                "full_closure_projected_status": full.get("closure_projected_status", "N/A"),
                "full_closure_local_mode0_status": full.get(
                    "closure_local_mode0_status", "N/A"
                ),
                "full_transport_closure_status": full.get("transport_closure_status", "N/A"),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
                "full_topology_status": full.get("topology_status", "N/A"),
                "full_mechanism_status": full.get("mechanism_status", "N/A"),
            }
        )

    rows.sort(key=lambda r: r["case_name"])

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Jw Activation Probe Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid cases found.\n")
        else:
            f.write(
                "| case | drift_momw | jw_l1(final/max) | jw_ew(final/max) | s_leak_abs(final/max) | momw0(final/max) | emf_vw_c | emf_cov_vxb | psi_w0-psi_proj | psi0 ratio | src_mode0_total | closure_proj | closure_local | transport | ledger | topology | mechanism |\n"
            )
            f.write(
                "|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|:---|:---|:---|\n"
            )
            for r in rows:
                f.write(
                    "| "
                    + " | ".join(
                        [
                            r["case_name"],
                            fmt(r["drift_momw_scale"]),
                            f"{fmt(r['full_final_jw_l1'])}/{fmt(r['full_max_abs_jw_l1'])}",
                            f"{fmt(r['full_final_jw_ew'])}/{fmt(r['full_max_abs_jw_ew'])}",
                            f"{fmt(r['full_final_s_leak_abs'])}/{fmt(r['full_max_abs_s_leak_abs'])}",
                            f"{fmt(r['full_final_momw_mode0'])}/{fmt(r['full_max_abs_momw_mode0'])}",
                            fmt(r["full_final_emf_vw_c_abs"]),
                            fmt(r["full_final_emf_cov_vxb_abs"]),
                            fmt(r["full_final_psi_w0_minus_psi_proj_span"]),
                            fmt(r["ratio_full_over_controlled_psi0"]),
                            fmt(r["full_final_int_src_mode0_total_abs"]),
                            str(r["full_closure_projected_status"]),
                            str(r["full_closure_local_mode0_status"]),
                            str(r["full_transport_closure_status"]),
                            str(r["full_em_bulk_ledger_status"]),
                            str(r["full_topology_status"]),
                            str(r["full_mechanism_status"]),
                        ]
                    )
                    + " |\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
