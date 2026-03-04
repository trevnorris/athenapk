#!/usr/bin/env python3
"""Summarize transverse/mixed-channel activation probe outputs."""

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


def series_max_abs(series):
    if series is None or len(series) == 0:
        return math.nan
    return max(abs(float(v)) for v in series)


def parse_arg_from_log(log_path: Path, key: str) -> float:
    if not log_path.exists():
        return math.nan
    pattern = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    try:
        text = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return math.nan
    matches = list(pattern.finditer(text))
    if not matches:
        return math.nan
    # Use the last match so explicit per-case overrides take precedence.
    raw = matches[-1].group(1)
    low = raw.strip().lower()
    if low in {"true", "on", "yes"}:
        return 1.0
    if low in {"false", "off", "no"}:
        return 0.0
    try:
        return float(raw)
    except ValueError:
        return math.nan
    return math.nan


def gate_status(value, threshold):
    if not math.isfinite(value):
        return "N/A"
    return "PASS" if value >= threshold else "FAIL"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        default="/projects/fluid-engine/out/step106_transverse_channel_activation_*/*",
        help="Glob for output case directories",
    )
    parser.add_argument(
        "--out-csv",
        default="/projects/fluid-engine/out/step106_transverse_channel_activation_summary.csv",
        help="Path to write summary CSV",
    )
    parser.add_argument(
        "--out-md",
        default="/projects/fluid-engine/out/step106_transverse_channel_activation_summary.md",
        help="Path to write summary markdown",
    )
    parser.add_argument(
        "--min-max-ew-l1",
        type=float,
        default=2.0e-2,
        help="Activation floor for max|Ew| proxy (m4d_ew_l1)",
    )
    parser.add_argument(
        "--min-max-mixed-c2",
        type=float,
        default=1.0e-1,
        help="Activation floor for max C_a^2 proxy (m4d_mixed_c2)",
    )
    parser.add_argument(
        "--min-max-jw-l1",
        type=float,
        default=1.0e-2,
        help="Activation floor for max|Jw| proxy (m4d_jw_l1)",
    )
    parser.add_argument(
        "--min-max-jw-ew",
        type=float,
        default=1.0e-12,
        help="Activation floor for max|int Jw*Ew|",
    )
    parser.add_argument(
        "--min-max-em-leak-w",
        type=float,
        default=1.0e-3,
        help="Activation floor for max|S_EM^w proxy| (m4d_em_leak_w)",
    )
    parser.add_argument(
        "--min-max-emf-vw-c-abs",
        type=float,
        default=1.0e-4,
        help="Activation floor for max|emf_vw_c_abs|",
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

        max_abs_jw_l1 = series_max_abs(full_cols.get("m4d_jw_l1"))
        max_abs_ew_l1 = series_max_abs(full_cols.get("m4d_ew_l1"))
        max_abs_mixed_c2 = series_max_abs(full_cols.get("m4d_mixed_c2"))
        max_abs_jw_ew = series_max_abs(full_cols.get("m4d_int_jw_ew"))
        max_abs_s_leak_abs = series_max_abs(full_cols.get("m4d_int_s_leak_abs"))
        max_abs_em_leak_w = series_max_abs(full_cols.get("m4d_em_leak_w"))
        max_abs_emf_vw_c_abs = series_max_abs(full_cols.get("m4d_emf_vw_c_abs"))
        max_abs_emf_cov_vxb_abs = series_max_abs(full_cols.get("m4d_emf_cov_vxb_abs"))
        max_abs_response_w0_drive_work = series_max_abs(
            full_cols.get("m4d_response_w0_drive_work")
        )
        max_abs_response_w0_drive_leak = series_max_abs(
            full_cols.get("m4d_response_w0_drive_leak")
        )
        max_abs_response_lambda_drive_work = series_max_abs(
            full_cols.get("m4d_response_lambda_drive_work")
        )
        max_abs_response_lambda_drive_leak = series_max_abs(
            full_cols.get("m4d_response_lambda_drive_leak")
        )
        max_abs_dynamic_mix_a_abs = series_max_abs(
            full_cols.get("m4d_int_response_w0_dynamic_mix_a_abs")
        )
        max_abs_dynamic_mix_pi_abs = series_max_abs(
            full_cols.get("m4d_int_response_w0_dynamic_mix_pi_abs")
        )
        max_abs_lambda_dynamic_mix_a_abs = series_max_abs(
            full_cols.get("m4d_int_response_lambda_dynamic_mix_a_abs")
        )
        max_abs_lambda_dynamic_mix_pi_abs = series_max_abs(
            full_cols.get("m4d_int_response_lambda_dynamic_mix_pi_abs")
        )

        gate_ew_active = gate_status(max_abs_ew_l1, args.min_max_ew_l1)
        gate_c_active = gate_status(max_abs_mixed_c2, args.min_max_mixed_c2)
        gate_jw_active = gate_status(max_abs_jw_l1, args.min_max_jw_l1)
        gate_jw_ew_active = gate_status(max_abs_jw_ew, args.min_max_jw_ew)
        gate_em_leak_w_active = gate_status(max_abs_em_leak_w, args.min_max_em_leak_w)
        gate_emf_vw_c_active = gate_status(
            max_abs_emf_vw_c_abs, args.min_max_emf_vw_c_abs
        )

        log_path = d / "run.log"
        rows.append(
            {
                "case_name": d.name,
                "aw_mode1_amp": parse_arg_from_log(log_path, "problem/harris_4d/aw_mode1_amp"),
                "piw_mode1_amp": parse_arg_from_log(
                    log_path, "problem/harris_4d/piw_mode1_amp"
                ),
                "ay_mode1_amp": parse_arg_from_log(log_path, "problem/harris_4d/ay_mode1_amp"),
                "piy_mode1_amp": parse_arg_from_log(
                    log_path, "problem/harris_4d/piy_mode1_amp"
                ),
                "drift_momw_scale": parse_arg_from_log(
                    log_path, "problem/harris_4d/drift_momw_scale"
                ),
                "em_source_current_gain": parse_arg_from_log(
                    log_path, "modes4d/em_source_current_gain"
                ),
                "em_source_timelike_gain": parse_arg_from_log(
                    log_path, "modes4d/em_source_timelike_gain"
                ),
                "em_source_damping_gain": parse_arg_from_log(
                    log_path, "modes4d/em_source_damping_gain"
                ),
                "plasma_momw_source_gain": parse_arg_from_log(
                    log_path, "modes4d/plasma_momw_source_gain"
                ),
                "plasma_momw_pressure_source_gain": parse_arg_from_log(
                    log_path, "modes4d/plasma_momw_pressure_source_gain"
                ),
                "response_w0_drive_from_work_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_drive_from_work_gain"
                ),
                "response_w0_drive_from_leak_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_drive_from_leak_gain"
                ),
                "response_w0_dynamic_mixing_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_dynamic_mixing_gain"
                ),
                "response_w0_local_mode1_amp": parse_arg_from_log(
                    log_path, "modes4d/response_w0_local_mode1_amp"
                ),
                "response_w0_local_mode1_omega": parse_arg_from_log(
                    log_path, "modes4d/response_w0_local_mode1_omega"
                ),
                "response_w0_local_mode2_amp": parse_arg_from_log(
                    log_path, "modes4d/response_w0_local_mode2_amp"
                ),
                "response_w0_local_mode2_omega": parse_arg_from_log(
                    log_path, "modes4d/response_w0_local_mode2_omega"
                ),
                "response_lambda_drive_from_work_gain": parse_arg_from_log(
                    log_path, "modes4d/response_lambda_drive_from_work_gain"
                ),
                "response_lambda_drive_from_leak_gain": parse_arg_from_log(
                    log_path, "modes4d/response_lambda_drive_from_leak_gain"
                ),
                "response_lambda_dynamic_mixing_gain": parse_arg_from_log(
                    log_path, "modes4d/response_lambda_dynamic_mixing_gain"
                ),
                "full_final_psi0_span": full.get("final_psi0_span", math.nan),
                "ratio_full_over_controlled_psi0": safe_ratio(
                    full.get("final_psi0_span", math.nan),
                    controlled.get("final_psi0_span", math.nan),
                ),
                "full_final_psi_w0_minus_psi_proj_span": full.get(
                    "final_psi_w0_minus_psi_proj_span", math.nan
                ),
                "full_final_jw_l1": full.get("final_jw_l1", math.nan),
                "full_max_abs_jw_l1": max_abs_jw_l1,
                "full_final_ew_l1": full.get("final_ew_l1", math.nan),
                "full_max_abs_ew_l1": max_abs_ew_l1,
                "full_final_mixed_c2": full.get("final_mixed_c2", math.nan),
                "full_max_abs_mixed_c2": max_abs_mixed_c2,
                "full_final_jw_ew": full.get("final_jw_ew", math.nan),
                "full_max_abs_jw_ew": max_abs_jw_ew,
                "full_final_s_leak_abs": full.get("final_s_leak_abs", math.nan),
                "full_max_abs_s_leak_abs": max_abs_s_leak_abs,
                "full_final_em_leak_w": full.get("final_em_leak_w", math.nan),
                "full_max_abs_em_leak_w": max_abs_em_leak_w,
                "full_final_emf_vw_c_abs": full.get("final_emf_vw_c_abs", math.nan),
                "full_max_abs_emf_vw_c_abs": max_abs_emf_vw_c_abs,
                "full_final_emf_cov_vxb_abs": full.get("final_emf_cov_vxb_abs", math.nan),
                "full_max_abs_emf_cov_vxb_abs": max_abs_emf_cov_vxb_abs,
                "full_final_int_src_mode0_total_abs": full.get(
                    "final_int_src_mode0_total_abs", math.nan
                ),
                "full_final_response_w0_drive_work": full.get(
                    "final_response_w0_drive_work", math.nan
                ),
                "full_max_abs_response_w0_drive_work": max_abs_response_w0_drive_work,
                "full_final_response_w0_drive_leak": full.get(
                    "final_response_w0_drive_leak", math.nan
                ),
                "full_max_abs_response_w0_drive_leak": max_abs_response_w0_drive_leak,
                "full_final_response_lambda_drive_work": full.get(
                    "final_response_lambda_drive_work", math.nan
                ),
                "full_max_abs_response_lambda_drive_work": max_abs_response_lambda_drive_work,
                "full_final_response_lambda_drive_leak": full.get(
                    "final_response_lambda_drive_leak", math.nan
                ),
                "full_max_abs_response_lambda_drive_leak": max_abs_response_lambda_drive_leak,
                "full_final_int_response_w0_dynamic_mix_a_abs": full.get(
                    "final_int_response_w0_dynamic_mix_a_abs", math.nan
                ),
                "full_max_abs_int_response_w0_dynamic_mix_a_abs": max_abs_dynamic_mix_a_abs,
                "full_final_int_response_w0_dynamic_mix_pi_abs": full.get(
                    "final_int_response_w0_dynamic_mix_pi_abs", math.nan
                ),
                "full_max_abs_int_response_w0_dynamic_mix_pi_abs": max_abs_dynamic_mix_pi_abs,
                "full_final_response_w0_local_abs": full.get(
                    "final_response_w0_local_abs", math.nan
                ),
                "full_final_response_w0_local_dot_abs": full.get(
                    "final_response_w0_local_dot_abs", math.nan
                ),
                "full_final_int_response_lambda_dynamic_mix_a_abs": full.get(
                    "final_int_response_lambda_dynamic_mix_a_abs", math.nan
                ),
                "full_max_abs_int_response_lambda_dynamic_mix_a_abs": max_abs_lambda_dynamic_mix_a_abs,
                "full_final_int_response_lambda_dynamic_mix_pi_abs": full.get(
                    "final_int_response_lambda_dynamic_mix_pi_abs", math.nan
                ),
                "full_max_abs_int_response_lambda_dynamic_mix_pi_abs": max_abs_lambda_dynamic_mix_pi_abs,
                "gate_ew_active": gate_ew_active,
                "gate_c_active": gate_c_active,
                "gate_jw_active": gate_jw_active,
                "gate_jw_ew_active": gate_jw_ew_active,
                "gate_em_leak_w_active": gate_em_leak_w_active,
                "gate_emf_vw_c_active": gate_emf_vw_c_active,
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
        f.write("# Transverse Channel Activation Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        f.write("Activation floors used:\n")
        f.write(
            f"- `max|Ew_l1| >= {args.min_max_ew_l1:.3e}`\n"
            f"- `max|mixed_c2| >= {args.min_max_mixed_c2:.3e}`\n"
            f"- `max|Jw_l1| >= {args.min_max_jw_l1:.3e}`\n"
            f"- `max|int_jw_ew| >= {args.min_max_jw_ew:.3e}`\n"
            f"- `max|em_leak_w| >= {args.min_max_em_leak_w:.3e}`\n"
            f"- `max|emf_vw_c_abs| >= {args.min_max_emf_vw_c_abs:.3e}`\n\n"
        )
        if not rows:
            f.write("No valid cases found.\n")
        else:
            f.write(
                "| case | aw/piw | ay/piy | drift_momw | em_cur/time/damp | momw/press gains | w0 drv work/leak/mix/local | lambda drv work/leak/mix | Ew(final/max) | C2(final/max) | Jw(final/max) | JwEw(final/max) | S_leak(final/max) | S_EMw(final/max) | emf_vw_c(final/max) | emf_cov(final/max) | gate Ew/C/Jw/JwEw/Semw/EMF | psi_w0-psi_proj | psi0 ratio | proj/local/transport/ledger |\n"
            )
            f.write(
                "|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|:---|\n"
            )
            for r in rows:
                gate_pack = (
                    f"{r['gate_ew_active']}/"
                    f"{r['gate_c_active']}/"
                    f"{r['gate_jw_active']}/"
                    f"{r['gate_jw_ew_active']}/"
                    f"{r['gate_em_leak_w_active']}/"
                    f"{r['gate_emf_vw_c_active']}"
                )
                status_pack = (
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
                            f"{fmt(r['aw_mode1_amp'])}/{fmt(r['piw_mode1_amp'])}",
                            f"{fmt(r['ay_mode1_amp'])}/{fmt(r['piy_mode1_amp'])}",
                            fmt(r["drift_momw_scale"]),
                            f"{fmt(r['em_source_current_gain'])}/{fmt(r['em_source_timelike_gain'])}/{fmt(r['em_source_damping_gain'])}",
                            f"{fmt(r['plasma_momw_source_gain'])}/{fmt(r['plasma_momw_pressure_source_gain'])}",
                            f"{fmt(r['response_w0_drive_from_work_gain'])}/{fmt(r['response_w0_drive_from_leak_gain'])}/{fmt(r['response_w0_dynamic_mixing_gain'])}/{fmt(r['response_w0_local_mode1_amp'])}@{fmt(r['response_w0_local_mode1_omega'])}",
                            f"{fmt(r['response_lambda_drive_from_work_gain'])}/{fmt(r['response_lambda_drive_from_leak_gain'])}/{fmt(r['response_lambda_dynamic_mixing_gain'])}",
                            f"{fmt(r['full_final_ew_l1'])}/{fmt(r['full_max_abs_ew_l1'])}",
                            f"{fmt(r['full_final_mixed_c2'])}/{fmt(r['full_max_abs_mixed_c2'])}",
                            f"{fmt(r['full_final_jw_l1'])}/{fmt(r['full_max_abs_jw_l1'])}",
                            f"{fmt(r['full_final_jw_ew'])}/{fmt(r['full_max_abs_jw_ew'])}",
                            f"{fmt(r['full_final_s_leak_abs'])}/{fmt(r['full_max_abs_s_leak_abs'])}",
                            f"{fmt(r['full_final_em_leak_w'])}/{fmt(r['full_max_abs_em_leak_w'])}",
                            f"{fmt(r['full_final_emf_vw_c_abs'])}/{fmt(r['full_max_abs_emf_vw_c_abs'])}",
                            f"{fmt(r['full_final_emf_cov_vxb_abs'])}/{fmt(r['full_max_abs_emf_cov_vxb_abs'])}",
                            gate_pack,
                            fmt(r["full_final_psi_w0_minus_psi_proj_span"]),
                            fmt(r["ratio_full_over_controlled_psi0"]),
                            status_pack,
                        ]
                    )
                    + " |\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
