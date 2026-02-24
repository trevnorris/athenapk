#!/usr/bin/env python3
"""Summarize Step-116 local warp gradient-mixing probe runs."""

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


def parse_arg_from_log(log_path: Path, key: str):
    if not log_path.exists():
        return math.nan
    token = f"{key}="
    last_value = None

    def capture_last_value(line: str):
        nonlocal last_value
        start = 0
        while True:
            idx = line.find(token, start)
            if idx < 0:
                break
            tail = line[idx + len(token) :].strip()
            value = tail.split()[0].split(",")[0]
            last_value = value
            start = idx + len(token)

    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        # Prefer the command line, where overrides can appear multiple times.
        for _ in range(96):
            line = f.readline()
            if not line:
                break
            if line.startswith("command,"):
                capture_last_value(line)
                break

        # Fallback: scan an initial chunk in case command line is unavailable.
        if last_value is None:
            f.seek(0)
            for _ in range(256):
                line = f.readline()
                if not line:
                    break
                capture_last_value(line)

    if last_value is not None:
        try:
            return float(last_value)
        except ValueError:
            return last_value
    return math.nan


def ratio(num, den):
    if not (math.isfinite(num) and math.isfinite(den)):
        return math.nan
    if den == 0.0:
        return math.nan
    return num / den


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
        psi0_ratio = ratio(psi0_full, psi0_ctrl)

        grad_mix_a = full.get("final_int_response_w0_local_gradient_mix_a_abs", math.nan)
        grad_mix_pi = full.get("final_int_response_w0_local_gradient_mix_pi_abs", math.nan)
        total_mix_a = full.get("final_int_response_w0_dynamic_mix_a_abs", math.nan)
        total_mix_pi = full.get("final_int_response_w0_dynamic_mix_pi_abs", math.nan)

        rows.append(
            {
                "case_name": case_dir.name,
                "mix_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_dynamic_mixing_gain"
                ),
                "grad_mix_enable": parse_arg_from_log(
                    log_path, "modes4d/response_w0_local_gradient_mixing_enable"
                ),
                "grad_mix_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_local_gradient_mixing_gain"
                ),
                "omega1": parse_arg_from_log(log_path, "modes4d/response_w0_local_mode1_omega"),
                "omega2": parse_arg_from_log(log_path, "modes4d/response_w0_local_mode2_omega"),
                "full_final_response_w0_local_abs": full.get(
                    "final_response_w0_local_abs", math.nan
                ),
                "full_final_response_w0_local_dot_abs": full.get(
                    "final_response_w0_local_dot_abs", math.nan
                ),
                "full_final_response_w0_local_grad_abs": full.get(
                    "final_response_w0_local_grad_abs", math.nan
                ),
                "full_final_int_response_w0_dynamic_mix_a_abs": total_mix_a,
                "full_final_int_response_w0_dynamic_mix_pi_abs": total_mix_pi,
                "full_final_int_response_w0_local_gradient_mix_a_abs": grad_mix_a,
                "full_final_int_response_w0_local_gradient_mix_pi_abs": grad_mix_pi,
                "grad_fraction_mix_a": ratio(grad_mix_a, total_mix_a),
                "grad_fraction_mix_pi": ratio(grad_mix_pi, total_mix_pi),
                "psi0_ratio_full_over_controlled": psi0_ratio,
                "full_final_psi_w0_minus_psi_proj_span": full.get(
                    "final_psi_w0_minus_psi_proj_span", math.nan
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
        f.write("# Local Warp Gradient-Mixing Probe Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid case outputs found.\n")
        else:
            f.write(
                "| case | mix_gain | grad_mix enable/gain | omega1/omega2 | local |w0|/|w0dot|/|grad w0| | mix_a total/grad/fraction | mix_pi total/grad/fraction | psi0 ratio | psi_w0-psi_proj | jw_ew | s_leak_abs | projected/local/transport/ledger |\n"
            )
            f.write(
                "|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|\n"
            )
            for row in rows:
                status = (
                    f"{row['full_closure_projected_status']}/"
                    f"{row['full_closure_local_mode0_status']}/"
                    f"{row['full_transport_closure_status']}/"
                    f"{row['full_em_bulk_ledger_status']}"
                )
                f.write(
                    "| "
                    + " | ".join(
                        [
                            row["case_name"],
                            fmt(row["mix_gain"]),
                            f"{row['grad_mix_enable']}/{fmt(row['grad_mix_gain'])}",
                            f"{fmt(row['omega1'])}/{fmt(row['omega2'])}",
                            f"{fmt(row['full_final_response_w0_local_abs'])}/"
                            f"{fmt(row['full_final_response_w0_local_dot_abs'])}/"
                            f"{fmt(row['full_final_response_w0_local_grad_abs'])}",
                            f"{fmt(row['full_final_int_response_w0_dynamic_mix_a_abs'])}/"
                            f"{fmt(row['full_final_int_response_w0_local_gradient_mix_a_abs'])}/"
                            f"{fmt(row['grad_fraction_mix_a'])}",
                            f"{fmt(row['full_final_int_response_w0_dynamic_mix_pi_abs'])}/"
                            f"{fmt(row['full_final_int_response_w0_local_gradient_mix_pi_abs'])}/"
                            f"{fmt(row['grad_fraction_mix_pi'])}",
                            fmt(row["psi0_ratio_full_over_controlled"]),
                            fmt(row["full_final_psi_w0_minus_psi_proj_span"]),
                            fmt(row["full_final_jw_ew"]),
                            fmt(row["full_final_s_leak_abs"]),
                            status,
                        ]
                    )
                    + " |\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
