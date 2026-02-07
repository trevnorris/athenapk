#!/usr/bin/env python3
"""Run the tuned 4D Harris regression with closure + activity checks."""

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument(
        "--workdir",
        default=".",
        help="Run directory used by harris_scan_matrix.py",
    )
    parser.add_argument(
        "--output-dir",
        default="harris_regression_outputs",
        help="Output directory for copied history files (relative to workdir)",
    )
    parser.add_argument(
        "--min-jw-ew-abs",
        type=float,
        default=1.0e-15,
        help="Minimum |final_jw_ew| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-s-leak-abs",
        type=float,
        default=1.0e-14,
        help="Minimum final_s_leak_abs for full-case activity PASS",
    )
    parser.add_argument(
        "--min-mixed-ew2",
        type=float,
        default=1.0e-3,
        help="Minimum final_mixed_ew2 for full-case activity PASS",
    )
    parser.add_argument(
        "--min-jw-mode-l2-1",
        type=float,
        default=1.0e-30,
        help="Minimum final_jw_mode_l2_1 for full-case activity PASS",
    )
    parser.add_argument(
        "--min-abs-em-leak-w",
        type=float,
        default=1.0e-4,
        help="Minimum |final_em_leak_w| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-abs-helicity-sub",
        type=float,
        default=1.0e-12,
        help="Minimum |final_helicity_sub| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-abs-edotb-sub",
        type=float,
        default=1.0e-11,
        help="Minimum |final_edotb_sub| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-src-em-laplacian-abs",
        type=float,
        default=1.0e-2,
        help="Minimum |final_int_src_em_laplacian_abs| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-src-em-mass-abs",
        type=float,
        default=1.0e-4,
        help="Minimum |final_int_src_em_mass_abs| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-src-em-current-abs",
        type=float,
        default=1.0e-13,
        help="Minimum |final_int_src_em_current_abs| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-src-em-spatial-mixed-abs",
        type=float,
        default=1.0e-2,
        help="Minimum |final_int_src_em_spatial_mixed_abs| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-src-em-timelike-abs",
        type=float,
        default=1.0e-3,
        help="Minimum |final_int_src_em_timelike_abs| for full-case activity PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-s-leak-abs",
        type=float,
        default=9.0e-1,
        help="Minimum |corr(psi0_span, int_s_leak_abs)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-jw-ew",
        type=float,
        default=9.0e-1,
        help="Minimum |corr(psi0_span, int_jw_ew)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-mixed-ew2",
        type=float,
        default=9.0e-1,
        help="Minimum |corr(psi0_span, mixed_ew2)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-mixed-c2",
        type=float,
        default=9.0e-1,
        help="Minimum |corr(psi0_span, mixed_c2)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-jw-mode-l2-1",
        type=float,
        default=9.0e-1,
        help="Minimum |corr(psi0_span, jw_mode_l2_1)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-em-leak-w",
        type=float,
        default=3.0e-1,
        help="Minimum |corr(psi0_span, em_leak_w)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-helicity-sub",
        type=float,
        default=7.0e-1,
        help="Minimum |corr(psi0_span, helicity_sub)| for full-case PASS",
    )
    parser.add_argument(
        "--min-abs-corr-psi0-edotb-sub",
        type=float,
        default=7.0e-1,
        help="Minimum |corr(psi0_span, edotb_sub)| for full-case PASS",
    )
    parser.add_argument(
        "--scan-athena-arg",
        action="append",
        default=[],
        help="Additional athena override forwarded to both controlled/full scan runs",
    )
    parser.add_argument(
        "--scan-controlled-arg",
        action="append",
        default=[],
        help="Additional athena override forwarded only to controlled scan run",
    )
    parser.add_argument(
        "--scan-full-arg",
        action="append",
        default=[],
        help="Additional athena override forwarded only to full scan run",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    scan_script = repo_root / "scripts" / "harris_scan_matrix.py"
    controlled_input = repo_root / "inputs" / "harris_4d_controlled.in"
    full_input = repo_root / "inputs" / "harris_4d_full.in"

    workdir = Path(args.workdir).resolve()
    output_dir = workdir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(scan_script),
        "--binary",
        str(Path(args.binary).resolve()),
        "--workdir",
        str(workdir),
        "--controlled-input",
        str(controlled_input),
        "--full-input",
        str(full_input),
        "--output-dir",
        str(output_dir),
        "--fail-on-check",
        "--check-transport-closure",
        "--full-min-jw-ew-abs",
        str(args.min_jw_ew_abs),
        "--full-min-s-leak-abs",
        str(args.min_s_leak_abs),
        "--full-min-mixed-ew2",
        str(args.min_mixed_ew2),
        "--full-min-jw-mode-l2-1",
        str(args.min_jw_mode_l2_1),
        "--full-min-abs-em-leak-w",
        str(args.min_abs_em_leak_w),
        "--full-min-abs-helicity-sub",
        str(args.min_abs_helicity_sub),
        "--full-min-abs-edotb-sub",
        str(args.min_abs_edotb_sub),
        "--full-min-src-em-laplacian-abs",
        str(args.min_src_em_laplacian_abs),
        "--full-min-src-em-mass-abs",
        str(args.min_src_em_mass_abs),
        "--full-min-src-em-current-abs",
        str(args.min_src_em_current_abs),
        "--full-min-src-em-spatial-mixed-abs",
        str(args.min_src_em_spatial_mixed_abs),
        "--full-min-src-em-timelike-abs",
        str(args.min_src_em_timelike_abs),
        "--full-min-abs-corr-psi0-s-leak-abs",
        str(args.min_abs_corr_psi0_s_leak_abs),
        "--full-min-abs-corr-psi0-jw-ew",
        str(args.min_abs_corr_psi0_jw_ew),
        "--full-min-abs-corr-psi0-mixed-ew2",
        str(args.min_abs_corr_psi0_mixed_ew2),
        "--full-min-abs-corr-psi0-mixed-c2",
        str(args.min_abs_corr_psi0_mixed_c2),
        "--full-min-abs-corr-psi0-jw-mode-l2-1",
        str(args.min_abs_corr_psi0_jw_mode_l2_1),
        "--full-min-abs-corr-psi0-em-leak-w",
        str(args.min_abs_corr_psi0_em_leak_w),
        "--full-min-abs-corr-psi0-helicity-sub",
        str(args.min_abs_corr_psi0_helicity_sub),
        "--full-min-abs-corr-psi0-edotb-sub",
        str(args.min_abs_corr_psi0_edotb_sub),
    ]
    for arg in args.scan_athena_arg:
        cmd.extend(["--athena-arg", str(arg)])
    for arg in args.scan_controlled_arg:
        cmd.extend(["--controlled-arg", str(arg)])
    for arg in args.scan_full_arg:
        cmd.extend(["--full-arg", str(arg)])
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
