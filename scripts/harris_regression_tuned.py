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
        "--full-min-jw-ew-abs",
        str(args.min_jw_ew_abs),
        "--full-min-s-leak-abs",
        str(args.min_s_leak_abs),
        "--full-min-mixed-ew2",
        str(args.min_mixed_ew2),
        "--full-min-jw-mode-l2-1",
        str(args.min_jw_mode_l2_1),
    ]
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
