#!/usr/bin/env python3
"""Regression check for strict-1D CPAW initialization/runtime."""

import argparse
import re
import subprocess
from pathlib import Path


def run_cpaw(
    binary,
    input_path,
    workdir,
    log_path,
    nx2,
    nx3,
    tlim,
    nlim,
):
    cmd = [
        str(binary),
        "-i",
        str(input_path),
        "parthenon/output0/dt=-1",
        f"parthenon/mesh/nx2={nx2}",
        f"parthenon/mesh/nx3={nx3}",
        f"parthenon/meshblock/nx2={nx2}",
        f"parthenon/meshblock/nx3={nx3}",
        f"parthenon/time/tlim={tlim}",
        f"parthenon/time/nlim={nlim}",
    ]
    proc = subprocess.run(
        cmd,
        cwd=workdir,
        check=False,
        capture_output=True,
        text=True,
    )
    text = proc.stdout + proc.stderr
    log_path.write_text(text, encoding="utf-8")
    return proc.returncode, text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--input", default=None, help="CPAW input deck")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_cpaw_strict1d_regression",
        help="Output directory (relative to workdir if not absolute)",
    )
    parser.add_argument(
        "--nx2",
        type=int,
        default=1,
        help="Strict-1D x2 mesh size override",
    )
    parser.add_argument(
        "--nx3",
        type=int,
        default=1,
        help="Strict-1D x3 mesh size override",
    )
    parser.add_argument(
        "--tlim",
        type=float,
        default=1.0e-2,
        help="Short runtime for crash-gating",
    )
    parser.add_argument(
        "--nlim",
        type=int,
        default=50,
        help="Cycle cap for crash-gating",
    )
    parser.add_argument(
        "--min-final-cycle",
        type=int,
        default=1,
        help="Minimum final cycle expected in a successful strict-1D run",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    binary = Path(args.binary).resolve()
    input_path = Path(args.input).resolve() if args.input else repo_root / "inputs" / "cpaw.in"
    workdir = Path(args.workdir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "cpaw_strict1d.log"
    returncode, text = run_cpaw(
        binary,
        input_path,
        workdir,
        log_path,
        args.nx2,
        args.nx3,
        args.tlim,
        args.nlim,
    )

    cycles = [int(m.group(1)) for m in re.finditer(r"cycle=(\d+)", text)]
    final_cycle = max(cycles) if cycles else -1
    driver_completed = "Driver completed." in text

    failures = []
    if returncode != 0:
        failures.append("nonzero_exit")
    if not driver_completed:
        failures.append("driver_not_completed")
    if final_cycle < args.min_final_cycle:
        failures.append("insufficient_cycles")

    status = "FAIL" if failures else "PASS"
    print("case,returncode,final_cycle,driver_completed,status,failures,log")
    print(
        f"cpaw_strict1d,{returncode},{final_cycle},{int(driver_completed)},"
        f"{status},{'|'.join(failures) if failures else 'none'},{log_path}"
    )

    if failures:
        raise SystemExit("cpaw strict-1d regression failed")


if __name__ == "__main__":
    main()
