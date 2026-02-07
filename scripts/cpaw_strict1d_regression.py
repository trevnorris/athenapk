#!/usr/bin/env python3
"""Regression check for strict-1D CPAW initialization/runtime and parity metrics."""

import argparse
import math
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


def parse_cpaw_errors(path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        rows.append(stripped.split())
    if not rows:
        raise RuntimeError(f"No data rows in {path}")
    vals = rows[-1]
    if len(vals) < 13:
        raise RuntimeError(f"Unexpected cpaw-errors format in {path}")
    return {
        "Nx1": int(vals[0]),
        "Nx2": int(vals[1]),
        "Nx3": int(vals[2]),
        "Ncycle": int(vals[3]),
        "RMS": float(vals[4]),
        "d": float(vals[5]),
        "M1": float(vals[6]),
        "M2": float(vals[7]),
        "M3": float(vals[8]),
        "E": float(vals[9]),
        "B1c": float(vals[10]),
        "B2c": float(vals[11]),
        "B3c": float(vals[12]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--input", default=None, help="CPAW input deck")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="outputs",
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
    parser.add_argument(
        "--max-rms-error",
        type=float,
        default=1.0e-2,
        help="Maximum allowed final RMS CPAW error from cpaw-errors.dat",
    )
    parser.add_argument(
        "--min-transverse-error",
        type=float,
        default=1.0e-6,
        help="Minimum required final B2/B3 error magnitude to confirm nontrivial wave activity",
    )
    parser.add_argument(
        "--max-transverse-rel-diff",
        type=float,
        default=1.0e-2,
        help="Maximum allowed relative mismatch between final B2 and B3 errors",
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
    cpaw_errors_path = workdir / "cpaw-errors.dat"
    try:
        cpaw_errors_path.unlink()
    except FileNotFoundError:
        pass
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

    cpaw_errors_present = cpaw_errors_path.exists()
    parsed = None
    if cpaw_errors_present:
        parsed = parse_cpaw_errors(cpaw_errors_path)
    rms_error = math.nan if parsed is None else parsed["RMS"]
    b2_error = math.nan if parsed is None else abs(parsed["B2c"])
    b3_error = math.nan if parsed is None else abs(parsed["B3c"])
    b23_rel_diff = math.nan
    if parsed is not None:
        denom = max(b2_error, b3_error, 1.0e-30)
        b23_rel_diff = abs(b2_error - b3_error) / denom

    failures = []
    if returncode != 0:
        failures.append("nonzero_exit")
    if not driver_completed:
        failures.append("driver_not_completed")
    if final_cycle < args.min_final_cycle:
        failures.append("insufficient_cycles")
    if not cpaw_errors_present:
        failures.append("missing_cpaw_errors")
    if parsed is not None:
        if parsed["Nx2"] != args.nx2 or parsed["Nx3"] != args.nx3:
            failures.append("mesh_override_mismatch")
        if parsed["Ncycle"] < args.min_final_cycle:
            failures.append("cpaw_errors_insufficient_cycle")
        for k in ("RMS", "d", "M1", "M2", "M3", "E", "B1c", "B2c", "B3c"):
            if not math.isfinite(parsed[k]):
                failures.append(f"nonfinite_{k.lower()}")
        if math.isfinite(rms_error) and rms_error > args.max_rms_error:
            failures.append("rms_error_high")
        if math.isfinite(b2_error) and b2_error < args.min_transverse_error:
            failures.append("b2_error_too_small")
        if math.isfinite(b3_error) and b3_error < args.min_transverse_error:
            failures.append("b3_error_too_small")
        if math.isfinite(b23_rel_diff) and b23_rel_diff > args.max_transverse_rel_diff:
            failures.append("b2_b3_mismatch")

    status = "FAIL" if failures else "PASS"
    print(
        "case,returncode,final_cycle,driver_completed,cpaw_errors_present,rms_error,"
        "b2_error,b3_error,b23_rel_diff,status,failures,log"
    )
    print(
        f"cpaw_strict1d,{returncode},{final_cycle},{int(driver_completed)},"
        f"{int(cpaw_errors_present)},{rms_error:.6e},{b2_error:.6e},{b3_error:.6e},"
        f"{b23_rel_diff:.6e},{status},{'|'.join(failures) if failures else 'none'},{log_path}"
    )

    if failures:
        raise SystemExit("cpaw strict-1d regression failed")


if __name__ == "__main__":
    main()
