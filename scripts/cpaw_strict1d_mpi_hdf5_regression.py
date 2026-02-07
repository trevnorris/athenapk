#!/usr/bin/env python3
"""MPI+HDF5 strict-1D CPAW regression with quantitative phdf checks."""

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np


def latest_phdf(workdir):
    candidates = sorted(
        workdir.glob("parthenon.out0*.phdf"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        return None
    return candidates[-1]


def run_cpaw_mpi(
    mpirun,
    ranks,
    binary,
    input_path,
    workdir,
    log_path,
    nx2,
    nx3,
    tlim,
    nlim,
    output_dt,
):
    cmd = [
        str(mpirun),
        "-n",
        str(ranks),
        str(binary),
        "-i",
        str(input_path),
        f"parthenon/time/tlim={tlim}",
        f"parthenon/time/nlim={nlim}",
        f"parthenon/output0/dt={output_dt}",
        f"parthenon/mesh/nx2={nx2}",
        f"parthenon/mesh/nx3={nx3}",
        f"parthenon/meshblock/nx2={nx2}",
        f"parthenon/meshblock/nx3={nx3}",
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


def parse_cons_metrics(phdf_path, b2_comp, b3_comp):
    with h5py.File(phdf_path, "r") as f:
        if "cons" not in f:
            raise RuntimeError(f"'cons' dataset not found in {phdf_path}")
        cons = np.asarray(f["cons"])
        if cons.ndim != 5:
            raise RuntimeError(f"Expected 5D cons dataset, got shape {cons.shape}")
        ncomp = cons.shape[1]
        if b2_comp < 0 or b2_comp >= ncomp or b3_comp < 0 or b3_comp >= ncomp:
            raise RuntimeError(
                f"Requested B component indices ({b2_comp},{b3_comp}) out of bounds for {ncomp} comps"
            )
        if "LogicalLocations" not in f:
            raise RuntimeError("Missing LogicalLocations dataset in phdf output")
        logical = np.asarray(f["LogicalLocations"])
        if logical.shape[1] != 3:
            raise RuntimeError(f"Unexpected LogicalLocations shape: {logical.shape}")

        finite_cons = bool(np.all(np.isfinite(cons)))
        b2 = cons[:, b2_comp, :, :, :]
        b3 = cons[:, b3_comp, :, :, :]
        b2_mean_abs = float(np.mean(np.abs(b2)))
        b3_mean_abs = float(np.mean(np.abs(b3)))
        b2_max_abs = float(np.max(np.abs(b2)))
        b3_max_abs = float(np.max(np.abs(b3)))
        denom = max(b2_mean_abs, b3_mean_abs, 1.0e-30)
        b23_rel_diff = abs(b2_mean_abs - b3_mean_abs) / denom

        # Strict-1D geometric checks from block-local shape + logical placement.
        strict_shape = (cons.shape[2] == 1 and cons.shape[3] == 1)
        logical_yz_single = bool(np.max(logical[:, 1]) == 0 and np.max(logical[:, 2]) == 0)

    return {
        "finite_cons": finite_cons,
        "b2_mean_abs": b2_mean_abs,
        "b3_mean_abs": b3_mean_abs,
        "b2_max_abs": b2_max_abs,
        "b3_max_abs": b3_max_abs,
        "b23_rel_diff": b23_rel_diff,
        "strict_shape": strict_shape,
        "logical_yz_single": logical_yz_single,
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
    parser.add_argument("--mpirun", default="mpirun", help="MPI launcher executable")
    parser.add_argument("--ranks", type=int, default=2, help="MPI ranks")
    parser.add_argument("--nx2", type=int, default=1, help="Strict-1D x2 mesh override")
    parser.add_argument("--nx3", type=int, default=1, help="Strict-1D x3 mesh override")
    parser.add_argument(
        "--tlim", type=float, default=1.0e-2, help="Short runtime for crash/metric gate"
    )
    parser.add_argument(
        "--nlim", type=int, default=50, help="Cycle cap for crash/metric gate"
    )
    parser.add_argument(
        "--output-dt",
        type=float,
        default=5.0e-3,
        help="HDF5 output cadence override; should be <= tlim",
    )
    parser.add_argument(
        "--min-final-cycle",
        type=int,
        default=1,
        help="Minimum expected final cycle",
    )
    parser.add_argument(
        "--min-transverse-mean-abs",
        type=float,
        default=1.0e-4,
        help="Minimum required mean(|B2|) and mean(|B3|) in final phdf",
    )
    parser.add_argument(
        "--max-transverse-rel-diff",
        type=float,
        default=1.0e-2,
        help="Maximum allowed relative mismatch between mean(|B2|) and mean(|B3|)",
    )
    parser.add_argument(
        "--b2-comp",
        type=int,
        default=6,
        help="Component index for B2 in cons dataset",
    )
    parser.add_argument(
        "--b3-comp",
        type=int,
        default=7,
        help="Component index for B3 in cons dataset",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    binary = Path(args.binary).resolve()
    input_path = Path(args.input).resolve() if args.input else repo_root / "inputs" / "cpaw.in"
    mpirun = Path(shutil.which(args.mpirun) or args.mpirun)
    workdir = Path(args.workdir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for p in workdir.glob("parthenon.out0*.phdf"):
        p.unlink()
    for p in workdir.glob("parthenon.out*.hst"):
        p.unlink()

    log_path = output_dir / "cpaw_strict1d_mpi_hdf5.log"
    returncode, text = run_cpaw_mpi(
        mpirun=mpirun,
        ranks=args.ranks,
        binary=binary,
        input_path=input_path,
        workdir=workdir,
        log_path=log_path,
        nx2=args.nx2,
        nx3=args.nx3,
        tlim=args.tlim,
        nlim=args.nlim,
        output_dt=args.output_dt,
    )

    cycles = [int(m.group(1)) for m in re.finditer(r"cycle=(\d+)", text)]
    final_cycle = max(cycles) if cycles else -1
    driver_completed = "Driver completed." in text

    phdf = latest_phdf(workdir)
    phdf_present = phdf is not None
    metrics = None
    if phdf_present:
        metrics = parse_cons_metrics(phdf, args.b2_comp, args.b3_comp)

    b2_mean_abs = math.nan if metrics is None else metrics["b2_mean_abs"]
    b3_mean_abs = math.nan if metrics is None else metrics["b3_mean_abs"]
    b23_rel_diff = math.nan if metrics is None else metrics["b23_rel_diff"]
    finite_cons = False if metrics is None else metrics["finite_cons"]
    strict_shape = False if metrics is None else metrics["strict_shape"]
    logical_yz_single = False if metrics is None else metrics["logical_yz_single"]

    failures = []
    if returncode != 0:
        failures.append("nonzero_exit")
    if not driver_completed:
        failures.append("driver_not_completed")
    if final_cycle < args.min_final_cycle:
        failures.append("insufficient_cycles")
    if not phdf_present:
        failures.append("missing_phdf")
    if metrics is not None:
        if not finite_cons:
            failures.append("nonfinite_cons")
        if not strict_shape:
            failures.append("not_strict1d_shape")
        if not logical_yz_single:
            failures.append("not_strict1d_logical_layout")
        if b2_mean_abs < args.min_transverse_mean_abs:
            failures.append("b2_mean_abs_too_small")
        if b3_mean_abs < args.min_transverse_mean_abs:
            failures.append("b3_mean_abs_too_small")
        if b23_rel_diff > args.max_transverse_rel_diff:
            failures.append("b2_b3_mean_mismatch")

    status = "FAIL" if failures else "PASS"
    print(
        "case,returncode,final_cycle,driver_completed,phdf_present,finite_cons,strict_shape,"
        "logical_yz_single,b2_mean_abs,b3_mean_abs,b23_rel_diff,status,failures,log,phdf"
    )
    print(
        f"cpaw_strict1d_mpi_hdf5,{returncode},{final_cycle},{int(driver_completed)},"
        f"{int(phdf_present)},{int(finite_cons)},{int(strict_shape)},{int(logical_yz_single)},"
        f"{b2_mean_abs:.6e},{b3_mean_abs:.6e},{b23_rel_diff:.6e},"
        f"{status},{'|'.join(failures) if failures else 'none'},{log_path},{phdf if phdf else 'none'}"
    )

    if failures:
        raise SystemExit("cpaw strict-1d mpi+hdf5 regression failed")


if __name__ == "__main__":
    main()
