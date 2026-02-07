#!/usr/bin/env python3
"""Conservative-path diagnostics checks for time-like EM mixed couplings."""

import argparse
import math
import shutil
import subprocess
from pathlib import Path


def parse_hst(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    header = None
    for line in lines:
        if line.startswith("# [1]="):
            header = line
            break
    if header is None:
        raise RuntimeError(f"Could not find history header in {path}")

    labels = []
    for token in header.replace("#", "").split():
        if token.startswith("[") and "]=" in token:
            labels.append(token.split("]=", 1)[1])

    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])

    if not rows:
        raise RuntimeError(f"No history data rows in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


def require_col(cols, key):
    values = cols.get(key)
    if values is None:
        raise RuntimeError(f"Missing required history column: {key}")
    for value in values:
        if not math.isfinite(value):
            raise RuntimeError(f"Non-finite value in {key}")
    return values


def run_case(binary, input_path, workdir, output_hst, extra_args):
    for filename in ("parthenon.out0.hst", "parthenon.out1.hst"):
        try:
            (workdir / filename).unlink()
        except FileNotFoundError:
            pass

    cmd = [str(binary), "-i", str(input_path)]
    cmd.extend(extra_args)
    subprocess.run(cmd, cwd=workdir, check=True)

    produced = None
    for candidate in ("parthenon.out1.hst", "parthenon.out0.hst"):
        path = workdir / candidate
        if path.exists():
            produced = path
            break
    if produced is None:
        raise RuntimeError("Expected history output not found (parthenon.out1.hst/out0.hst)")

    shutil.copyfile(produced, output_hst)
    return parse_hst(output_hst)


def fmt(value):
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--input", default=None, help="EM pulse input deck")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_em_conservative_timelike_diag_regression",
        help="Output directory (relative to workdir if not absolute)",
    )
    parser.add_argument("--sigma", type=float, default=1.0e9)
    parser.add_argument("--tlim", type=float, default=2.0e-2)
    parser.add_argument("--nlim", type=int, default=400)
    parser.add_argument(
        "--min-aw0-a0-from-piw-abs-delta",
        type=float,
        default=0.0,
        help="Optional weak floor for aw_mode0 case; defaults to disabled (0.0)",
    )
    parser.add_argument(
        "--min-a0mode1-a0-from-piw-abs-delta",
        type=float,
        default=1.0e-25,
        help="Weak floor for back-coupled a0<-piw diagnostic in a0_mode1 case",
    )
    parser.add_argument(
        "--min-a0mode1-aw-from-pi0-abs-delta",
        type=float,
        default=1.0e-50,
        help="Weak floor for targeted aw<-pi0 diagnostic in a0_mode1 case",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    binary = Path(args.binary).resolve()
    input_path = Path(args.input).resolve() if args.input else repo_root / "inputs" / "em4d_pulse.in"
    workdir = Path(args.workdir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    common = [
        "modes4d/em_conservative_transport=true",
        f"problem/em4d_pulse/sigma={args.sigma}",
        "problem/em4d_pulse/pi_pulse_scale=1.0",
        f"parthenon/time/tlim={args.tlim}",
        f"parthenon/time/nlim={args.nlim}",
    ]

    cases = []

    cols_aw0 = run_case(
        binary,
        input_path,
        workdir,
        output_dir / "aw0_conservative.out.hst",
        common + ["problem/em4d_pulse/component=aw_mode0"],
    )
    aw0_a0_from_piw = require_col(cols_aw0, "m4d_int_src_timelike_a0_from_piw")
    aw0_a0_from_piw_abs = require_col(cols_aw0, "m4d_int_src_timelike_a0_from_piw_abs")
    aw0_aw_from_pi0 = require_col(cols_aw0, "m4d_int_src_timelike_aw_from_pi0")
    aw0_aw_from_pi0_abs = require_col(cols_aw0, "m4d_int_src_timelike_aw_from_pi0_abs")
    aw0_a0_from_piw_abs_delta = aw0_a0_from_piw_abs[-1] - aw0_a0_from_piw_abs[0]
    aw0_failures = []
    if (
        args.min_aw0_a0_from_piw_abs_delta > 0.0
        and aw0_a0_from_piw_abs_delta < args.min_aw0_a0_from_piw_abs_delta
    ):
        aw0_failures.append("aw0_timelike_a0_from_piw_abs_delta")
    cases.append(
        {
            "case": "aw0_conservative_timelike_diag",
            "final_int_src_timelike_a0_from_piw": aw0_a0_from_piw[-1],
            "final_int_src_timelike_a0_from_piw_abs": aw0_a0_from_piw_abs[-1],
            "delta_int_src_timelike_a0_from_piw_abs": aw0_a0_from_piw_abs_delta,
            "final_int_src_timelike_aw_from_pi0": aw0_aw_from_pi0[-1],
            "final_int_src_timelike_aw_from_pi0_abs": aw0_aw_from_pi0_abs[-1],
            "status": "FAIL" if aw0_failures else "PASS",
            "failures": "|".join(aw0_failures) if aw0_failures else "none",
        }
    )

    cols_a0m1 = run_case(
        binary,
        input_path,
        workdir,
        output_dir / "a0mode1_conservative.out.hst",
        common + ["problem/em4d_pulse/component=a0_mode1"],
    )
    a0m1_a0_from_piw = require_col(cols_a0m1, "m4d_int_src_timelike_a0_from_piw")
    a0m1_a0_from_piw_abs = require_col(cols_a0m1, "m4d_int_src_timelike_a0_from_piw_abs")
    a0m1_aw_from_pi0 = require_col(cols_a0m1, "m4d_int_src_timelike_aw_from_pi0")
    a0m1_aw_from_pi0_abs = require_col(cols_a0m1, "m4d_int_src_timelike_aw_from_pi0_abs")
    a0m1_a0_from_piw_abs_delta = a0m1_a0_from_piw_abs[-1] - a0m1_a0_from_piw_abs[0]
    a0m1_aw_from_pi0_abs_delta = a0m1_aw_from_pi0_abs[-1] - a0m1_aw_from_pi0_abs[0]
    a0m1_failures = []
    if (
        args.min_a0mode1_a0_from_piw_abs_delta > 0.0
        and a0m1_a0_from_piw_abs_delta < args.min_a0mode1_a0_from_piw_abs_delta
    ):
        a0m1_failures.append("a0mode1_timelike_a0_from_piw_abs_delta")
    if (
        args.min_a0mode1_aw_from_pi0_abs_delta > 0.0
        and a0m1_aw_from_pi0_abs_delta < args.min_a0mode1_aw_from_pi0_abs_delta
    ):
        a0m1_failures.append("a0mode1_timelike_aw_from_pi0_abs_delta")
    cases.append(
        {
            "case": "a0mode1_conservative_timelike_diag",
            "final_int_src_timelike_a0_from_piw": a0m1_a0_from_piw[-1],
            "final_int_src_timelike_a0_from_piw_abs": a0m1_a0_from_piw_abs[-1],
            "delta_int_src_timelike_a0_from_piw_abs": a0m1_a0_from_piw_abs[-1]
            - a0m1_a0_from_piw_abs[0],
            "final_int_src_timelike_aw_from_pi0": a0m1_aw_from_pi0[-1],
            "final_int_src_timelike_aw_from_pi0_abs": a0m1_aw_from_pi0_abs[-1],
            "delta_int_src_timelike_aw_from_pi0_abs": a0m1_aw_from_pi0_abs_delta,
            "status": "FAIL" if a0m1_failures else "PASS",
            "failures": "|".join(a0m1_failures) if a0m1_failures else "none",
        }
    )

    print(
        "case,final_int_src_timelike_a0_from_piw,final_int_src_timelike_a0_from_piw_abs,delta_int_src_timelike_a0_from_piw_abs,final_int_src_timelike_aw_from_pi0,final_int_src_timelike_aw_from_pi0_abs,delta_int_src_timelike_aw_from_pi0_abs,status,failures"
    )
    for case in cases:
        print(
            ",".join(
                [
                    case["case"],
                    fmt(case["final_int_src_timelike_a0_from_piw"]),
                    fmt(case["final_int_src_timelike_a0_from_piw_abs"]),
                    fmt(case["delta_int_src_timelike_a0_from_piw_abs"]),
                    fmt(case["final_int_src_timelike_aw_from_pi0"]),
                    fmt(case["final_int_src_timelike_aw_from_pi0_abs"]),
                    fmt(case.get("delta_int_src_timelike_aw_from_pi0_abs", math.nan)),
                    case["status"],
                    case["failures"],
                ]
            )
        )

    failures = [c for c in cases if c["status"] != "PASS"]
    if failures:
        failed = ", ".join(c["case"] for c in failures)
        raise SystemExit(f"conservative timelike diagnostics regression failed: {failed}")


if __name__ == "__main__":
    main()
