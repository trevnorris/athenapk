#!/usr/bin/env python3
"""Targeted regression checks for time-like EM mixed couplings."""

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--input", default=None, help="EM pulse input deck")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_em_timelike_regression",
        help="Output directory (relative to workdir if not absolute)",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=1.0e9,
        help="Pulse sigma override to suppress spatial-gradient couplings",
    )
    parser.add_argument(
        "--tlim",
        type=float,
        default=2.0e-2,
        help="Short runtime for the targeted coupling checks",
    )
    parser.add_argument(
        "--nlim",
        type=int,
        default=400,
        help="Iteration cap for targeted checks",
    )
    parser.add_argument(
        "--min-aw0-to-a0mode1-pi2-mode1",
        type=float,
        default=1.0e-10,
        help="Minimum max(m4d_em_pi2_mode_1) for the aw_mode0 -> a0_mode1 coupling check",
    )
    parser.add_argument(
        "--min-a0mode1-to-aw0-piw0",
        type=float,
        default=1.0e-5,
        help="Minimum max(|m4d_piw_mode_0|) for the a0_mode1 -> aw_mode0 coupling check",
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
        output_dir / "aw0_to_a0mode1.out.hst",
        common + ["problem/em4d_pulse/component=aw_mode0"],
    )
    aw0_pi2_mode1 = max(require_col(cols_aw0, "m4d_em_pi2_mode_1"))
    aw0_piw0 = max(abs(v) for v in require_col(cols_aw0, "m4d_piw_mode_0"))
    aw0_mixed_c2 = max(require_col(cols_aw0, "m4d_mixed_c2"))
    aw0_failures = []
    if aw0_pi2_mode1 < args.min_aw0_to_a0mode1_pi2_mode1:
        aw0_failures.append("aw0_to_a0mode1_em_pi2_mode_1")
    cases.append(
        {
            "case": "aw0_to_a0mode1_timelike",
            "max_em_pi2_mode_1": aw0_pi2_mode1,
            "max_abs_piw_mode_0": aw0_piw0,
            "max_mixed_c2": aw0_mixed_c2,
            "status": "FAIL" if aw0_failures else "PASS",
            "failures": "|".join(aw0_failures) if aw0_failures else "none",
        }
    )

    cols_a0m1 = run_case(
        binary,
        input_path,
        workdir,
        output_dir / "a0mode1_to_aw0.out.hst",
        common + ["problem/em4d_pulse/component=a0_mode1"],
    )
    a0m1_a2_mode0 = max(require_col(cols_a0m1, "m4d_em_a2_mode_0"))
    a0m1_piw0 = max(abs(v) for v in require_col(cols_a0m1, "m4d_piw_mode_0"))
    a0m1_mixed_ew2 = max(require_col(cols_a0m1, "m4d_mixed_ew2"))
    a0m1_failures = []
    if a0m1_piw0 < args.min_a0mode1_to_aw0_piw0:
        a0m1_failures.append("a0mode1_to_aw0_abs_piw_mode_0")
    cases.append(
        {
            "case": "a0mode1_to_aw0_timelike",
            "max_em_a2_mode_0": a0m1_a2_mode0,
            "max_abs_piw_mode_0": a0m1_piw0,
            "max_mixed_ew2": a0m1_mixed_ew2,
            "status": "FAIL" if a0m1_failures else "PASS",
            "failures": "|".join(a0m1_failures) if a0m1_failures else "none",
        }
    )

    print(
        "case,max_em_pi2_mode_1,max_em_a2_mode_0,max_abs_piw_mode_0,max_mixed_c2,max_mixed_ew2,status,failures"
    )
    for case in cases:
        print(
            f"{case['case']},"
            f"{case.get('max_em_pi2_mode_1', math.nan):.6e},"
            f"{case.get('max_em_a2_mode_0', math.nan):.6e},"
            f"{case.get('max_abs_piw_mode_0', math.nan):.6e},"
            f"{case.get('max_mixed_c2', math.nan):.6e},"
            f"{case.get('max_mixed_ew2', math.nan):.6e},"
            f"{case['status']},{case['failures']}"
        )

    failures = [c for c in cases if c["status"] != "PASS"]
    if failures:
        failed = ", ".join(c["case"] for c in failures)
        raise SystemExit(f"time-like EM coupling regression failed: {failed}")


if __name__ == "__main__":
    main()
