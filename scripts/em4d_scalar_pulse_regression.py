#!/usr/bin/env python3
"""Run scalar-channel pulse regression checks for modes4d."""

import argparse
import math
import shutil
import subprocess
import sys
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


def run_case(binary, input_path, workdir, output_hst):
    for filename in ("parthenon.out0.hst", "parthenon.out1.hst"):
        try:
            (workdir / filename).unlink()
        except FileNotFoundError:
            pass

    cmd = [str(binary), "-i", str(input_path)]
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


def require_col(cols, key):
    values = cols.get(key)
    if values is None:
        raise RuntimeError(f"Missing required history column: {key}")
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument(
        "--input",
        default=None,
        help="Scalar-channel pulse input deck",
    )
    parser.add_argument(
        "--workdir",
        default=".",
        help="Run directory",
    )
    parser.add_argument(
        "--output-hst",
        default="em4d_scalar_pulse.out.hst",
        help="Path for copied history output (relative to workdir if not absolute)",
    )
    parser.add_argument(
        "--min-max-mixed-ew2",
        type=float,
        default=1.0e-6,
        help="Minimum max(m4d_mixed_ew2) for PASS",
    )
    parser.add_argument(
        "--min-max-mixed-c2",
        type=float,
        default=1.0e-6,
        help="Minimum max(m4d_mixed_c2) for PASS",
    )
    parser.add_argument(
        "--min-max-em-a2-mode1",
        type=float,
        default=1.0e-8,
        help="Minimum max(m4d_em_a2_mode_1) for PASS",
    )
    parser.add_argument(
        "--min-pulse-xc-shift",
        type=float,
        default=1.0e-3,
        help="Minimum |m4d_pulse_xc(final)-m4d_pulse_xc(initial)| for PASS",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]

    binary = Path(args.binary).resolve()
    if args.input is None:
        input_path = repo_root / "inputs" / "em4d_scalar_pulse.in"
    else:
        input_path = Path(args.input).resolve()
    workdir = Path(args.workdir).resolve()
    output_hst = Path(args.output_hst)
    if not output_hst.is_absolute():
        output_hst = workdir / output_hst
    output_hst.parent.mkdir(parents=True, exist_ok=True)

    cols = run_case(binary, input_path, workdir, output_hst)
    mixed_ew2 = require_col(cols, "m4d_mixed_ew2")
    mixed_c2 = require_col(cols, "m4d_mixed_c2")
    em_a2_mode1 = require_col(cols, "m4d_em_a2_mode_1")
    pulse_xc = require_col(cols, "m4d_pulse_xc")
    jw_ew = cols.get("m4d_int_jw_ew", [math.nan])
    s_leak_abs = cols.get("m4d_int_s_leak_abs", [math.nan])

    max_mixed_ew2 = max(mixed_ew2)
    max_mixed_c2 = max(mixed_c2)
    max_em_a2_mode1 = max(em_a2_mode1)
    pulse_xc_shift = abs(pulse_xc[-1] - pulse_xc[0])
    final_jw_ew = jw_ew[-1]
    final_s_leak_abs = s_leak_abs[-1]

    failures = []
    if max_mixed_ew2 < args.min_max_mixed_ew2:
        failures.append("mixed_ew2")
    if max_mixed_c2 < args.min_max_mixed_c2:
        failures.append("mixed_c2")
    if max_em_a2_mode1 < args.min_max_em_a2_mode1:
        failures.append("em_a2_mode_1")
    if pulse_xc_shift < args.min_pulse_xc_shift:
        failures.append("pulse_xc_shift")

    print(
        "case,max_mixed_ew2,max_mixed_c2,max_em_a2_mode_1,pulse_xc_shift,final_jw_ew,final_s_leak_abs,status,failures"
    )
    status = "FAIL" if failures else "PASS"
    fail_text = "|".join(failures) if failures else "none"
    print(
        f"scalar_pulse,{max_mixed_ew2:.6e},{max_mixed_c2:.6e},{max_em_a2_mode1:.6e},"
        f"{pulse_xc_shift:.6e},{final_jw_ew:.6e},{final_s_leak_abs:.6e},{status},{fail_text}"
    )

    if failures:
        raise SystemExit(f"scalar pulse checks failed: {', '.join(failures)}")


if __name__ == "__main__":
    main()
