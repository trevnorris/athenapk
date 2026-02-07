#!/usr/bin/env python3
"""Regression checks for conservative EM source-channel decomposition and gain gating."""

import argparse
import math
import shutil
import subprocess
from pathlib import Path


CHANNEL_KEYS = (
    "m4d_int_src_em_laplacian_abs",
    "m4d_int_src_em_mass_abs",
    "m4d_int_src_em_current_abs",
    "m4d_int_src_em_damping_abs",
    "m4d_int_src_em_spatial_mixed_abs",
    "m4d_int_src_em_timelike_abs",
)


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
        raise RuntimeError(f"No history rows found in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for idx, name in enumerate(labels):
            cols[name].append(row[idx])
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
    cols = parse_hst(output_hst)
    return {key: require_col(cols, key)[-1] for key in CHANNEL_KEYS}


def fmt(value):
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_em_conservative_source_channels_regression",
        help="Output directory (relative to workdir if not absolute)",
    )
    parser.add_argument(
        "--harris-input",
        default=None,
        help="Harris full input deck (default: inputs/harris_4d_full.in)",
    )
    parser.add_argument(
        "--harris-tlim",
        type=float,
        default=2.0e-2,
        help="Short Harris runtime limit for this regression",
    )
    parser.add_argument(
        "--harris-nlim",
        type=int,
        default=200,
        help="Short Harris cycle limit for this regression",
    )
    parser.add_argument(
        "--damping-probe-value",
        type=float,
        default=1.0e-1,
        help="Positive em_damping used to exercise damping-source gating",
    )
    parser.add_argument("--max-cons-laplacian-abs", type=float, default=1.0e-30)
    parser.add_argument("--max-cons-spatial-mixed-abs", type=float, default=1.0e-30)
    parser.add_argument("--min-cons-mass-abs", type=float, default=1.0e-8)
    parser.add_argument("--min-cons-current-abs", type=float, default=1.0e-14)
    parser.add_argument("--min-cons-timelike-abs", type=float, default=1.0e-5)
    parser.add_argument("--min-noncons-laplacian-abs", type=float, default=1.0e-3)
    parser.add_argument("--min-noncons-spatial-mixed-abs", type=float, default=1.0e-4)
    parser.add_argument("--min-damping-on-abs", type=float, default=1.0e-6)
    parser.add_argument("--max-zero-channel-abs", type=float, default=1.0e-30)
    args = parser.parse_args()

    if args.damping_probe_value <= 0.0:
        raise SystemExit("--damping-probe-value must be > 0")

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    binary = Path(args.binary).resolve()
    input_path = (
        Path(args.harris_input).resolve()
        if args.harris_input
        else (repo_root / "inputs" / "harris_4d_full.in")
    )
    workdir = Path(args.workdir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    common_runtime = [
        f"parthenon/time/tlim={args.harris_tlim}",
        f"parthenon/time/nlim={args.harris_nlim}",
    ]

    cases = {}
    case_meta = {}

    def add_case(case_key, *, conservative, mass_gain, current_gain, damping_gain, timelike_gain,
                 damping):
        case_args = [
            f"modes4d/em_conservative_transport={'true' if conservative else 'false'}",
            f"modes4d/em_source_mass_gain={mass_gain}",
            f"modes4d/em_source_current_gain={current_gain}",
            f"modes4d/em_source_damping_gain={damping_gain}",
            f"modes4d/em_source_timelike_gain={timelike_gain}",
            f"modes4d/em_damping={damping}",
            *common_runtime,
        ]
        cases[case_key] = run_case(
            binary=binary,
            input_path=input_path,
            workdir=workdir,
            output_hst=output_dir / f"{case_key}.out.hst",
            extra_args=case_args,
        )
        case_meta[case_key] = {
            "em_conservative_transport": "true" if conservative else "false",
            "mass_gain": mass_gain,
            "current_gain": current_gain,
            "damping_gain": damping_gain,
            "timelike_gain": timelike_gain,
            "em_damping": damping,
        }

    add_case(
        "conservative_baseline",
        conservative=True,
        mass_gain=1.0,
        current_gain=1.0,
        damping_gain=1.0,
        timelike_gain=1.0,
        damping=0.0,
    )
    add_case(
        "nonconservative_baseline",
        conservative=False,
        mass_gain=1.0,
        current_gain=1.0,
        damping_gain=1.0,
        timelike_gain=1.0,
        damping=0.0,
    )
    add_case(
        "conservative_mass_gain0",
        conservative=True,
        mass_gain=0.0,
        current_gain=1.0,
        damping_gain=1.0,
        timelike_gain=1.0,
        damping=0.0,
    )
    add_case(
        "conservative_current_gain0",
        conservative=True,
        mass_gain=1.0,
        current_gain=0.0,
        damping_gain=1.0,
        timelike_gain=1.0,
        damping=0.0,
    )
    add_case(
        "conservative_timelike_gain0",
        conservative=True,
        mass_gain=1.0,
        current_gain=1.0,
        damping_gain=1.0,
        timelike_gain=0.0,
        damping=0.0,
    )
    add_case(
        "conservative_damping_on",
        conservative=True,
        mass_gain=1.0,
        current_gain=1.0,
        damping_gain=1.0,
        timelike_gain=1.0,
        damping=args.damping_probe_value,
    )
    add_case(
        "conservative_damping_gain0",
        conservative=True,
        mass_gain=1.0,
        current_gain=1.0,
        damping_gain=0.0,
        timelike_gain=1.0,
        damping=args.damping_probe_value,
    )

    case_failures = {name: [] for name in cases}
    global_failures = []

    def check_max(case_name, label, value, bound):
        if value > bound:
            case_failures[case_name].append(label)
            global_failures.append(f"{case_name}:{label}")

    def check_min(case_name, label, value, bound):
        if value < bound:
            case_failures[case_name].append(label)
            global_failures.append(f"{case_name}:{label}")

    cons = cases["conservative_baseline"]
    noncons = cases["nonconservative_baseline"]
    mass0 = cases["conservative_mass_gain0"]
    current0 = cases["conservative_current_gain0"]
    timelike0 = cases["conservative_timelike_gain0"]
    damp_on = cases["conservative_damping_on"]
    damp0 = cases["conservative_damping_gain0"]

    check_max(
        "conservative_baseline",
        "laplacian_not_suppressed",
        cons["m4d_int_src_em_laplacian_abs"],
        args.max_cons_laplacian_abs,
    )
    check_max(
        "conservative_baseline",
        "spatial_mixed_not_suppressed",
        cons["m4d_int_src_em_spatial_mixed_abs"],
        args.max_cons_spatial_mixed_abs,
    )
    check_min(
        "conservative_baseline",
        "mass_not_active",
        cons["m4d_int_src_em_mass_abs"],
        args.min_cons_mass_abs,
    )
    check_min(
        "conservative_baseline",
        "current_not_active",
        cons["m4d_int_src_em_current_abs"],
        args.min_cons_current_abs,
    )
    check_min(
        "conservative_baseline",
        "timelike_not_active",
        cons["m4d_int_src_em_timelike_abs"],
        args.min_cons_timelike_abs,
    )
    check_min(
        "nonconservative_baseline",
        "laplacian_not_active",
        noncons["m4d_int_src_em_laplacian_abs"],
        args.min_noncons_laplacian_abs,
    )
    check_min(
        "nonconservative_baseline",
        "spatial_mixed_not_active",
        noncons["m4d_int_src_em_spatial_mixed_abs"],
        args.min_noncons_spatial_mixed_abs,
    )
    check_max(
        "conservative_mass_gain0",
        "mass_not_suppressed_by_gain",
        mass0["m4d_int_src_em_mass_abs"],
        args.max_zero_channel_abs,
    )
    check_max(
        "conservative_current_gain0",
        "current_not_suppressed_by_gain",
        current0["m4d_int_src_em_current_abs"],
        args.max_zero_channel_abs,
    )
    check_max(
        "conservative_timelike_gain0",
        "timelike_not_suppressed_by_gain",
        timelike0["m4d_int_src_em_timelike_abs"],
        args.max_zero_channel_abs,
    )
    check_min(
        "conservative_damping_on",
        "damping_not_active",
        damp_on["m4d_int_src_em_damping_abs"],
        args.min_damping_on_abs,
    )
    check_max(
        "conservative_damping_gain0",
        "damping_not_suppressed_by_gain",
        damp0["m4d_int_src_em_damping_abs"],
        args.max_zero_channel_abs,
    )

    header = [
        "case",
        "em_conservative_transport",
        "mass_gain",
        "current_gain",
        "damping_gain",
        "timelike_gain",
        "em_damping",
        *CHANNEL_KEYS,
        "status",
        "failures",
    ]
    print(",".join(header))
    for case_name, values in cases.items():
        failures = case_failures[case_name]
        status = "FAIL" if failures else "PASS"
        print(
            ",".join(
                [
                    case_name,
                    case_meta[case_name]["em_conservative_transport"],
                    fmt(case_meta[case_name]["mass_gain"]),
                    fmt(case_meta[case_name]["current_gain"]),
                    fmt(case_meta[case_name]["damping_gain"]),
                    fmt(case_meta[case_name]["timelike_gain"]),
                    fmt(case_meta[case_name]["em_damping"]),
                    *(fmt(values[key]) for key in CHANNEL_KEYS),
                    status,
                    "|".join(failures) if failures else "none",
                ]
            )
        )

    if global_failures:
        raise SystemExit(
            "conservative source-channel regression failed: "
            + ", ".join(global_failures)
        )


if __name__ == "__main__":
    main()
