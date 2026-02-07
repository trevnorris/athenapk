#!/usr/bin/env python3
"""Run tuned Harris regression with optional w-flux source enabled and gated."""

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def parse_scan_csv(stdout_text):
    lines = [line.strip() for line in stdout_text.splitlines() if line.strip()]
    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("case,"):
            header_idx = idx
            break
    if header_idx is None:
        raise RuntimeError("Could not find CSV header from harris_scan_matrix.py output.")

    rows = []
    for line in lines[header_idx + 1 :]:
        if line.startswith("controlled,") or line.startswith("full,"):
            rows.append(line)
        if len(rows) == 2:
            break
    if len(rows) < 2:
        raise RuntimeError("Could not find controlled/full CSV rows from scan output.")

    reader = csv.DictReader([lines[header_idx], *rows])
    parsed = {row["case"]: row for row in reader}
    if "controlled" not in parsed or "full" not in parsed:
        raise RuntimeError("Scan CSV rows missing controlled/full cases.")
    return parsed


def parse_float(row, key):
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return math.nan


def fmt(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def run_tuned_case(
    *,
    script,
    binary,
    workdir,
    case_dir,
    gain,
    scan_athena_args,
    scan_controlled_args,
    scan_full_args,
):
    case_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(script),
        "--binary",
        str(binary),
        "--workdir",
        str(case_dir),
        "--output-dir",
        str(case_dir / "outputs"),
        "--scan-full-arg",
        f"modes4d/plasma_w_flux_source_gain={gain}",
    ]
    for arg in scan_athena_args:
        cmd.extend(["--scan-athena-arg", arg])
    for arg in scan_controlled_args:
        cmd.extend(["--scan-controlled-arg", arg])
    for arg in scan_full_args:
        cmd.extend(["--scan-full-arg", arg])

    proc = subprocess.run(
        cmd,
        cwd=workdir,
        text=True,
        capture_output=True,
        check=False,
    )
    (case_dir / "scan_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (case_dir / "scan_stderr.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            f"tuned regression failed for gain={gain} (rc={proc.returncode}); "
            f"see {case_dir}/scan_stderr.log"
        )
    parsed = parse_scan_csv(proc.stdout)
    return parsed["full"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Directory to run simulations in")
    parser.add_argument(
        "--output-dir",
        default="harris_wflux_regression_outputs",
        help="Directory for w-flux regression outputs",
    )
    parser.add_argument(
        "--baseline-gain",
        type=float,
        default=0.0,
        help="Baseline full-case plasma_w_flux_source_gain (default: 0.0)",
    )
    parser.add_argument(
        "--enabled-gain",
        type=float,
        default=1.0e-1,
        help="Enabled full-case plasma_w_flux_source_gain (default: 0.1)",
    )
    parser.add_argument(
        "--min-enabled-jw-mode1-abs",
        type=float,
        default=1.0e-18,
        help="Minimum |final_jw_mode_1| with enabled-gain",
    )
    parser.add_argument(
        "--min-enabled-jw-mode1-delta-abs",
        type=float,
        default=1.0e-18,
        help="Minimum |delta final_jw_mode_1| between enabled and baseline",
    )
    parser.add_argument(
        "--min-enabled-jw-mode1-growth-factor",
        type=float,
        default=10.0,
        help="Minimum growth factor in |final_jw_mode_1| from baseline to enabled (if baseline != 0)",
    )
    parser.add_argument(
        "--scan-athena-arg",
        action="append",
        default=[],
        help="Additional athena override forwarded to both controlled/full runs",
    )
    parser.add_argument(
        "--scan-controlled-arg",
        action="append",
        default=[],
        help="Additional athena override forwarded only to controlled run",
    )
    parser.add_argument(
        "--scan-full-arg",
        action="append",
        default=[],
        help="Additional athena override forwarded only to full run",
    )
    args = parser.parse_args()

    if args.enabled_gain <= args.baseline_gain:
        raise SystemExit("enabled-gain must be greater than baseline-gain")
    if args.min_enabled_jw_mode1_growth_factor <= 0.0:
        raise SystemExit("min-enabled-jw-mode1-growth-factor must be > 0")

    repo_root = Path(__file__).resolve().parents[1]
    tuned_script = repo_root / "scripts" / "harris_regression_tuned.py"

    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    output_dir_arg = Path(args.output_dir)
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (workdir / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        full_baseline = run_tuned_case(
            script=tuned_script,
            binary=binary,
            workdir=workdir,
            case_dir=output_dir / "baseline",
            gain=args.baseline_gain,
            scan_athena_args=args.scan_athena_arg,
            scan_controlled_args=args.scan_controlled_arg,
            scan_full_args=args.scan_full_arg,
        )
        full_enabled = run_tuned_case(
            script=tuned_script,
            binary=binary,
            workdir=workdir,
            case_dir=output_dir / "enabled",
            gain=args.enabled_gain,
            scan_athena_args=args.scan_athena_arg,
            scan_controlled_args=args.scan_controlled_arg,
            scan_full_args=args.scan_full_arg,
        )
    except Exception as err:
        raise SystemExit(str(err))

    baseline_abs = abs(parse_float(full_baseline, "final_jw_mode_1"))
    enabled_abs = abs(parse_float(full_enabled, "final_jw_mode_1"))
    delta_abs = abs(enabled_abs - baseline_abs)
    growth_factor = (
        math.inf
        if baseline_abs == 0.0
        else (enabled_abs / baseline_abs)
    )

    failures = []
    if not math.isfinite(enabled_abs):
        failures.append("enabled final_jw_mode_1 is nan/inf")
    if enabled_abs < args.min_enabled_jw_mode1_abs:
        failures.append(
            f"enabled |final_jw_mode_1|={enabled_abs:.6e} < {args.min_enabled_jw_mode1_abs:.6e}"
        )
    if delta_abs < args.min_enabled_jw_mode1_delta_abs:
        failures.append(
            f"|delta final_jw_mode_1|={delta_abs:.6e} < {args.min_enabled_jw_mode1_delta_abs:.6e}"
        )
    if math.isfinite(growth_factor) and growth_factor < args.min_enabled_jw_mode1_growth_factor:
        failures.append(
            f"growth factor={growth_factor:.6e} < {args.min_enabled_jw_mode1_growth_factor:.6e}"
        )

    rows = [
        {
            "case": "baseline",
            "gain": args.baseline_gain,
            "final_jw_mode_1_abs": baseline_abs,
            "final_jw_mode_l2_1": parse_float(full_baseline, "final_jw_mode_l2_1"),
            "closure_status": full_baseline.get("closure_status", "N/A"),
            "transport_closure_status": full_baseline.get("transport_closure_status", "N/A"),
            "correlation_status": full_baseline.get("correlation_status", "N/A"),
            "activity_status": full_baseline.get("activity_status", "N/A"),
        },
        {
            "case": "enabled",
            "gain": args.enabled_gain,
            "final_jw_mode_1_abs": enabled_abs,
            "final_jw_mode_l2_1": parse_float(full_enabled, "final_jw_mode_l2_1"),
            "closure_status": full_enabled.get("closure_status", "N/A"),
            "transport_closure_status": full_enabled.get("transport_closure_status", "N/A"),
            "correlation_status": full_enabled.get("correlation_status", "N/A"),
            "activity_status": full_enabled.get("activity_status", "N/A"),
        },
    ]
    header = [
        "case",
        "gain",
        "final_jw_mode_1_abs",
        "final_jw_mode_l2_1",
        "closure_status",
        "transport_closure_status",
        "correlation_status",
        "activity_status",
    ]
    print(",".join(header))
    for row in rows:
        print(",".join(fmt(row[key]) for key in header))
    print(
        "delta,"
        + ",".join(
            [
                f"enabled_minus_baseline_final_jw_mode_1_abs={delta_abs:.6e}",
                f"growth_factor_final_jw_mode_1_abs={growth_factor:.6e}",
            ]
        )
    )

    if failures:
        raise SystemExit("w-flux regression failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
