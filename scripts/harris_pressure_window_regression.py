#!/usr/bin/env python3
"""Run tuned Harris regression over a pressure-transport gain window."""

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


def parse_gains(raw):
    gains = []
    for token in raw.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        value = float(stripped)
        if value < 0.0:
            raise ValueError(f"pressure gain must be >= 0. Got {value}")
        gains.append(stripped)
    if not gains:
        raise ValueError("At least one pressure gain must be provided.")
    return gains


def make_gain_tag(raw):
    return raw.replace("+", "p").replace("-", "m")


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Directory to run simulations in")
    parser.add_argument(
        "--output-dir",
        default="harris_pressure_window_regression_outputs",
        help="Directory for pressure-window outputs",
    )
    parser.add_argument(
        "--gains",
        default="1e-8,1e-7,1e-6,1e-3",
        help="Comma-separated pressure gains to test in full-case overrides",
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

    gains = parse_gains(args.gains)
    required_finite = [
        "final_psi0_span",
        "final_s_leak_abs",
        "final_jw_ew",
        "final_mixed_ew2",
        "final_mixed_c2",
        "final_em_leak_w",
        "final_helicity_sub",
        "final_edotb_sub",
        "closure_max_norm",
        "closure_max_abs_rate",
        "closure_local_mode0_max_abs_rate",
    ]
    required_status = [
        "closure_status",
        "closure_local_mode0_status",
        "correlation_status",
        "activity_status",
        "em_transport_closure_status",
        "transport_closure_status",
    ]

    failures = []
    summary_rows = []
    for gain_raw in gains:
        gain_tag = make_gain_tag(gain_raw)
        case_dir = output_dir / f"gain_{gain_tag}"
        case_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable,
            str(tuned_script),
            "--binary",
            str(binary),
            "--workdir",
            str(case_dir),
            "--output-dir",
            str(case_dir / "outputs"),
            "--scan-full-arg",
            f"modes4d/plasma_pressure_transport_gain={gain_raw}",
        ]
        for arg in args.scan_athena_arg:
            cmd.extend(["--scan-athena-arg", arg])
        for arg in args.scan_controlled_arg:
            cmd.extend(["--scan-controlled-arg", arg])
        for arg in args.scan_full_arg:
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
            failures.append(
                f"gain={gain_raw}: tuned regression failed "
                f"(rc={proc.returncode}); see {case_dir}/scan_stderr.log"
            )
            continue

        try:
            parsed = parse_scan_csv(proc.stdout)
        except Exception as err:
            failures.append(f"gain={gain_raw}: could not parse scan CSV ({err})")
            continue

        full = parsed["full"]
        finite_failures = []
        for key in required_finite:
            value = parse_float(full, key)
            if not math.isfinite(value):
                finite_failures.append(f"{key}=nan/inf")
        status_failures = []
        for key in required_status:
            value = full.get(key, "N/A")
            if value != "PASS":
                status_failures.append(f"{key}={value}")
        if finite_failures or status_failures:
            joined = "|".join(finite_failures + status_failures)
            failures.append(f"gain={gain_raw}: {joined}")
            continue

        summary_rows.append(
            {
                "gain": gain_raw,
                "closure_status": full.get("closure_status", "N/A"),
                "transport_closure_status": full.get("transport_closure_status", "N/A"),
                "correlation_status": full.get("correlation_status", "N/A"),
                "activity_status": full.get("activity_status", "N/A"),
                "final_mixed_ew2": parse_float(full, "final_mixed_ew2"),
                "final_s_leak_abs": parse_float(full, "final_s_leak_abs"),
                "final_jw_ew": parse_float(full, "final_jw_ew"),
                "full_corr_psi0_mixed_ew2": parse_float(full, "corr_psi0_mixed_ew2"),
                "full_corr_psi0_s_leak_abs": parse_float(full, "corr_psi0_s_leak_abs"),
            }
        )

    header = [
        "gain",
        "closure_status",
        "transport_closure_status",
        "correlation_status",
        "activity_status",
        "final_mixed_ew2",
        "final_s_leak_abs",
        "final_jw_ew",
        "full_corr_psi0_mixed_ew2",
        "full_corr_psi0_s_leak_abs",
    ]
    print(",".join(header))
    for row in summary_rows:
        print(",".join(fmt(row[key]) for key in header))

    if failures:
        raise SystemExit("pressure-window regression failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
