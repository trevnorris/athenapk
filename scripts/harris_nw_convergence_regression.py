#!/usr/bin/env python3
"""Run Harris scans across N_w and gate high-mode stabilization."""

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


ABS_METRICS = {
    "final_jw_ew",
    "final_em_leak_w",
    "final_helicity_sub",
    "final_edotb_sub",
}

CONVERGENCE_METRICS = [
    "final_psi0_span",
    "final_s_leak_abs",
    "final_jw_ew",
    "final_mixed_ew2",
    "final_mixed_c2",
    "final_em_leak_w",
    "final_helicity_sub",
    "final_edotb_sub",
]


def parse_modes(raw_value):
    values = []
    for token in raw_value.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        value = int(stripped)
        if value < 2:
            raise ValueError(f"N_w values must be >= 2, got {value}")
        values.append(value)
    if len(values) < 2:
        raise ValueError("Provide at least two N_w values for convergence checks.")
    return sorted(set(values))


def resolve_input_path(raw_value, repo_root, workdir):
    path = Path(raw_value)
    if path.is_absolute():
        return path.resolve()
    for base in (Path.cwd(), workdir, repo_root):
        candidate = (base / path).resolve()
        if candidate.exists():
            return candidate
    return (repo_root / path).resolve()


def parse_scan_csv(stdout_text):
    lines = [line.strip() for line in stdout_text.splitlines() if line.strip()]
    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("case,"):
            header_idx = idx
            break
    if header_idx is None:
        raise RuntimeError("Could not find scan CSV header in regression output.")

    rows = []
    for line in lines[header_idx + 1 :]:
        if line.startswith("controlled,") or line.startswith("full,"):
            rows.append(line)
        if len(rows) == 2:
            break
    if len(rows) < 2:
        raise RuntimeError("Could not find controlled/full scan rows in regression output.")

    reader = csv.DictReader([lines[header_idx], *rows])
    parsed = {row["case"]: row for row in reader}
    if "controlled" not in parsed or "full" not in parsed:
        raise RuntimeError("Regression output missing controlled/full rows.")
    return parsed


def parse_float(row, key):
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return math.nan


def metric_value(row, metric):
    value = parse_float(row, metric)
    if metric in ABS_METRICS and not math.isnan(value):
        return abs(value)
    return value


def rel_delta(current, previous, min_scale):
    if math.isnan(current) or math.isnan(previous):
        return math.nan
    denom = max(abs(current), abs(previous), min_scale)
    return abs(current - previous) / denom


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument(
        "--workdir",
        default=".",
        help="Run directory used by harris_regression_tuned.py",
    )
    parser.add_argument(
        "--output-dir",
        default="harris_nw_convergence_outputs",
        help="Output directory (relative to workdir)",
    )
    parser.add_argument(
        "--controlled-input",
        default="inputs/harris_4d_controlled.in",
        help="Controlled-limit input deck",
    )
    parser.add_argument(
        "--full-input",
        default="inputs/harris_4d_full.in",
        help="Full-channel input deck",
    )
    parser.add_argument(
        "--n-modes",
        default="2,3,4",
        help="Comma-separated N_w values to check for high-mode stabilization",
    )
    parser.add_argument(
        "--nq-offset",
        type=int,
        default=2,
        help="Use n_quadrature = N_w + nq-offset for each scan run",
    )
    parser.add_argument(
        "--convergence-rel-tol",
        type=float,
        default=7.0e-1,
        help="Maximum allowed relative delta between the two highest N_w runs",
    )
    parser.add_argument(
        "--convergence-min-scale",
        type=float,
        default=1.0e-10,
        help="Minimum denominator scale for relative-delta normalization",
    )
    parser.add_argument(
        "--scan-arg",
        action="append",
        default=[],
        help="Extra athena override forwarded to all scan runs",
    )
    parser.add_argument(
        "--scan-controlled-arg",
        action="append",
        default=[],
        help="Extra athena override forwarded only to controlled runs",
    )
    parser.add_argument(
        "--scan-full-arg",
        action="append",
        default=[],
        help="Extra athena override forwarded only to full runs",
    )
    parser.add_argument(
        "--check-transport-closure",
        action="store_true",
        help="Require transport closure PASS in each N_w run",
    )
    parser.add_argument(
        "--highest-min-jw-ew-abs",
        type=float,
        default=1.0e-15,
        help="Minimum |final_jw_ew| required at highest N_w",
    )
    parser.add_argument(
        "--highest-min-s-leak-abs",
        type=float,
        default=1.0e-14,
        help="Minimum final_s_leak_abs required at highest N_w",
    )
    parser.add_argument(
        "--highest-min-mixed-ew2",
        type=float,
        default=1.0e-3,
        help="Minimum final_mixed_ew2 required at highest N_w",
    )
    parser.add_argument(
        "--highest-min-abs-em-leak-w",
        type=float,
        default=1.0e-4,
        help="Minimum |final_em_leak_w| required at highest N_w",
    )
    parser.add_argument(
        "--highest-min-abs-helicity-sub",
        type=float,
        default=1.0e-12,
        help="Minimum |final_helicity_sub| required at highest N_w",
    )
    parser.add_argument(
        "--highest-min-abs-edotb-sub",
        type=float,
        default=1.0e-11,
        help="Minimum |final_edotb_sub| required at highest N_w",
    )
    parser.add_argument(
        "--fail-on-check",
        action="store_true",
        help="Exit nonzero if N_w closure or convergence checks fail",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    scan_script = repo_root / "scripts" / "harris_scan_matrix.py"

    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    controlled_input = resolve_input_path(args.controlled_input, repo_root, workdir)
    full_input = resolve_input_path(args.full_input, repo_root, workdir)
    output_dir_arg = Path(args.output_dir)
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (workdir / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    n_modes_values = parse_modes(args.n_modes)
    if args.nq_offset < 0:
        raise ValueError("--nq-offset must be >= 0")
    if args.convergence_rel_tol < 0.0:
        raise ValueError("--convergence-rel-tol must be >= 0")
    if args.convergence_min_scale <= 0.0:
        raise ValueError("--convergence-min-scale must be > 0")

    rows = []
    for n_modes in n_modes_values:
        n_quadrature = n_modes + args.nq_offset
        case_dir = output_dir / f"n_modes_{n_modes}"
        case_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            sys.executable,
            str(scan_script),
            "--binary",
            str(binary),
            "--workdir",
            str(workdir),
            "--controlled-input",
            str(controlled_input),
            "--full-input",
            str(full_input),
            "--output-dir",
            str(case_dir / "scan_outputs"),
            "--controlled-arg",
            f"modes4d/n_modes=1",
            "--controlled-arg",
            "modes4d/n_quadrature=3",
            "--full-arg",
            f"modes4d/n_modes={n_modes}",
            "--full-arg",
            f"modes4d/n_quadrature={n_quadrature}",
            "--full-arg",
            "modes4d/plasma_pseudospectral_transport=true",
        ]
        if args.check_transport_closure:
            cmd.append("--check-transport-closure")
        for scan_arg in args.scan_arg:
            cmd.extend(["--athena-arg", str(scan_arg)])
        for scan_arg in args.scan_controlled_arg:
            cmd.extend(["--controlled-arg", str(scan_arg)])
        for scan_arg in args.scan_full_arg:
            cmd.extend(["--full-arg", str(scan_arg)])

        proc = subprocess.run(cmd, cwd=workdir, text=True, capture_output=True)
        (case_dir / "scan_stdout.log").write_text(proc.stdout, encoding="utf-8")
        (case_dir / "scan_stderr.log").write_text(proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise RuntimeError(
                f"Harris scan failed for N_w={n_modes}, "
                f"returncode={proc.returncode}. See {case_dir}/scan_stderr.log"
            )

        parsed = parse_scan_csv(proc.stdout)
        full_row = parsed["full"]
        rows.append(
            {
                "n_modes": n_modes,
                "n_quadrature": n_quadrature,
                "full_closure_status": full_row.get("closure_status", "N/A"),
                "full_local_closure_status": full_row.get(
                    "closure_local_mode0_status", "N/A"
                ),
                "full_transport_status": full_row.get("transport_closure_status", "N/A"),
                **{metric: metric_value(full_row, metric) for metric in CONVERGENCE_METRICS},
            }
        )

    convergence_rows = []
    failures = []
    for metric in CONVERGENCE_METRICS:
        deltas = []
        for i in range(1, len(rows)):
            current = rows[i][metric]
            previous = rows[i - 1][metric]
            delta = rel_delta(current, previous, args.convergence_min_scale)
            if not math.isnan(delta):
                deltas.append(delta)
        final_delta = deltas[-1] if deltas else math.nan
        max_delta = max(deltas) if deltas else math.nan
        status = (
            "PASS"
            if (not math.isnan(final_delta) and final_delta <= args.convergence_rel_tol)
            else "FAIL"
        )
        convergence_rows.append(
            {
                "metric": metric,
                "final_rel_delta": final_delta,
                "max_rel_delta": max_delta,
                "status": status,
            }
        )
        if status != "PASS":
            failures.append(f"convergence:{metric}")

    for row in rows:
        case_failures = []
        if row["full_closure_status"] != "PASS":
            case_failures.append("closure")
        if row["full_local_closure_status"] != "PASS":
            case_failures.append("local_closure")
        if args.check_transport_closure and row["full_transport_status"] != "PASS":
            case_failures.append("transport")
        if case_failures:
            failures.append(f"N_w={row['n_modes']}:{'|'.join(case_failures)}")

    highest_row = rows[-1]
    highest_activity_checks = [
        ("final_jw_ew", args.highest_min_jw_ew_abs),
        ("final_s_leak_abs", args.highest_min_s_leak_abs),
        ("final_mixed_ew2", args.highest_min_mixed_ew2),
        ("final_em_leak_w", args.highest_min_abs_em_leak_w),
        ("final_helicity_sub", args.highest_min_abs_helicity_sub),
        ("final_edotb_sub", args.highest_min_abs_edotb_sub),
    ]
    for metric, threshold in highest_activity_checks:
        if threshold <= 0.0:
            continue
        value = highest_row[metric]
        if math.isnan(value) or value < threshold:
            failures.append(
                f"N_w={highest_row['n_modes']}:activity:{metric}<{threshold:.3e}"
            )

    result_status = "PASS" if not failures else "FAIL"

    scan_keys = [
        "n_modes",
        "n_quadrature",
        *CONVERGENCE_METRICS,
        "full_closure_status",
        "full_local_closure_status",
        "full_transport_status",
    ]
    print(",".join(scan_keys))
    for row in rows:
        print(",".join(fmt(row[key]) for key in scan_keys))

    print("metric,final_rel_delta,max_rel_delta,status")
    for row in convergence_rows:
        print(",".join(fmt(row[key]) for key in ("metric", "final_rel_delta", "max_rel_delta", "status")))
    print(f"scan_check_status,{result_status}")
    print(f"scan_check_failures,{ '|'.join(failures) if failures else 'none' }")

    summary_path = output_dir / "nw_convergence_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(scan_keys)
        for row in rows:
            writer.writerow([row[key] for key in scan_keys])
        writer.writerow([])
        writer.writerow(["metric", "final_rel_delta", "max_rel_delta", "status"])
        for row in convergence_rows:
            writer.writerow(
                [row["metric"], row["final_rel_delta"], row["max_rel_delta"], row["status"]]
            )
        writer.writerow([])
        writer.writerow(["scan_check_status", result_status])
        writer.writerow(["scan_check_failures", "|".join(failures) if failures else "none"])

    if args.fail_on_check and failures:
        raise SystemExit(
            "N_w convergence checks failed: " + ", ".join(failures)
        )


if __name__ == "__main__":
    main()
