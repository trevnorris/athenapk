#!/usr/bin/env python3
"""Sweep full-case mode count and summarize spillback diagnostics."""

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
        raise RuntimeError("Could not find controlled/full rows from scan output.")

    reader = csv.DictReader([lines[header_idx], *rows])
    parsed = {row["case"]: row for row in reader}
    if "controlled" not in parsed or "full" not in parsed:
        raise RuntimeError("Missing controlled/full case rows in scan output.")
    return parsed


def parse_float(row, key):
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return math.nan


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "nan"
    return f"{value:.6e}"


def ratio_or_nan(numer, denom):
    if not math.isfinite(numer) or not math.isfinite(denom) or denom == 0.0:
        return math.nan
    return numer / denom


def parse_mode_list(raw):
    values = []
    for token in raw.split(","):
        s = token.strip()
        if not s:
            continue
        n = int(s)
        if n <= 0:
            raise ValueError(f"n_modes must be > 0 (got {n})")
        values.append(n)
    if not values:
        raise ValueError("At least one n_modes value is required.")
    return values


def default_nquad(n_modes):
    # Matches existing defaults: n=1 -> 3, n=4 -> 6.
    return max(3, n_modes + 2)


def run_case(
    scan_script,
    binary,
    workdir,
    case_dir,
    controlled_input,
    full_input,
    n_modes_full,
    n_quad_full,
    controlled_n_modes,
    controlled_n_quad,
    scan_athena_args,
    scan_controlled_args,
    scan_full_args,
):
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
        f"modes4d/n_modes={controlled_n_modes}",
        "--controlled-arg",
        f"modes4d/n_quadrature={controlled_n_quad}",
        "--full-arg",
        f"modes4d/n_modes={n_modes_full}",
        "--full-arg",
        f"modes4d/n_quadrature={n_quad_full}",
    ]
    for arg in scan_athena_args:
        cmd.extend(["--athena-arg", arg])
    for arg in scan_controlled_args:
        cmd.extend(["--controlled-arg", arg])
    for arg in scan_full_args:
        cmd.extend(["--full-arg", arg])

    proc = subprocess.run(cmd, cwd=workdir, text=True, capture_output=True)
    (case_dir / "scan_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (case_dir / "scan_stderr.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            f"scan failed (rc={proc.returncode}); see {case_dir / 'scan_stderr.log'}"
        )
    return parse_scan_csv(proc.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument("--controlled-input", default="inputs/harris_4d_controlled.in")
    parser.add_argument("--full-input", default="inputs/harris_4d_full_symbreak.in")
    parser.add_argument("--output-dir", default="harris_modecount_ablation_outputs")
    parser.add_argument(
        "--mode-counts",
        default="2,4,6",
        help="Comma-separated full-case n_modes values to test",
    )
    parser.add_argument(
        "--controlled-n-modes",
        type=int,
        default=1,
        help="Controlled-case n_modes override (default: 1)",
    )
    parser.add_argument(
        "--controlled-n-quadrature",
        type=int,
        default=3,
        help="Controlled-case n_quadrature override (default: 3)",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue sweep when one mode-count fails",
    )
    parser.add_argument(
        "--athena-arg",
        action="append",
        default=[],
        help="Additional athena override passed to both controlled/full runs",
    )
    parser.add_argument(
        "--controlled-arg",
        action="append",
        default=[],
        help="Additional athena override passed only to controlled run",
    )
    parser.add_argument(
        "--full-arg",
        action="append",
        default=[],
        help="Additional athena override passed only to full run",
    )
    args = parser.parse_args()

    mode_counts = parse_mode_list(args.mode_counts)
    if args.controlled_n_modes <= 0:
        raise SystemExit("--controlled-n-modes must be > 0")
    if args.controlled_n_quadrature <= 0:
        raise SystemExit("--controlled-n-quadrature must be > 0")

    repo_root = Path(__file__).resolve().parents[1]
    scan_script = repo_root / "scripts" / "harris_scan_matrix.py"
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    controlled_input = Path(args.controlled_input).resolve()
    full_input = Path(args.full_input).resolve()
    output_dir_arg = Path(args.output_dir).expanduser()
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (Path.cwd() / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for n_modes in mode_counts:
        n_quad = default_nquad(n_modes)
        case_dir = output_dir / f"nmodes_{n_modes}"
        try:
            parsed = run_case(
                scan_script=scan_script,
                binary=binary,
                workdir=workdir,
                case_dir=case_dir,
                controlled_input=controlled_input,
                full_input=full_input,
                n_modes_full=n_modes,
                n_quad_full=n_quad,
                controlled_n_modes=args.controlled_n_modes,
                controlled_n_quad=args.controlled_n_quadrature,
                scan_athena_args=args.athena_arg,
                scan_controlled_args=args.controlled_arg,
                scan_full_args=args.full_arg,
            )
            controlled = parsed["controlled"]
            full = parsed["full"]
            full_psi0 = parse_float(full, "final_psi0_span")
            controlled_psi0 = parse_float(controlled, "final_psi0_span")
            full_psi_proj = parse_float(full, "final_psi_proj_span")
            controlled_psi_proj = parse_float(controlled, "final_psi_proj_span")
            full_psi_w0 = parse_float(full, "final_psi_w0_span")
            controlled_psi_w0 = parse_float(controlled, "final_psi_w0_span")
            full_w0_minus_proj = parse_float(full, "final_psi_w0_minus_psi_proj_span")
            controlled_w0_minus_proj = parse_float(
                controlled, "final_psi_w0_minus_psi_proj_span"
            )
            row = {
                "n_modes_full": n_modes,
                "n_quadrature_full": n_quad,
                "status": "PASS",
                "ratio_psi0": ratio_or_nan(full_psi0, controlled_psi0),
                "ratio_psi_proj": ratio_or_nan(full_psi_proj, controlled_psi_proj),
                "ratio_psi_w0": ratio_or_nan(full_psi_w0, controlled_psi_w0),
                "ratio_w0_minus_proj": ratio_or_nan(
                    full_w0_minus_proj, controlled_w0_minus_proj
                ),
                "full_final_int_src_mode0_total_abs": parse_float(
                    full, "final_int_src_mode0_total_abs"
                ),
                "full_final_int_div_mode0_total_abs": parse_float(
                    full, "final_int_div_mode0_total_abs"
                ),
                "full_final_jw_ew_abs": abs(parse_float(full, "final_jw_ew")),
                "full_final_s_leak_abs": parse_float(full, "final_s_leak_abs"),
                "full_final_mixed_ew2": parse_float(full, "final_mixed_ew2"),
                "full_corr_psi0_mixed_ew2": parse_float(full, "corr_psi0_mixed_ew2"),
                "full_corr_psiw0_mixed_ew2": parse_float(full, "corr_psiw0_mixed_ew2"),
                "full_closure_status": full.get("closure_status", "N/A"),
                "full_transport_status": full.get("transport_closure_status", "N/A"),
                "full_activity_status": full.get("activity_status", "N/A"),
                "full_correlation_status": full.get("correlation_status", "N/A"),
                "case_dir": str(case_dir),
            }
        except Exception as err:  # noqa: BLE001
            row = {
                "n_modes_full": n_modes,
                "n_quadrature_full": n_quad,
                "status": "RUN_FAIL",
                "ratio_psi0": math.nan,
                "ratio_psi_proj": math.nan,
                "ratio_psi_w0": math.nan,
                "ratio_w0_minus_proj": math.nan,
                "full_final_int_src_mode0_total_abs": math.nan,
                "full_final_int_div_mode0_total_abs": math.nan,
                "full_final_jw_ew_abs": math.nan,
                "full_final_s_leak_abs": math.nan,
                "full_final_mixed_ew2": math.nan,
                "full_corr_psi0_mixed_ew2": math.nan,
                "full_corr_psiw0_mixed_ew2": math.nan,
                "full_closure_status": "N/A",
                "full_transport_status": "N/A",
                "full_activity_status": "N/A",
                "full_correlation_status": "N/A",
                "case_dir": str(case_dir),
                "error": str(err),
            }
            if not args.continue_on_error:
                rows.append(row)
                break
        rows.append(row)

    if not rows:
        raise SystemExit("No mode-count results were produced.")

    fieldnames = [
        "n_modes_full",
        "n_quadrature_full",
        "status",
        "ratio_psi0",
        "ratio_psi_proj",
        "ratio_psi_w0",
        "ratio_w0_minus_proj",
        "full_final_int_src_mode0_total_abs",
        "full_final_int_div_mode0_total_abs",
        "full_final_jw_ew_abs",
        "full_final_s_leak_abs",
        "full_final_mixed_ew2",
        "full_corr_psi0_mixed_ew2",
        "full_corr_psiw0_mixed_ew2",
        "full_closure_status",
        "full_transport_status",
        "full_activity_status",
        "full_correlation_status",
        "case_dir",
        "error",
    ]
    summary_csv = output_dir / "modecount_ablation_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(",".join(fieldnames))
    for row in rows:
        print(",".join(fmt(row.get(k, "")) for k in fieldnames))
    print(f"summary_csv,{summary_csv}")

    any_fail = any(row.get("status") != "PASS" for row in rows)
    print(f"overall_status,{'FAIL' if any_fail else 'PASS'}")


if __name__ == "__main__":
    main()
