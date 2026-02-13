#!/usr/bin/env python3
"""Run matched/point/gaussian projection-kernel sweeps for Harris scans."""

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class KernelCase:
    name: str
    kernel: str
    sigma_factor: float | None


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
        raise RuntimeError("Could not find controlled/full rows in scan output.")

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


def parse_csv_floats(raw):
    values = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(float(token))
    return values


def resolve_input_path(raw_value, repo_root, workdir):
    path = Path(raw_value)
    if path.is_absolute():
        return path.resolve()
    for base in (Path.cwd(), workdir, repo_root):
        candidate = (base / path).resolve()
        if candidate.exists():
            return candidate
    return (repo_root / path).resolve()


def run_case(
    case,
    scan_script,
    binary,
    workdir,
    controlled_input,
    full_input,
    case_dir,
    set_controlled_kernel,
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
    ]
    for arg in scan_athena_args:
        cmd.extend(["--athena-arg", arg])

    # Enforce full-case kernel selection.
    cmd.extend(["--full-arg", f"modes4d/diag_projection_kernel={case.kernel}"])
    if case.sigma_factor is not None:
        cmd.extend(
            ["--full-arg", f"modes4d/diag_projection_sigma_factor={case.sigma_factor}"]
        )

    if set_controlled_kernel:
        cmd.extend(["--controlled-arg", f"modes4d/diag_projection_kernel={case.kernel}"])
        if case.sigma_factor is not None:
            cmd.extend(
                [
                    "--controlled-arg",
                    f"modes4d/diag_projection_sigma_factor={case.sigma_factor}",
                ]
            )

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


def make_cases(include_kernels, gaussian_sigma_factors):
    cases = []
    for kernel in include_kernels:
        k = kernel.strip().lower()
        if k == "matched":
            cases.append(KernelCase(name="matched", kernel="matched", sigma_factor=None))
        elif k == "point":
            cases.append(KernelCase(name="point", kernel="point", sigma_factor=None))
        elif k == "gaussian":
            for sigma in gaussian_sigma_factors:
                label = f"gaussian_sigma_{sigma:g}".replace(".", "p")
                cases.append(KernelCase(name=label, kernel="gaussian", sigma_factor=sigma))
        else:
            raise ValueError(
                f"Unsupported kernel '{kernel}'. Supported: matched, point, gaussian"
            )
    if not cases:
        raise ValueError("No sweep cases were generated.")
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--controlled-input",
        default="inputs/harris_4d_controlled.in",
        help="Controlled input deck",
    )
    parser.add_argument(
        "--full-input",
        default="inputs/harris_4d_full_symbreak.in",
        help="Full input deck",
    )
    parser.add_argument(
        "--output-dir",
        default="harris_projection_kernel_sweep_outputs",
        help="Directory for case subdirs + summary CSV",
    )
    parser.add_argument(
        "--kernels",
        default="matched,point,gaussian",
        help="Comma-separated kernels to include (matched,point,gaussian)",
    )
    parser.add_argument(
        "--gaussian-sigma-factors",
        default="0.5,0.75,1.0,1.5,2.0",
        help="Comma-separated sigma_factor values for gaussian cases",
    )
    parser.add_argument(
        "--set-controlled-kernel",
        action="store_true",
        help="Also set projection kernel/sigma on controlled run (default: full-only)",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue sweep when one case fails",
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

    repo_root = Path(__file__).resolve().parents[1]
    scan_script = repo_root / "scripts" / "harris_scan_matrix.py"
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    controlled_input = resolve_input_path(args.controlled_input, repo_root, workdir)
    full_input = resolve_input_path(args.full_input, repo_root, workdir)
    output_dir_arg = Path(args.output_dir).expanduser()
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (Path.cwd() / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    kernels = [k.strip() for k in args.kernels.split(",") if k.strip()]
    gaussian_sigmas = parse_csv_floats(args.gaussian_sigma_factors)
    cases = make_cases(kernels, gaussian_sigmas)

    rows = []
    for case in cases:
        case_dir = output_dir / case.name
        try:
            parsed = run_case(
                case=case,
                scan_script=scan_script,
                binary=binary,
                workdir=workdir,
                controlled_input=controlled_input,
                full_input=full_input,
                case_dir=case_dir,
                set_controlled_kernel=args.set_controlled_kernel,
                scan_athena_args=args.athena_arg,
                scan_controlled_args=args.controlled_arg,
                scan_full_args=args.full_arg,
            )
            controlled = parsed["controlled"]
            full = parsed["full"]

            controlled_psi0 = parse_float(controlled, "final_psi0_span")
            controlled_psi_proj = parse_float(controlled, "final_psi_proj_span")
            controlled_psi_w0 = parse_float(controlled, "final_psi_w0_span")
            full_psi0 = parse_float(full, "final_psi0_span")
            full_psi_proj = parse_float(full, "final_psi_proj_span")
            full_psi_w0 = parse_float(full, "final_psi_w0_span")

            row = {
                "case_name": case.name,
                "kernel": case.kernel,
                "sigma_factor": case.sigma_factor if case.sigma_factor is not None else math.nan,
                "status": "PASS",
                "controlled_final_psi0_span": controlled_psi0,
                "controlled_final_psi_proj_span": controlled_psi_proj,
                "controlled_final_psi_w0_span": controlled_psi_w0,
                "full_final_psi0_span": full_psi0,
                "full_final_psi_proj_span": full_psi_proj,
                "full_final_psi_proj_even_span": parse_float(
                    full, "final_psi_proj_even_span"
                ),
                "full_final_psi_proj_odd_span": parse_float(
                    full, "final_psi_proj_odd_span"
                ),
                "full_final_psi_proj_matched_span": parse_float(
                    full, "final_psi_proj_matched_span"
                ),
                "full_final_psi_proj_point_span": parse_float(
                    full, "final_psi_proj_point_span"
                ),
                "full_final_psi_proj_gaussian_span": parse_float(
                    full, "final_psi_proj_gaussian_span"
                ),
                "full_final_psi_w0_span": full_psi_w0,
                "full_final_psi_w0_even_span": parse_float(
                    full, "final_psi_w0_even_span"
                ),
                "full_final_psi_w0_odd_span": parse_float(
                    full, "final_psi_w0_odd_span"
                ),
                "full_final_psi_w0_minus_psi_proj_span": parse_float(
                    full, "final_psi_w0_minus_psi_proj_span"
                ),
                "full_final_psi_proj_odd_over_even_span": parse_float(
                    full, "final_psi_proj_odd_over_even_span"
                ),
                "full_final_psi_w0_odd_over_even_span": parse_float(
                    full, "final_psi_w0_odd_over_even_span"
                ),
                "ratio_full_over_controlled_psi0": ratio_or_nan(full_psi0, controlled_psi0),
                "ratio_full_over_controlled_psiproj": ratio_or_nan(
                    full_psi_proj, controlled_psi_proj
                ),
                "ratio_full_over_controlled_psiw0": ratio_or_nan(
                    full_psi_w0, controlled_psi_w0
                ),
                "full_final_em_a2_even": parse_float(full, "final_em_a2_even"),
                "full_final_em_a2_odd": parse_float(full, "final_em_a2_odd"),
                "full_final_em_a2_odd_fraction": parse_float(
                    full, "final_em_a2_odd_fraction"
                ),
                "full_final_jw_l2_even": parse_float(full, "final_jw_l2_even"),
                "full_final_jw_l2_odd": parse_float(full, "final_jw_l2_odd"),
                "full_final_jw_l2_odd_fraction": parse_float(
                    full, "final_jw_l2_odd_fraction"
                ),
                "full_final_emf_vw_c_abs": parse_float(full, "final_emf_vw_c_abs"),
                "full_final_emf_cov_vxb_abs": parse_float(
                    full, "final_emf_cov_vxb_abs"
                ),
                "full_final_jw_l1": parse_float(full, "final_jw_l1"),
                "full_final_ew_l1": parse_float(full, "final_ew_l1"),
                "full_corr_psi0_emf_vw_c_abs": parse_float(
                    full, "corr_psi0_emf_vw_c_abs"
                ),
                "full_corr_psi0_em_a2_odd_fraction": parse_float(
                    full, "corr_psi0_em_a2_odd_fraction"
                ),
                "full_corr_psi0_jw_l2_odd_fraction": parse_float(
                    full, "corr_psi0_jw_l2_odd_fraction"
                ),
                "full_corr_psiw0_emf_vw_c_abs": parse_float(
                    full, "corr_psiw0_emf_vw_c_abs"
                ),
                "full_corr_psiw0_minus_psiproj_mixed_ew2": parse_float(
                    full, "corr_psiw0_minus_psiproj_mixed_ew2"
                ),
                "full_corr_psiw0_minus_psiproj_emf_vw_c_abs": parse_float(
                    full, "corr_psiw0_minus_psiproj_emf_vw_c_abs"
                ),
                "full_corr_psiw0_minus_psiproj_emf_cov_vxb_abs": parse_float(
                    full, "corr_psiw0_minus_psiproj_emf_cov_vxb_abs"
                ),
                "full_closure_status": full.get("closure_status", "N/A"),
                "full_closure_local_mode0_status": full.get(
                    "closure_local_mode0_status", "N/A"
                ),
                "full_transport_status": full.get("transport_closure_status", "N/A"),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
                "full_correlation_status": full.get("correlation_status", "N/A"),
                "full_activity_status": full.get("activity_status", "N/A"),
                "full_topology_status": full.get("topology_status", "N/A"),
                "full_topology_failures": full.get("topology_failures", "n/a"),
                "case_dir": str(case_dir),
            }
        except Exception as err:  # noqa: BLE001
            row = {
                "case_name": case.name,
                "kernel": case.kernel,
                "sigma_factor": case.sigma_factor if case.sigma_factor is not None else math.nan,
                "status": "RUN_FAIL",
                "controlled_final_psi0_span": math.nan,
                "controlled_final_psi_proj_span": math.nan,
                "controlled_final_psi_w0_span": math.nan,
                "full_final_psi0_span": math.nan,
                "full_final_psi_proj_span": math.nan,
                "full_final_psi_proj_even_span": math.nan,
                "full_final_psi_proj_odd_span": math.nan,
                "full_final_psi_proj_matched_span": math.nan,
                "full_final_psi_proj_point_span": math.nan,
                "full_final_psi_proj_gaussian_span": math.nan,
                "full_final_psi_w0_span": math.nan,
                "full_final_psi_w0_even_span": math.nan,
                "full_final_psi_w0_odd_span": math.nan,
                "full_final_psi_w0_minus_psi_proj_span": math.nan,
                "full_final_psi_proj_odd_over_even_span": math.nan,
                "full_final_psi_w0_odd_over_even_span": math.nan,
                "ratio_full_over_controlled_psi0": math.nan,
                "ratio_full_over_controlled_psiproj": math.nan,
                "ratio_full_over_controlled_psiw0": math.nan,
                "full_final_em_a2_even": math.nan,
                "full_final_em_a2_odd": math.nan,
                "full_final_em_a2_odd_fraction": math.nan,
                "full_final_jw_l2_even": math.nan,
                "full_final_jw_l2_odd": math.nan,
                "full_final_jw_l2_odd_fraction": math.nan,
                "full_final_emf_vw_c_abs": math.nan,
                "full_final_emf_cov_vxb_abs": math.nan,
                "full_final_jw_l1": math.nan,
                "full_final_ew_l1": math.nan,
                "full_corr_psi0_emf_vw_c_abs": math.nan,
                "full_corr_psi0_em_a2_odd_fraction": math.nan,
                "full_corr_psi0_jw_l2_odd_fraction": math.nan,
                "full_corr_psiw0_emf_vw_c_abs": math.nan,
                "full_corr_psiw0_minus_psiproj_mixed_ew2": math.nan,
                "full_corr_psiw0_minus_psiproj_emf_vw_c_abs": math.nan,
                "full_corr_psiw0_minus_psiproj_emf_cov_vxb_abs": math.nan,
                "full_closure_status": "N/A",
                "full_closure_local_mode0_status": "N/A",
                "full_transport_status": "N/A",
                "full_em_bulk_ledger_status": "N/A",
                "full_correlation_status": "N/A",
                "full_activity_status": "N/A",
                "full_topology_status": "N/A",
                "full_topology_failures": "n/a",
                "case_dir": str(case_dir),
                "error": str(err),
            }
            if not args.continue_on_error:
                rows.append(row)
                break
        rows.append(row)

    if not rows:
        raise SystemExit("No projection-kernel sweep results were produced.")

    fieldnames = [
        "case_name",
        "kernel",
        "sigma_factor",
        "status",
        "controlled_final_psi0_span",
        "controlled_final_psi_proj_span",
        "controlled_final_psi_w0_span",
        "full_final_psi0_span",
        "full_final_psi_proj_span",
        "full_final_psi_proj_even_span",
        "full_final_psi_proj_odd_span",
        "full_final_psi_proj_matched_span",
        "full_final_psi_proj_point_span",
        "full_final_psi_proj_gaussian_span",
        "full_final_psi_w0_span",
        "full_final_psi_w0_even_span",
        "full_final_psi_w0_odd_span",
        "full_final_psi_w0_minus_psi_proj_span",
        "full_final_psi_proj_odd_over_even_span",
        "full_final_psi_w0_odd_over_even_span",
        "ratio_full_over_controlled_psi0",
        "ratio_full_over_controlled_psiproj",
        "ratio_full_over_controlled_psiw0",
        "full_final_em_a2_even",
        "full_final_em_a2_odd",
        "full_final_em_a2_odd_fraction",
        "full_final_jw_l2_even",
        "full_final_jw_l2_odd",
        "full_final_jw_l2_odd_fraction",
        "full_final_emf_vw_c_abs",
        "full_final_emf_cov_vxb_abs",
        "full_final_jw_l1",
        "full_final_ew_l1",
        "full_corr_psi0_emf_vw_c_abs",
        "full_corr_psi0_em_a2_odd_fraction",
        "full_corr_psi0_jw_l2_odd_fraction",
        "full_corr_psiw0_emf_vw_c_abs",
        "full_corr_psiw0_minus_psiproj_mixed_ew2",
        "full_corr_psiw0_minus_psiproj_emf_vw_c_abs",
        "full_corr_psiw0_minus_psiproj_emf_cov_vxb_abs",
        "full_closure_status",
        "full_closure_local_mode0_status",
        "full_transport_status",
        "full_em_bulk_ledger_status",
        "full_correlation_status",
        "full_activity_status",
        "full_topology_status",
        "full_topology_failures",
        "case_dir",
        "error",
    ]
    summary_csv = output_dir / "projection_kernel_sweep_summary.csv"
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
