#!/usr/bin/env python3
"""Run observer-kernel and N_w convergence sweeps in one matrix."""

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class KernelCase:
    kernel: str
    sigma_factor: float | None


def parse_csv_row_pair(stdout_text):
    lines = [line.strip() for line in stdout_text.splitlines() if line.strip()]
    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("case,"):
            header_idx = idx
            break
    if header_idx is None:
        raise RuntimeError("Could not find harris_scan_matrix CSV header in stdout.")

    rows = []
    for line in lines[header_idx + 1 :]:
        if line.startswith("controlled,") or line.startswith("full,"):
            rows.append(line)
        if len(rows) == 2:
            break
    if len(rows) != 2:
        raise RuntimeError("Could not find controlled/full rows in scan stdout.")

    reader = csv.DictReader([lines[header_idx], *rows])
    parsed = {row["case"]: row for row in reader}
    if "controlled" not in parsed or "full" not in parsed:
        raise RuntimeError("Missing controlled/full rows in parsed CSV output.")
    return parsed


def parse_float(row, key):
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return math.nan


def ratio_or_nan(numer, denom):
    if not math.isfinite(numer) or not math.isfinite(denom) or denom == 0.0:
        return math.nan
    return numer / denom


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "nan"
    return f"{value:.6e}"


def parse_resolution_list(raw):
    values = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        pieces = token.lower().split("x")
        if len(pieces) != 3:
            raise ValueError(f"resolution '{token}' must be in nx1xnx2xnx3 form")
        nx = tuple(int(p) for p in pieces)
        if any(v <= 0 for v in nx):
            raise ValueError(f"resolution '{token}' must use positive integers")
        values.append(nx)
    if not values:
        raise ValueError("at least one resolution is required")
    return values


def parse_meshblock(raw):
    pieces = raw.lower().split("x")
    if len(pieces) != 3:
        raise ValueError("--meshblock must be in nx1xnx2xnx3 form")
    mb = tuple(int(p) for p in pieces)
    if any(v <= 0 for v in mb):
        raise ValueError("--meshblock values must be positive")
    return mb


def parse_modes(raw):
    values = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        val = int(token)
        if val < 2:
            raise ValueError("full n_modes values must be >= 2")
        values.append(val)
    if not values:
        raise ValueError("at least one full n_modes value is required")
    return sorted(set(values))


def parse_kernel_cases(raw_kernels, raw_sigmas):
    kernels = []
    sigma_values = []
    for token in raw_sigmas.split(","):
        token = token.strip()
        if token:
            sigma_values.append(float(token))
    if not sigma_values:
        sigma_values = [1.0]

    for token in raw_kernels.split(","):
        k = token.strip().lower()
        if not k:
            continue
        if k == "matched":
            kernels.append(KernelCase(kernel="matched", sigma_factor=None))
        elif k == "point":
            kernels.append(KernelCase(kernel="point", sigma_factor=None))
        elif k == "gaussian":
            for sigma in sigma_values:
                kernels.append(KernelCase(kernel="gaussian", sigma_factor=sigma))
        else:
            raise ValueError(
                f"unsupported kernel '{k}', expected matched/point/gaussian"
            )
    if not kernels:
        raise ValueError("at least one kernel case is required")
    return kernels


def case_label(kernel_case):
    if kernel_case.sigma_factor is None:
        return kernel_case.kernel
    sigma = str(kernel_case.sigma_factor).replace(".", "p")
    return f"{kernel_case.kernel}_s{sigma}"


def run_scan(
    scan_script,
    binary,
    workdir,
    controlled_input,
    full_input,
    output_dir,
    nx,
    mb,
    full_n_modes,
    controlled_n_modes,
    nq_offset,
    kernel_case,
    tlim,
    nlim,
    history_dt,
    set_controlled_kernel,
    athena_args,
    controlled_args,
    full_args,
    scan_flags,
):
    full_nq = full_n_modes + nq_offset
    controlled_nq = controlled_n_modes + nq_offset
    if controlled_nq < controlled_n_modes:
        raise ValueError("controlled n_quadrature must be >= controlled n_modes")
    if full_nq < full_n_modes:
        raise ValueError("full n_quadrature must be >= full n_modes")

    case_dir = (
        output_dir
        / f"r{nx[0]}x{nx[1]}x{nx[2]}"
        / f"mb{mb[0]}x{mb[1]}x{mb[2]}"
        / f"nw{full_n_modes}"
        / case_label(kernel_case)
    )
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
        "--athena-arg",
        f"parthenon/mesh/nx1={nx[0]}",
        "--athena-arg",
        f"parthenon/mesh/nx2={nx[1]}",
        "--athena-arg",
        f"parthenon/mesh/nx3={nx[2]}",
        "--athena-arg",
        f"parthenon/meshblock/nx1={mb[0]}",
        "--athena-arg",
        f"parthenon/meshblock/nx2={mb[1]}",
        "--athena-arg",
        f"parthenon/meshblock/nx3={mb[2]}",
        "--athena-arg",
        f"parthenon/time/tlim={tlim}",
        "--athena-arg",
        f"parthenon/time/nlim={nlim}",
        "--athena-arg",
        "parthenon/output0/dt=-1",
        "--athena-arg",
        f"parthenon/output1/dt={history_dt}",
        "--controlled-arg",
        f"modes4d/n_modes={controlled_n_modes}",
        "--controlled-arg",
        f"modes4d/n_quadrature={controlled_nq}",
        "--full-arg",
        f"modes4d/n_modes={full_n_modes}",
        "--full-arg",
        f"modes4d/n_quadrature={full_nq}",
        "--full-arg",
        f"modes4d/diag_projection_kernel={kernel_case.kernel}",
    ]
    if kernel_case.sigma_factor is not None:
        cmd.extend(
            [
                "--full-arg",
                f"modes4d/diag_projection_sigma_factor={kernel_case.sigma_factor}",
            ]
        )
    if set_controlled_kernel:
        cmd.extend(
            [
                "--controlled-arg",
                f"modes4d/diag_projection_kernel={kernel_case.kernel}",
            ]
        )
        if kernel_case.sigma_factor is not None:
            cmd.extend(
                [
                    "--controlled-arg",
                    f"modes4d/diag_projection_sigma_factor={kernel_case.sigma_factor}",
                ]
            )

    for arg in athena_args:
        cmd.extend(["--athena-arg", arg])
    for arg in controlled_args:
        cmd.extend(["--controlled-arg", arg])
    for arg in full_args:
        cmd.extend(["--full-arg", arg])
    cmd.extend(scan_flags)

    proc = subprocess.run(cmd, cwd=workdir, text=True, capture_output=True)
    (case_dir / "scan_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (case_dir / "scan_stderr.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            f"scan failed (rc={proc.returncode}); see {case_dir / 'scan_stderr.log'}"
        )
    parsed = parse_csv_row_pair(proc.stdout)
    return case_dir, parsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Runtime working directory")
    parser.add_argument(
        "--controlled-input",
        default="inputs/harris_4d_controlled.in",
        help="Controlled-limit input deck",
    )
    parser.add_argument(
        "--full-input",
        default="inputs/harris_4d_full_symbreak.in",
        help="Full-channel input deck",
    )
    parser.add_argument(
        "--output-dir",
        default="harris_observer_convergence_outputs",
        help="Output directory",
    )
    parser.add_argument(
        "--resolutions",
        default="128x64x64,160x80x80",
        help="Comma-separated list of nx1xnx2xnx3",
    )
    parser.add_argument(
        "--meshblock",
        default="32x16x16",
        help="Meshblock shape nx1xnx2xnx3",
    )
    parser.add_argument(
        "--full-n-modes",
        default="2,4",
        help="Comma-separated N_w values for full runs",
    )
    parser.add_argument(
        "--controlled-n-modes",
        type=int,
        default=1,
        help="N_w for controlled runs",
    )
    parser.add_argument(
        "--nq-offset",
        type=int,
        default=2,
        help="Use n_quadrature = n_modes + nq_offset",
    )
    parser.add_argument(
        "--kernels",
        default="matched,point,gaussian",
        help="Comma-separated kernels to run",
    )
    parser.add_argument(
        "--gaussian-sigma-factors",
        default="0.5,0.75,1.0,1.5",
        help="Sigma factors for gaussian kernel",
    )
    parser.add_argument(
        "--set-controlled-kernel",
        action="store_true",
        help="Apply observer kernel to controlled run as well",
    )
    parser.add_argument("--tlim", type=float, default=0.35, help="Simulation tlim")
    parser.add_argument("--nlim", type=int, default=10000, help="Simulation nlim")
    parser.add_argument(
        "--history-dt",
        type=float,
        default=0.002,
        help="History output cadence (parthenon/output1/dt)",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue after failed cases and include them in summary",
    )
    parser.add_argument(
        "--athena-arg",
        action="append",
        default=[],
        help="Extra athena override forwarded to both runs",
    )
    parser.add_argument(
        "--controlled-arg",
        action="append",
        default=[],
        help="Extra override forwarded to controlled run",
    )
    parser.add_argument(
        "--full-arg",
        action="append",
        default=[],
        help="Extra override forwarded to full run",
    )
    parser.add_argument(
        "--scan-flag",
        action="append",
        default=[],
        help="Raw flag forwarded directly to harris_scan_matrix.py (e.g. --check-topology)",
    )
    args = parser.parse_args()

    if args.controlled_n_modes <= 0:
        raise SystemExit("--controlled-n-modes must be > 0")
    if args.nq_offset < 0:
        raise SystemExit("--nq-offset must be >= 0")

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

    resolutions = parse_resolution_list(args.resolutions)
    meshblock = parse_meshblock(args.meshblock)
    full_n_modes_values = parse_modes(args.full_n_modes)
    kernel_cases = parse_kernel_cases(args.kernels, args.gaussian_sigma_factors)

    rows = []
    for nx in resolutions:
        for nw in full_n_modes_values:
            for kernel_case in kernel_cases:
                try:
                    case_dir, parsed = run_scan(
                        scan_script=scan_script,
                        binary=binary,
                        workdir=workdir,
                        controlled_input=controlled_input,
                        full_input=full_input,
                        output_dir=output_dir,
                        nx=nx,
                        mb=meshblock,
                        full_n_modes=nw,
                        controlled_n_modes=args.controlled_n_modes,
                        nq_offset=args.nq_offset,
                        kernel_case=kernel_case,
                        tlim=args.tlim,
                        nlim=args.nlim,
                        history_dt=args.history_dt,
                        set_controlled_kernel=args.set_controlled_kernel,
                        athena_args=args.athena_arg,
                        controlled_args=args.controlled_arg,
                        full_args=args.full_arg,
                        scan_flags=args.scan_flag,
                    )
                    controlled = parsed["controlled"]
                    full = parsed["full"]
                    row = {
                        "resolution": f"{nx[0]}x{nx[1]}x{nx[2]}",
                        "meshblock": f"{meshblock[0]}x{meshblock[1]}x{meshblock[2]}",
                        "full_n_modes": nw,
                        "kernel": kernel_case.kernel,
                        "sigma_factor": kernel_case.sigma_factor,
                        "status": "PASS",
                        "controlled_final_psi0_span": parse_float(controlled, "final_psi0_span"),
                        "full_final_psi0_span": parse_float(full, "final_psi0_span"),
                        "ratio_full_over_controlled_psi0": ratio_or_nan(
                            parse_float(full, "final_psi0_span"),
                            parse_float(controlled, "final_psi0_span"),
                        ),
                        "controlled_final_psi_proj_span": parse_float(
                            controlled, "final_psi_proj_span"
                        ),
                        "full_final_psi_proj_span": parse_float(full, "final_psi_proj_span"),
                        "ratio_full_over_controlled_psi_proj": ratio_or_nan(
                            parse_float(full, "final_psi_proj_span"),
                            parse_float(controlled, "final_psi_proj_span"),
                        ),
                        "controlled_final_psi_w0_span": parse_float(
                            controlled, "final_psi_w0_span"
                        ),
                        "full_final_psi_w0_span": parse_float(full, "final_psi_w0_span"),
                        "ratio_full_over_controlled_psi_w0": ratio_or_nan(
                            parse_float(full, "final_psi_w0_span"),
                            parse_float(controlled, "final_psi_w0_span"),
                        ),
                        "full_final_psi_w0_minus_psi_proj_span": parse_float(
                            full, "final_psi_w0_minus_psi_proj_span"
                        ),
                        "full_final_mixed_ew2": parse_float(full, "final_mixed_ew2"),
                        "full_final_emf_vw_c_abs": abs(parse_float(full, "final_emf_vw_c_abs")),
                        "full_final_emf_cov_vxb_abs": abs(
                            parse_float(full, "final_emf_cov_vxb_abs")
                        ),
                        "full_final_s_leak_abs": parse_float(full, "final_s_leak_abs"),
                        "full_final_helicity_sub_abs": abs(
                            parse_float(full, "final_helicity_sub")
                        ),
                        "full_final_edotb_sub_abs": abs(parse_float(full, "final_edotb_sub")),
                        "full_corr_psiw0_mixed_ew2": parse_float(full, "corr_psiw0_mixed_ew2"),
                        "full_corr_psiw0_emf_vw_c_abs": parse_float(
                            full, "corr_psiw0_emf_vw_c_abs"
                        ),
                        "full_corr_psiw0_emf_cov_vxb_abs": parse_float(
                            full, "corr_psiw0_emf_cov_vxb_abs"
                        ),
                        "full_corr_psiw0_minus_psiproj_src_mode0_total_abs": parse_float(
                            full, "corr_psiw0_minus_psiproj_src_mode0_total_abs"
                        ),
                        "full_correlation_status": full.get("correlation_status", "N/A"),
                        "full_activity_status": full.get("activity_status", "N/A"),
                        "full_mechanism_status": full.get("mechanism_status", "N/A"),
                        "full_topology_status": full.get("topology_status", "N/A"),
                        "full_closure_status": full.get("closure_status", "N/A"),
                        "full_closure_local_mode0_status": full.get(
                            "closure_local_mode0_status", "N/A"
                        ),
                        "full_transport_status": full.get("transport_closure_status", "N/A"),
                        "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
                        "full_history_nonfinite_status": full.get(
                            "history_nonfinite_status", "N/A"
                        ),
                        "case_dir": str(case_dir),
                    }
                except Exception as exc:  # noqa: BLE001
                    if not args.continue_on_error:
                        raise
                    row = {
                        "resolution": f"{nx[0]}x{nx[1]}x{nx[2]}",
                        "meshblock": f"{meshblock[0]}x{meshblock[1]}x{meshblock[2]}",
                        "full_n_modes": nw,
                        "kernel": kernel_case.kernel,
                        "sigma_factor": kernel_case.sigma_factor,
                        "status": "FAIL",
                        "error": str(exc),
                    }
                rows.append(row)
                if row["status"] == "PASS":
                    print(
                        "case,"
                        f"{row['resolution']},"
                        f"nw={row['full_n_modes']},"
                        f"{row['kernel']},"
                        f"{fmt(row.get('sigma_factor', math.nan))},"
                        f"{row['full_correlation_status']},"
                        f"{row['full_activity_status']},"
                        f"{row['full_mechanism_status']},"
                        f"{row['full_topology_status']},"
                        f"{fmt(row['ratio_full_over_controlled_psi0'])},"
                        f"{fmt(row['ratio_full_over_controlled_psi_w0'])}"
                    )
                else:
                    print(
                        "case,"
                        f"{row['resolution']},"
                        f"nw={row['full_n_modes']},"
                        f"{row['kernel']},"
                        f"{fmt(row.get('sigma_factor', math.nan))},"
                        "FAIL,"
                        f"{row.get('error', 'unknown error')}"
                    )

    fieldnames = sorted({k for row in rows for k in row.keys()})
    summary_csv = output_dir / "observer_convergence_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    report_md = output_dir / "observer_convergence_report.md"
    lines = []
    lines.append("# Observer + Convergence Sweep")
    lines.append("")
    lines.append(f"- summary_csv: `{summary_csv}`")
    lines.append("")
    lines.append(
        "| resolution | meshblock | full_n_modes | kernel | sigma | status | "
        "psi0_ratio | psi_w0_ratio | corr | activity | mechanism | topology | "
        "closure | closure_local | transport | ledger | nonfinite |"
    )
    lines.append(
        "| :--- | :--- | ---: | :--- | ---: | :---: | ---: | ---: | :---: | "
        ":---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |"
    )
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row.get("resolution", "n/a"),
                    row.get("meshblock", "n/a"),
                    str(row.get("full_n_modes", "n/a")),
                    row.get("kernel", "n/a"),
                    fmt(row.get("sigma_factor", math.nan)),
                    row.get("status", "FAIL"),
                    fmt(row.get("ratio_full_over_controlled_psi0", math.nan)),
                    fmt(row.get("ratio_full_over_controlled_psi_w0", math.nan)),
                    row.get("full_correlation_status", "n/a"),
                    row.get("full_activity_status", "n/a"),
                    row.get("full_mechanism_status", "n/a"),
                    row.get("full_topology_status", "n/a"),
                    row.get("full_closure_status", "n/a"),
                    row.get("full_closure_local_mode0_status", "n/a"),
                    row.get("full_transport_status", "n/a"),
                    row.get("full_em_bulk_ledger_status", "n/a"),
                    row.get("full_history_nonfinite_status", "n/a"),
                ]
            )
            + " |"
        )
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    overall_status = "PASS"
    for row in rows:
        if row.get("status") != "PASS":
            overall_status = "FAIL"
            break
    print(f"summary_csv,{summary_csv}")
    print(f"report_md,{report_md}")
    print(f"overall_status,{overall_status}")


if __name__ == "__main__":
    main()
