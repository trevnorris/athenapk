#!/usr/bin/env python3
"""Run controlled/full Harris scans over Lundquist number S."""

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def parse_input_sections(path):
    sections = {}
    current = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            current = line[1:-1].strip()
            sections.setdefault(current, {})
            continue
        if "=" in line and current is not None:
            key, value = line.split("=", 1)
            sections[current][key.strip()] = value.strip()
    return sections


def get_float(sections, section, key, default):
    value = sections.get(section, {}).get(key)
    if value is None:
        return default
    return float(value)


def parse_s_values(raw):
    values = []
    for token in raw.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        value = float(stripped)
        if value <= 0.0:
            raise ValueError(f"Lundquist values must be > 0. Got {value}")
        values.append(value)
    if not values:
        raise ValueError("At least one Lundquist value must be provided.")
    return values


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


def pearson(xs, ys):
    if len(xs) != len(ys) or len(xs) < 2:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) * (x - mx) for x in xs)
    vy = sum((y - my) * (y - my) for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def loglog_slope(xs, ys):
    samples = []
    for x, y in zip(xs, ys):
        if not (math.isfinite(x) and math.isfinite(y)):
            continue
        if x <= 0.0 or y <= 0.0:
            continue
        samples.append((math.log10(x), math.log10(y)))
    if len(samples) < 2:
        return math.nan
    mean_x = sum(x for x, _ in samples) / len(samples)
    mean_y = sum(y for _, y in samples) / len(samples)
    var_x = sum((x - mean_x) * (x - mean_x) for x, _ in samples)
    if var_x <= 0.0:
        return math.nan
    cov_xy = sum((x - mean_x) * (y - mean_y) for x, y in samples)
    return cov_xy / var_x


def finite_values(values):
    return [value for value in values if math.isfinite(value)]


def series_min(values):
    vals = finite_values(values)
    if not vals:
        return math.nan
    return min(vals)


def series_max(values):
    vals = finite_values(values)
    if not vals:
        return math.nan
    return max(vals)


def series_span(values):
    vals = finite_values(values)
    if not vals:
        return math.nan
    return max(vals) - min(vals)


def series_std(values):
    vals = finite_values(values)
    if not vals:
        return math.nan
    mean = sum(vals) / len(vals)
    return math.sqrt(sum((value - mean) * (value - mean) for value in vals) / len(vals))


def safe_ratio(numerator, denominator):
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return math.nan
    if abs(denominator) <= 0.0:
        return math.nan
    return numerator / denominator


def mode_activity_proxy(jw_mode_l2, jw_ew):
    candidates = []
    if math.isfinite(jw_mode_l2):
        candidates.append(jw_mode_l2)
    if math.isfinite(jw_ew):
        candidates.append(abs(jw_ew))
    if not candidates:
        return math.nan
    return max(candidates)


def parse_threshold_specs(specs, label, require_non_negative=False):
    parsed = {}
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"{label} entries must be metric=threshold. Got '{spec}'")
        metric, raw_value = spec.split("=", 1)
        metric = metric.strip()
        if not metric:
            raise ValueError(f"{label} entry has empty metric in '{spec}'")
        threshold = float(raw_value.strip())
        if require_non_negative and threshold < 0.0:
            raise ValueError(f"{label} threshold must be >= 0 for '{metric}'")
        parsed[metric] = threshold
    return parsed


def evaluate_scan_summary(
    scan_summary,
    min_abs_corr,
    min_corr,
    max_corr,
    min_metric,
    min_abs_metric,
    max_metric,
    min_corr_span,
    min_corr_std,
):
    failures = []
    for metric, threshold in min_abs_corr.items():
        value = scan_summary.get(metric, math.nan)
        if math.isnan(value) or abs(value) < threshold:
            failures.append(
                f"{metric}: |{value:.6e}| < {threshold:.6e}"
                if not math.isnan(value)
                else f"{metric}: nan < {threshold:.6e}"
            )
    for metric, threshold in min_corr.items():
        value = scan_summary.get(metric, math.nan)
        if math.isnan(value) or value < threshold:
            failures.append(
                f"{metric}: {value:.6e} < {threshold:.6e}"
                if not math.isnan(value)
                else f"{metric}: nan < {threshold:.6e}"
            )
    for metric, threshold in max_corr.items():
        value = scan_summary.get(metric, math.nan)
        if math.isnan(value) or value > threshold:
            failures.append(
                f"{metric}: {value:.6e} > {threshold:.6e}"
                if not math.isnan(value)
                else f"{metric}: nan > {threshold:.6e}"
            )
    for metric, threshold in min_metric.items():
        value = scan_summary.get(metric, math.nan)
        if math.isnan(value) or value < threshold:
            failures.append(
                f"{metric}: {value:.6e} < {threshold:.6e}"
                if not math.isnan(value)
                else f"{metric}: nan < {threshold:.6e}"
            )
    for metric, threshold in min_abs_metric.items():
        value = scan_summary.get(metric, math.nan)
        if math.isnan(value) or abs(value) < threshold:
            failures.append(
                f"{metric}: |{value:.6e}| < {threshold:.6e}"
                if not math.isnan(value)
                else f"{metric}: |nan| < {threshold:.6e}"
            )
    for metric, threshold in max_metric.items():
        value = scan_summary.get(metric, math.nan)
        if math.isnan(value) or value > threshold:
            failures.append(
                f"{metric}: {value:.6e} > {threshold:.6e}"
                if not math.isnan(value)
                else f"{metric}: nan > {threshold:.6e}"
            )
    for metric, threshold in min_corr_span.items():
        key = f"span_for_{metric}"
        value = scan_summary.get(key, math.nan)
        if math.isnan(value) or value < threshold:
            failures.append(
                f"{key}: {value:.6e} < {threshold:.6e}"
                if not math.isnan(value)
                else f"{key}: nan < {threshold:.6e}"
            )
    for metric, threshold in min_corr_std.items():
        key = f"std_for_{metric}"
        value = scan_summary.get(key, math.nan)
        if math.isnan(value) or value < threshold:
            failures.append(
                f"{key}: {value:.6e} < {threshold:.6e}"
                if not math.isnan(value)
                else f"{key}: nan < {threshold:.6e}"
            )
    if failures:
        return ("FAIL", "|".join(failures))
    return ("PASS", "none")


def resolve_input_path(raw_value, repo_root, workdir):
    path = Path(raw_value)
    if path.is_absolute():
        return path.resolve()
    for base in (Path.cwd(), workdir, repo_root):
        candidate = (base / path).resolve()
        if candidate.exists():
            return candidate
    return (repo_root / path).resolve()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Directory to run simulations in")
    parser.add_argument(
        "--controlled-input",
        default="inputs/harris_4d_controlled.in",
        help="Controlled-limit Harris input deck",
    )
    parser.add_argument(
        "--full-input",
        default="inputs/harris_4d_full.in",
        help="Full-channel Harris input deck",
    )
    parser.add_argument(
        "--output-dir",
        default="harris_lundquist_scan_outputs",
        help="Directory for scan outputs",
    )
    parser.add_argument(
        "--summary-csv",
        default="lundquist_scan_summary.csv",
        help="Filename for summary CSV in output-dir",
    )
    parser.add_argument(
        "--s-values",
        default="250,500,1000,2000",
        help="Comma-separated Lundquist values to scan",
    )
    parser.add_argument(
        "--length-scale",
        type=float,
        default=-1.0,
        help="Characteristic length for S=L*Va/eta; <=0 infers from x1max-x1min",
    )
    parser.add_argument(
        "--alfven-speed",
        type=float,
        default=-1.0,
        help="Alfven speed for S=L*Va/eta; <=0 infers from b0/sqrt(mu0*n_bg)",
    )
    parser.add_argument(
        "--fail-on-check",
        action="store_true",
        help="Propagate fail-on-check to each controlled/full scan run",
    )
    parser.add_argument(
        "--check-transport-closure",
        action="store_true",
        help="Propagate transport-closure gating to each controlled/full scan run",
    )
    parser.add_argument(
        "--scan-arg",
        action="append",
        default=[],
        help="Additional argument forwarded to harris_scan_matrix.py",
    )
    parser.add_argument(
        "--scan-min-abs-corr",
        action="append",
        default=[],
        help="Scan-level gate metric=threshold requiring |corr| >= threshold",
    )
    parser.add_argument(
        "--scan-min-corr",
        action="append",
        default=[],
        help="Scan-level gate metric=threshold requiring corr >= threshold",
    )
    parser.add_argument(
        "--scan-max-corr",
        action="append",
        default=[],
        help="Scan-level gate metric=threshold requiring corr <= threshold",
    )
    parser.add_argument(
        "--scan-min-metric",
        action="append",
        default=[],
        help="Scan-level gate metric=threshold requiring metric >= threshold",
    )
    parser.add_argument(
        "--scan-min-abs-metric",
        action="append",
        default=[],
        help="Scan-level gate metric=threshold requiring |metric| >= threshold",
    )
    parser.add_argument(
        "--scan-max-metric",
        action="append",
        default=[],
        help="Scan-level gate metric=threshold requiring metric <= threshold",
    )
    parser.add_argument(
        "--scan-min-corr-span",
        action="append",
        default=[],
        help="Scan-level gate corr_metric=threshold requiring span_for_<corr_metric> >= threshold",
    )
    parser.add_argument(
        "--scan-min-corr-std",
        action="append",
        default=[],
        help="Scan-level gate corr_metric=threshold requiring std_for_<corr_metric> >= threshold",
    )
    parser.add_argument(
        "--fail-on-scan-check",
        action="store_true",
        help="Exit nonzero when any scan-level correlation gate fails",
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

    s_values = parse_s_values(args.s_values)
    min_abs_corr = parse_threshold_specs(
        args.scan_min_abs_corr, "--scan-min-abs-corr", require_non_negative=True
    )
    min_corr = parse_threshold_specs(args.scan_min_corr, "--scan-min-corr")
    max_corr = parse_threshold_specs(args.scan_max_corr, "--scan-max-corr")
    min_metric = parse_threshold_specs(args.scan_min_metric, "--scan-min-metric")
    min_abs_metric = parse_threshold_specs(
        args.scan_min_abs_metric, "--scan-min-abs-metric", require_non_negative=True
    )
    max_metric = parse_threshold_specs(args.scan_max_metric, "--scan-max-metric")
    min_corr_span = parse_threshold_specs(
        args.scan_min_corr_span, "--scan-min-corr-span", require_non_negative=True
    )
    min_corr_std = parse_threshold_specs(
        args.scan_min_corr_std, "--scan-min-corr-std", require_non_negative=True
    )

    sections = parse_input_sections(full_input)
    x1min = get_float(sections, "parthenon/mesh", "x1min", -1.0)
    x1max = get_float(sections, "parthenon/mesh", "x1max", 1.0)
    inferred_l = x1max - x1min
    if inferred_l <= 0.0:
        raise RuntimeError("Could not infer positive length scale from full input deck.")
    length_scale = args.length_scale if args.length_scale > 0.0 else inferred_l

    b0 = get_float(sections, "problem/harris_4d", "b0", 1.0)
    n_bg = get_float(sections, "problem/harris_4d", "n_bg", 0.2)
    mu0 = get_float(sections, "modes4d", "em_mu0", 1.0)
    if n_bg <= 0.0 or mu0 <= 0.0:
        raise RuntimeError("Cannot infer Alfven speed with non-positive n_bg or em_mu0.")
    inferred_va = abs(b0) / math.sqrt(mu0 * n_bg)
    alfven_speed = args.alfven_speed if args.alfven_speed > 0.0 else inferred_va

    rows = []
    for s_value in s_values:
        eta = (length_scale * alfven_speed) / s_value
        case_tag = str(s_value).replace(".", "p")
        case_outdir = output_dir / f"S_{case_tag}"
        case_outdir.mkdir(parents=True, exist_ok=True)

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
            str(case_outdir),
            "--athena-arg",
            "diffusion/resistivity=ohmic",
            "--athena-arg",
            "diffusion/resistivity_coeff=fixed",
            "--athena-arg",
            "diffusion/integrator=rkl2",
            "--athena-arg",
            "diffusion/rkl2_max_dt_ratio=100.0",
            "--athena-arg",
            f"diffusion/ohm_diff_coeff_code={eta:.12e}",
        ]
        if args.fail_on_check:
            cmd.append("--fail-on-check")
        if args.check_transport_closure:
            cmd.append("--check-transport-closure")
        for scan_arg in args.scan_arg:
            cmd.append(scan_arg)

        proc = subprocess.run(
            cmd,
            cwd=workdir,
            text=True,
            capture_output=True,
        )
        (case_outdir / "scan_stdout.log").write_text(proc.stdout, encoding="utf-8")
        (case_outdir / "scan_stderr.log").write_text(proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise RuntimeError(
                f"Lundquist scan failed for S={s_value} (eta={eta:.6e}), "
                f"returncode={proc.returncode}. See {case_outdir}/scan_stderr.log"
            )

        parsed = parse_scan_csv(proc.stdout)
        full_row = parsed["full"]
        controlled_row = parsed["controlled"]
        controlled_final_jw_ew = parse_float(controlled_row, "final_jw_ew")
        controlled_final_s_leak_abs = parse_float(controlled_row, "final_s_leak_abs")
        controlled_final_mixed_ew2 = parse_float(controlled_row, "final_mixed_ew2")
        controlled_final_mixed_c2 = parse_float(controlled_row, "final_mixed_c2")
        controlled_final_em_leak_w = parse_float(controlled_row, "final_em_leak_w")
        controlled_final_helicity_sub = parse_float(controlled_row, "final_helicity_sub")
        controlled_final_edotb_sub = parse_float(controlled_row, "final_edotb_sub")
        controlled_final_jw_mode_l2_1 = parse_float(controlled_row, "final_jw_mode_l2_1")
        controlled_final_jw_mode_activity_proxy = mode_activity_proxy(
            controlled_final_jw_mode_l2_1, controlled_final_jw_ew
        )
        controlled_max_abs_dpsi0_dt = parse_float(controlled_row, "max_abs_dpsi0_dt")
        controlled_max_dpsi0_dt = parse_float(controlled_row, "max_dpsi0_dt")
        controlled_final_psi_proj_span = parse_float(controlled_row, "final_psi_proj_span")
        controlled_max_abs_dpsi_proj_dt = parse_float(controlled_row, "max_abs_dpsi_proj_dt")
        controlled_max_dpsi_proj_dt = parse_float(controlled_row, "max_dpsi_proj_dt")
        controlled_final_psi_w0_span = parse_float(controlled_row, "final_psi_w0_span")
        controlled_max_abs_dpsi_w0_dt = parse_float(controlled_row, "max_abs_dpsi_w0_dt")
        controlled_max_dpsi_w0_dt = parse_float(controlled_row, "max_dpsi_w0_dt")
        full_final_jw_ew = parse_float(full_row, "final_jw_ew")
        full_final_s_leak_abs = parse_float(full_row, "final_s_leak_abs")
        full_final_mixed_ew2 = parse_float(full_row, "final_mixed_ew2")
        full_final_mixed_c2 = parse_float(full_row, "final_mixed_c2")
        full_final_em_leak_w = parse_float(full_row, "final_em_leak_w")
        full_final_helicity_sub = parse_float(full_row, "final_helicity_sub")
        full_final_edotb_sub = parse_float(full_row, "final_edotb_sub")
        full_final_jw_mode_l2_1 = parse_float(full_row, "final_jw_mode_l2_1")
        full_final_jw_mode_activity_proxy = mode_activity_proxy(
            full_final_jw_mode_l2_1, full_final_jw_ew
        )
        full_max_abs_dpsi0_dt = parse_float(full_row, "max_abs_dpsi0_dt")
        full_max_dpsi0_dt = parse_float(full_row, "max_dpsi0_dt")
        full_final_psi_proj_span = parse_float(full_row, "final_psi_proj_span")
        full_max_abs_dpsi_proj_dt = parse_float(full_row, "max_abs_dpsi_proj_dt")
        full_max_dpsi_proj_dt = parse_float(full_row, "max_dpsi_proj_dt")
        full_final_psi_w0_span = parse_float(full_row, "final_psi_w0_span")
        full_max_abs_dpsi_w0_dt = parse_float(full_row, "max_abs_dpsi_w0_dt")
        full_max_dpsi_w0_dt = parse_float(full_row, "max_dpsi_w0_dt")
        controlled_final_psi0_span = parse_float(controlled_row, "final_psi0_span")
        full_final_psi0_span = parse_float(full_row, "final_psi0_span")
        rows.append(
            {
                "S": s_value,
                "eta": eta,
                "controlled_final_psi0_span": controlled_final_psi0_span,
                "controlled_final_jw_ew": controlled_final_jw_ew,
                "controlled_final_s_leak_abs": controlled_final_s_leak_abs,
                "controlled_final_mixed_ew2": controlled_final_mixed_ew2,
                "controlled_final_mixed_c2": controlled_final_mixed_c2,
                "controlled_final_em_leak_w": controlled_final_em_leak_w,
                "controlled_final_helicity_sub": controlled_final_helicity_sub,
                "controlled_final_edotb_sub": controlled_final_edotb_sub,
                "controlled_final_jw_mode_l2_1": controlled_final_jw_mode_l2_1,
                "controlled_final_jw_mode_activity_proxy": controlled_final_jw_mode_activity_proxy,
                "controlled_max_abs_dpsi0_dt": controlled_max_abs_dpsi0_dt,
                "controlled_max_dpsi0_dt": controlled_max_dpsi0_dt,
                "controlled_final_psi_proj_span": controlled_final_psi_proj_span,
                "controlled_max_abs_dpsi_proj_dt": controlled_max_abs_dpsi_proj_dt,
                "controlled_max_dpsi_proj_dt": controlled_max_dpsi_proj_dt,
                "controlled_final_psi_w0_span": controlled_final_psi_w0_span,
                "controlled_max_abs_dpsi_w0_dt": controlled_max_abs_dpsi_w0_dt,
                "controlled_max_dpsi_w0_dt": controlled_max_dpsi_w0_dt,
                "full_final_psi0_span": full_final_psi0_span,
                "full_final_psi_proj_span": full_final_psi_proj_span,
                "full_final_psi_w0_span": full_final_psi_w0_span,
                "full_final_jw_ew": full_final_jw_ew,
                "full_final_s_leak_abs": full_final_s_leak_abs,
                "full_final_mixed_ew2": full_final_mixed_ew2,
                "full_final_mixed_c2": full_final_mixed_c2,
                "full_final_em_leak_w": full_final_em_leak_w,
                "full_final_helicity_sub": full_final_helicity_sub,
                "full_final_edotb_sub": full_final_edotb_sub,
                "full_final_em_a2_mode_1": parse_float(full_row, "final_em_a2_mode_1"),
                "full_final_em_pi2_mode_1": parse_float(full_row, "final_em_pi2_mode_1"),
                "full_final_jw_mode_l2_1": full_final_jw_mode_l2_1,
                "full_final_jw_mode_activity_proxy": full_final_jw_mode_activity_proxy,
                "full_max_abs_dpsi0_dt": full_max_abs_dpsi0_dt,
                "full_max_dpsi0_dt": full_max_dpsi0_dt,
                "full_max_abs_dpsi_proj_dt": full_max_abs_dpsi_proj_dt,
                "full_max_dpsi_proj_dt": full_max_dpsi_proj_dt,
                "full_max_abs_dpsi_w0_dt": full_max_abs_dpsi_w0_dt,
                "full_max_dpsi_w0_dt": full_max_dpsi_w0_dt,
                "ratio_full_over_controlled_psi0_span": safe_ratio(
                    full_final_psi0_span, controlled_final_psi0_span
                ),
                "ratio_full_over_controlled_psi_proj_span": safe_ratio(
                    full_final_psi_proj_span, controlled_final_psi_proj_span
                ),
                "ratio_full_over_controlled_psi_w0_span": safe_ratio(
                    full_final_psi_w0_span, controlled_final_psi_w0_span
                ),
                "ratio_full_over_controlled_jw_ew_abs": safe_ratio(
                    abs(full_final_jw_ew), abs(controlled_final_jw_ew)
                ),
                "ratio_full_over_controlled_s_leak_abs": safe_ratio(
                    full_final_s_leak_abs, controlled_final_s_leak_abs
                ),
                "ratio_full_over_controlled_mixed_ew2": safe_ratio(
                    full_final_mixed_ew2, controlled_final_mixed_ew2
                ),
                "ratio_full_over_controlled_mixed_c2": safe_ratio(
                    full_final_mixed_c2, controlled_final_mixed_c2
                ),
                "ratio_full_over_controlled_max_abs_dpsi0_dt": safe_ratio(
                    full_max_abs_dpsi0_dt, controlled_max_abs_dpsi0_dt
                ),
                "ratio_full_over_controlled_max_abs_dpsi_proj_dt": safe_ratio(
                    full_max_abs_dpsi_proj_dt, controlled_max_abs_dpsi_proj_dt
                ),
                "ratio_full_over_controlled_max_abs_dpsi_w0_dt": safe_ratio(
                    full_max_abs_dpsi_w0_dt, controlled_max_abs_dpsi_w0_dt
                ),
                "full_corr_psi0_s_leak_abs": parse_float(full_row, "corr_psi0_s_leak_abs"),
                "full_corr_psi0_jw_ew": parse_float(full_row, "corr_psi0_jw_ew"),
                "full_corr_psi0_mixed_ew2": parse_float(full_row, "corr_psi0_mixed_ew2"),
                "full_corr_psi0_em_leak_w": parse_float(full_row, "corr_psi0_em_leak_w"),
                "full_corr_psi0_helicity_sub": parse_float(full_row, "corr_psi0_helicity_sub"),
                "full_corr_psi0_edotb_sub": parse_float(full_row, "corr_psi0_edotb_sub"),
                "full_corr_psi0_jw_mode_l2_1": parse_float(full_row, "corr_psi0_jw_mode_l2_1"),
                "full_corr_psiproj_s_leak_abs": parse_float(full_row, "corr_psiproj_s_leak_abs"),
                "full_corr_psiproj_jw_ew": parse_float(full_row, "corr_psiproj_jw_ew"),
                "full_corr_psiproj_mixed_ew2": parse_float(full_row, "corr_psiproj_mixed_ew2"),
                "full_closure_status": full_row.get("closure_status", "N/A"),
                "full_transport_closure_status": full_row.get(
                    "transport_closure_status", "N/A"
                ),
                "full_correlation_status": full_row.get("correlation_status", "N/A"),
                "full_activity_status": full_row.get("activity_status", "N/A"),
            }
        )

    log_s = [math.log10(r["S"]) for r in rows]
    corr_series = {
        "corr_logS_full_final_psi0_span": [r["full_final_psi0_span"] for r in rows],
        "corr_logS_full_final_psi_proj_span": [r["full_final_psi_proj_span"] for r in rows],
        "corr_logS_full_final_psi_w0_span": [r["full_final_psi_w0_span"] for r in rows],
        "corr_logS_full_final_s_leak_abs": [r["full_final_s_leak_abs"] for r in rows],
        "corr_logS_full_final_jw_ew_abs": [abs(r["full_final_jw_ew"]) for r in rows],
        "corr_logS_full_final_mixed_ew2": [r["full_final_mixed_ew2"] for r in rows],
        "corr_logS_full_final_em_leak_w_abs": [abs(r["full_final_em_leak_w"]) for r in rows],
        "corr_logS_full_final_helicity_sub_abs": [
            abs(r["full_final_helicity_sub"]) for r in rows
        ],
        "corr_logS_full_final_edotb_sub_abs": [abs(r["full_final_edotb_sub"]) for r in rows],
        "corr_logS_full_final_em_a2_mode_1": [r["full_final_em_a2_mode_1"] for r in rows],
        "corr_logS_full_final_em_pi2_mode_1": [r["full_final_em_pi2_mode_1"] for r in rows],
        "corr_logS_full_final_jw_mode_l2_1": [r["full_final_jw_mode_l2_1"] for r in rows],
        "corr_logS_full_final_jw_mode_activity_proxy": [
            r["full_final_jw_mode_activity_proxy"] for r in rows
        ],
        "corr_logS_full_max_abs_dpsi0_dt": [r["full_max_abs_dpsi0_dt"] for r in rows],
        "corr_logS_full_max_dpsi0_dt": [r["full_max_dpsi0_dt"] for r in rows],
        "corr_logS_full_max_abs_dpsi_proj_dt": [r["full_max_abs_dpsi_proj_dt"] for r in rows],
        "corr_logS_full_max_dpsi_proj_dt": [r["full_max_dpsi_proj_dt"] for r in rows],
        "corr_logS_full_max_abs_dpsi_w0_dt": [r["full_max_abs_dpsi_w0_dt"] for r in rows],
        "corr_logS_full_max_dpsi_w0_dt": [r["full_max_dpsi_w0_dt"] for r in rows],
    }
    scan_summary = {}
    for metric, values in corr_series.items():
        scan_summary[metric] = pearson(log_s, values)
        scan_summary[f"span_for_{metric}"] = series_span(values)
        scan_summary[f"std_for_{metric}"] = series_std(values)
    metric_series = {
        "full_final_psi0_span": [r["full_final_psi0_span"] for r in rows],
        "full_final_psi_proj_span": [r["full_final_psi_proj_span"] for r in rows],
        "full_final_psi_w0_span": [r["full_final_psi_w0_span"] for r in rows],
        "full_final_s_leak_abs": [r["full_final_s_leak_abs"] for r in rows],
        "full_final_jw_ew_abs": [abs(r["full_final_jw_ew"]) for r in rows],
        "full_final_mixed_ew2": [r["full_final_mixed_ew2"] for r in rows],
        "full_final_mixed_c2": [r["full_final_mixed_c2"] for r in rows],
        "full_final_em_leak_w_abs": [abs(r["full_final_em_leak_w"]) for r in rows],
        "full_final_helicity_sub_abs": [abs(r["full_final_helicity_sub"]) for r in rows],
        "full_final_edotb_sub_abs": [abs(r["full_final_edotb_sub"]) for r in rows],
        "full_final_jw_mode_activity_proxy": [
            r["full_final_jw_mode_activity_proxy"] for r in rows
        ],
        "full_max_abs_dpsi0_dt": [r["full_max_abs_dpsi0_dt"] for r in rows],
        "full_max_dpsi0_dt": [r["full_max_dpsi0_dt"] for r in rows],
        "full_max_abs_dpsi_proj_dt": [r["full_max_abs_dpsi_proj_dt"] for r in rows],
        "full_max_dpsi_proj_dt": [r["full_max_dpsi_proj_dt"] for r in rows],
        "full_max_abs_dpsi_w0_dt": [r["full_max_abs_dpsi_w0_dt"] for r in rows],
        "full_max_dpsi_w0_dt": [r["full_max_dpsi_w0_dt"] for r in rows],
    }
    for metric, values in metric_series.items():
        scan_summary[f"min_{metric}"] = series_min(values)
        scan_summary[f"max_{metric}"] = series_max(values)
        scan_summary[f"span_{metric}"] = series_span(values)
        scan_summary[f"std_{metric}"] = series_std(values)
    slope_metrics = {
        "slope_logS_controlled_max_abs_dpsi0_dt": loglog_slope(
            [r["S"] for r in rows],
            [r["controlled_max_abs_dpsi0_dt"] for r in rows],
        ),
        "slope_logS_full_max_abs_dpsi0_dt": loglog_slope(
            [r["S"] for r in rows],
            [r["full_max_abs_dpsi0_dt"] for r in rows],
        ),
        "slope_logS_controlled_max_abs_dpsi_w0_dt": loglog_slope(
            [r["S"] for r in rows],
            [r["controlled_max_abs_dpsi_w0_dt"] for r in rows],
        ),
        "slope_logS_full_max_abs_dpsi_w0_dt": loglog_slope(
            [r["S"] for r in rows],
            [r["full_max_abs_dpsi_w0_dt"] for r in rows],
        ),
    }
    scan_summary.update(slope_metrics)
    scan_summary["slope_delta_full_minus_controlled_max_abs_dpsi0_dt"] = (
        slope_metrics["slope_logS_full_max_abs_dpsi0_dt"]
        - slope_metrics["slope_logS_controlled_max_abs_dpsi0_dt"]
        if (
            math.isfinite(slope_metrics["slope_logS_full_max_abs_dpsi0_dt"])
            and math.isfinite(slope_metrics["slope_logS_controlled_max_abs_dpsi0_dt"])
        )
        else math.nan
    )
    scan_summary["slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt"] = (
        slope_metrics["slope_logS_full_max_abs_dpsi_w0_dt"]
        - slope_metrics["slope_logS_controlled_max_abs_dpsi_w0_dt"]
        if (
            math.isfinite(slope_metrics["slope_logS_full_max_abs_dpsi_w0_dt"])
            and math.isfinite(slope_metrics["slope_logS_controlled_max_abs_dpsi_w0_dt"])
        )
        else math.nan
    )
    scan_check_status, scan_check_failures = evaluate_scan_summary(
        scan_summary,
        min_abs_corr,
        min_corr,
        max_corr,
        min_metric,
        min_abs_metric,
        max_metric,
        min_corr_span,
        min_corr_std,
    )

    summary_csv = output_dir / args.summary_csv
    fieldnames = list(rows[0].keys())
    with summary_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(",".join(fieldnames))
    for row in rows:
        print(",".join(f"{row[k]:.6e}" if isinstance(row[k], float) else str(row[k]) for k in fieldnames))
    print(
        "scan_summary,"
        + ",".join(f"{k}={v:.6e}" for k, v in scan_summary.items())
    )
    print(f"scan_check_status,{scan_check_status}")
    print(f"scan_check_failures,{scan_check_failures}")
    print(f"summary_csv,{summary_csv}")
    print(
        "scan_settings,"
        f"L={length_scale:.6e},Va={alfven_speed:.6e},b0={b0:.6e},n_bg={n_bg:.6e},mu0={mu0:.6e}"
    )
    if args.fail_on_scan_check and scan_check_status == "FAIL":
        raise RuntimeError(
            "Lundquist scan-level checks failed: "
            f"{scan_check_failures}"
        )


if __name__ == "__main__":
    main()
