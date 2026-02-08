#!/usr/bin/env python3
"""Run deterministic phase-offset ensembles for Harris Lundquist production scans."""

import argparse
import csv
import math
import shlex
import statistics
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class SeedCase:
    seed: int
    phase_x: float
    phase_z: float

    @property
    def case_id(self) -> str:
        return f"seed_{self.seed:04d}"


def parse_int_csv(raw: str) -> list[int]:
    values = []
    for token in raw.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        values.append(int(stripped))
    if not values:
        raise ValueError("At least one seed must be provided.")
    return values


def phase_from_seed(seed: int, axis: int) -> float:
    # SplitMix64 hashing gives well-separated phases for nearby integer seeds.
    mask = (1 << 64) - 1
    z = (seed + ((axis + 1) * 0x9E3779B97F4A7C15)) & mask
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9 & mask
    z = (z ^ (z >> 27)) * 0x94D049BB133111EB & mask
    z = z ^ (z >> 31)
    return (z / float(1 << 64)) * (2.0 * math.pi)


def parse_kv_lines(text: str) -> dict[str, str]:
    out = {}
    for line in text.splitlines():
        if "," not in line:
            continue
        key, value = line.split(",", 1)
        out[key.strip()] = value.strip()
    return out


def parse_failure_reason(stdout: str, stderr: str, kv: dict[str, str]) -> str:
    stderr_lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    if stderr_lines:
        return stderr_lines[-1]
    scan_failure = kv.get("scan_check_failures", "")
    if scan_failure:
        return scan_failure
    stdout_lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if stdout_lines:
        return stdout_lines[-1]
    return "unknown_failure"


def parse_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value is None:
        return math.nan
    text = str(value).strip()
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def series_min(values: list[float]) -> float:
    finite = [v for v in values if math.isfinite(v)]
    return min(finite) if finite else math.nan


def slope_loglog(xs: list[float], ys: list[float]) -> float:
    pts = []
    for x, y in zip(xs, ys):
        if not (math.isfinite(x) and math.isfinite(y)):
            continue
        if x <= 0.0 or y <= 0.0:
            continue
        pts.append((math.log10(x), math.log10(y)))
    if len(pts) < 2:
        return math.nan
    mean_x = sum(x for x, _ in pts) / len(pts)
    mean_y = sum(y for _, y in pts) / len(pts)
    var_x = sum((x - mean_x) * (x - mean_x) for x, _ in pts)
    if var_x <= 0.0:
        return math.nan
    cov_xy = sum((x - mean_x) * (y - mean_y) for x, y in pts)
    return cov_xy / var_x


def fmt(v: float) -> str:
    return f"{v:.6e}" if math.isfinite(v) else "nan"


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise RuntimeError("No rows to write.")
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    pass_rows = [r for r in rows if r["status"] == "PASS"]

    def floats(key: str, subset: list[dict[str, str]]) -> list[float]:
        return [float(r[key]) for r in subset if r.get(key, "nan") not in ("", "nan")]

    ratio_vals = floats("min_ratio_psi_w0_span", pass_rows)
    rate_ratio_vals = floats("min_ratio_max_abs_dpsi_w0_dt", pass_rows)
    slope_vals = floats("slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt", pass_rows)
    jw_vals = floats("min_full_jw_ew_abs", pass_rows)

    lines = []
    lines.append("# Harris Lundquist Phase Ensemble")
    lines.append("")
    lines.append(f"- Generated (UTC): {datetime.now(timezone.utc).isoformat()}")
    lines.append(f"- Total seeds: **{len(rows)}**")
    lines.append(f"- Passing seeds: **{len(pass_rows)}**")
    lines.append("")
    if pass_rows:
        lines.append("## Passing-Seed Aggregate")
        lines.append("")
        lines.append(f"- min(full/control psi_w0 span) mean: `{fmt(statistics.fmean(ratio_vals))}`")
        lines.append(
            f"- min(full/control max|dpsi_w0/dt|) mean: `{fmt(statistics.fmean(rate_ratio_vals))}`"
        )
        lines.append(
            "- slope delta (full-controlled, dpsi_w0) mean: "
            f"`{fmt(statistics.fmean(slope_vals))}`"
        )
        lines.append(f"- min(|JwEw|) mean: `{fmt(statistics.fmean(jw_vals))}`")
        lines.append("")
    lines.append("## Seed Matrix")
    lines.append("")
    lines.append(
        "| seed | status | phase_x | phase_z | min psi_w0 ratio | min w0 rate ratio | "
        "slope delta dpsi_w0 | min |JwEw| | failure |"
    )
    lines.append(
        "| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |"
    )
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["seed"],
                    row["status"],
                    row["phase_x"],
                    row["phase_z"],
                    row["min_ratio_psi_w0_span"],
                    row["min_ratio_max_abs_dpsi_w0_dt"],
                    row["slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt"],
                    row["min_full_jw_ew_abs"],
                    row["failure_reason"],
                ]
            )
            + " |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--controlled-input", default="inputs/harris_4d_controlled.in")
    parser.add_argument("--full-input", default="inputs/harris_4d_full.in")
    parser.add_argument("--output-dir", default="phase_ensemble_outputs")
    parser.add_argument("--s-values", default="250,500,1000,2000")
    parser.add_argument("--tlim", type=float, default=0.20)
    parser.add_argument("--nlim", type=int, default=2000)
    parser.add_argument("--output-dt", type=float, default=0.02)
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--scan-arg", action="append", default=[])
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--skip-topology", action="store_true")
    parser.add_argument("--fail-on-any-fail", action="store_true")
    args = parser.parse_args()

    if args.tlim <= 0.0:
        raise ValueError("--tlim must be > 0")
    if args.nlim <= 0:
        raise ValueError("--nlim must be > 0")
    if args.output_dt <= 0.0:
        raise ValueError("--output-dt must be > 0")

    repo_root = Path(__file__).resolve().parents[1]
    production_script = repo_root / "scripts" / "run_harris_lundquist_production.sh"
    binary = Path(args.binary).resolve()
    controlled_input = Path(args.controlled_input).resolve()
    full_input = Path(args.full_input).resolve()
    workdir = Path(args.workdir).resolve()
    output_dir_arg = Path(args.output_dir)
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (workdir / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    (workdir / "runs").mkdir(parents=True, exist_ok=True)

    seeds = parse_int_csv(args.seeds)
    cases = [SeedCase(seed=s, phase_x=phase_from_seed(s, 0), phase_z=phase_from_seed(s, 1)) for s in seeds]
    rows = []
    for case in cases:
        case_workdir = workdir / "runs" / case.case_id
        case_workdir.mkdir(parents=True, exist_ok=True)
        cmd = [
            str(production_script),
            "--binary",
            str(binary),
            "--workdir",
            str(case_workdir),
            "--controlled-input",
            str(controlled_input),
            "--full-input",
            str(full_input),
            "--s-values",
            args.s_values,
            "--tlim",
            f"{args.tlim:.6g}",
            "--nlim",
            str(args.nlim),
            "--output-dt",
            f"{args.output_dt:.6g}",
            "--scan-arg",
            "--athena-arg",
            "--scan-arg",
            f"problem/harris_4d/perturb_phase_x={case.phase_x:.12g}",
            "--scan-arg",
            "--athena-arg",
            "--scan-arg",
            f"problem/harris_4d/perturb_phase_z={case.phase_z:.12g}",
        ]
        for extra in args.scan_arg:
            cmd.extend(["--scan-arg", extra])
        if args.skip_plots:
            cmd.append("--skip-plots")
        if args.skip_topology:
            cmd.append("--skip-topology")

        proc = subprocess.run(cmd, cwd=workdir, text=True, capture_output=True)
        stdout_log = case_workdir / "phase_ensemble_stdout.log"
        stderr_log = case_workdir / "phase_ensemble_stderr.log"
        stdout_log.write_text(proc.stdout, encoding="utf-8")
        stderr_log.write_text(proc.stderr, encoding="utf-8")

        kv = parse_kv_lines(proc.stdout)
        summary_csv = Path(kv.get("summary_csv", case_workdir / "outputs" / "lundquist_scan_summary.csv"))
        if not summary_csv.is_absolute():
            summary_csv = (workdir / summary_csv).resolve()

        min_ratio_psi_w0_span = math.nan
        min_ratio_rate_w0 = math.nan
        min_full_jw_ew_abs = math.nan
        slope_delta_w0 = math.nan
        if summary_csv.exists():
            summary_rows = list(csv.DictReader(summary_csv.open("r", encoding="utf-8")))
            min_ratio_psi_w0_span = series_min(
                [parse_float(r, "ratio_full_over_controlled_psi_w0_span") for r in summary_rows]
            )
            min_ratio_rate_w0 = series_min(
                [parse_float(r, "ratio_full_over_controlled_max_abs_dpsi_w0_dt") for r in summary_rows]
            )
            min_full_jw_ew_abs = series_min(
                [abs(parse_float(r, "full_final_jw_ew")) for r in summary_rows]
            )
            s_vals = [parse_float(r, "S") for r in summary_rows]
            ctrl_vals = [parse_float(r, "controlled_max_abs_dpsi_w0_dt") for r in summary_rows]
            full_vals = [parse_float(r, "full_max_abs_dpsi_w0_dt") for r in summary_rows]
            slope_delta_w0 = slope_loglog(s_vals, full_vals) - slope_loglog(s_vals, ctrl_vals)

        overall = kv.get("overall_status", "UNKNOWN")
        status = "PASS" if (proc.returncode == 0 and overall == "PASS") else "FAIL"
        failure_reason = "none" if status == "PASS" else parse_failure_reason(proc.stdout, proc.stderr, kv)
        rows.append(
            {
                "seed": str(case.seed),
                "status": status,
                "returncode": str(proc.returncode),
                "overall_status": overall,
                "failure_reason": failure_reason,
                "phase_x": f"{case.phase_x:.6e}",
                "phase_z": f"{case.phase_z:.6e}",
                "min_ratio_psi_w0_span": fmt(min_ratio_psi_w0_span),
                "min_ratio_max_abs_dpsi_w0_dt": fmt(min_ratio_rate_w0),
                "min_full_jw_ew_abs": fmt(min_full_jw_ew_abs),
                "slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt": fmt(slope_delta_w0),
                "summary_csv": str(summary_csv),
                "production_report": kv.get("production_report", ""),
                "stdout_log": str(stdout_log),
                "stderr_log": str(stderr_log),
                "cmd": " ".join(shlex.quote(piece) for piece in cmd),
            }
        )
        print(
            ",".join(
                [
                    "phase_case",
                    str(case.seed),
                    status,
                    overall,
                    f"{case.phase_x:.6e}",
                    f"{case.phase_z:.6e}",
                    fmt(min_ratio_psi_w0_span),
                ]
            )
        )

    summary_path = output_dir / "phase_ensemble_summary.csv"
    report_path = output_dir / "phase_ensemble_report.md"
    write_csv(summary_path, rows)
    write_report(report_path, rows)
    pass_count = sum(1 for row in rows if row["status"] == "PASS")

    print(f"phase_ensemble_summary_csv,{summary_path}")
    print(f"phase_ensemble_report,{report_path}")
    print(f"phase_ensemble_pass_count,{pass_count}/{len(rows)}")
    print("phase_ensemble_status,PASS" if pass_count == len(rows) else "phase_ensemble_status,FAIL")

    if args.fail_on_any_fail and pass_count != len(rows):
        raise SystemExit("At least one phase-seed case failed.")


if __name__ == "__main__":
    main()
