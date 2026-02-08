#!/usr/bin/env python3
"""Run a Harris 4D ablation campaign and generate a compact report bundle."""

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt


@dataclass
class Scenario:
    name: str
    description: str
    full_args: list[str]


SCENARIOS = [
    Scenario("baseline_full4d", "Tuned full-channel baseline", []),
    Scenario(
        "ablate_mixed_seed",
        "Disable mixed-sector seed (aw_mode1_amp, piw_mode1_amp)",
        [
            "problem/harris_4d/aw_mode1_amp=0.0",
            "problem/harris_4d/piw_mode1_amp=0.0",
        ],
    ),
    Scenario(
        "ablate_timelike_em",
        "Disable time-like EM mixed coupling",
        ["modes4d/em_source_timelike_gain=0.0"],
    ),
    Scenario(
        "ablate_transverse_force",
        "Disable plasma transverse momw forcing channels",
        [
            "modes4d/plasma_momw_source_gain=0.0",
            "modes4d/plasma_momw_pressure_source_gain=0.0",
        ],
    ),
    Scenario(
        "ablate_em_current_coupling",
        "Disable EM current source coupling",
        ["modes4d/em_source_current_gain=0.0"],
    ),
]


def parse_scan_csv(stdout_text):
    lines = [line.strip() for line in stdout_text.splitlines() if line.strip()]
    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("case,"):
            header_idx = idx
            break
    if header_idx is None:
        raise RuntimeError("Could not find CSV header in scan output.")

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
        raise RuntimeError("Missing controlled/full case rows.")
    return parsed


def parse_float(row, key):
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return math.nan


def fmt(value):
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def run_scenario(
    scenario,
    scan_script,
    binary,
    workdir,
    controlled_input,
    full_input,
    output_dir,
    tlim,
    nlim,
    output_dt,
    athena_args,
    scan_args,
    check_transport_closure,
    check_em_bulk_ledger,
):
    scenario_dir = output_dir / scenario.name
    scenario_dir.mkdir(parents=True, exist_ok=True)

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
        str(scenario_dir / "scan_outputs"),
        "--athena-arg",
        f"parthenon/time/tlim={tlim}",
        "--athena-arg",
        f"parthenon/time/nlim={nlim}",
        "--athena-arg",
        "parthenon/output0/dt=-1",
        "--athena-arg",
        f"parthenon/output1/dt={output_dt}",
    ]
    for scan_arg in scan_args:
        cmd.append(scan_arg)
    for athena_arg in athena_args:
        cmd.extend(["--athena-arg", athena_arg])
    if check_transport_closure:
        cmd.append("--check-transport-closure")
    if check_em_bulk_ledger:
        cmd.append("--check-em-bulk-ledger")
    for full_arg in scenario.full_args:
        cmd.extend(["--full-arg", full_arg])

    proc = subprocess.run(cmd, cwd=workdir, text=True, capture_output=True)
    (scenario_dir / "scan_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (scenario_dir / "scan_stderr.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            f"Scenario {scenario.name} failed with returncode={proc.returncode}. "
            f"See {scenario_dir / 'scan_stderr.log'}"
        )
    return parse_scan_csv(proc.stdout)


def safe_ratio(value, baseline):
    if not math.isfinite(value) or not math.isfinite(baseline) or baseline == 0.0:
        return math.nan
    return value / baseline


def write_csv(summary_csv, rows):
    if not rows:
        raise RuntimeError("No rows to write.")
    fieldnames = list(rows[0].keys())
    with summary_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(report_path, controlled_ref, rows, ablation_checks, overall_status):
    lines = []
    lines.append("# Harris 4D Ablation Report")
    lines.append("")
    lines.append(f"- Generated (UTC): {datetime.now(timezone.utc).isoformat()}")
    lines.append(f"- Overall status: **{overall_status}**")
    lines.append("")
    lines.append("## Controlled Reference")
    lines.append("")
    lines.append(
        "| psi0_span | abs_jw_ew | S_leak_abs | mixed_ew2 | mixed_c2 | closure | transport |"
    )
    lines.append("| ---: | ---: | ---: | ---: | ---: | :---: | :---: |")
    lines.append(
        "| "
        + " | ".join(
            [
                f"`{fmt(controlled_ref['full_final_psi0_span'])}`",
                f"`{fmt(abs(controlled_ref['full_final_jw_ew']))}`",
                f"`{fmt(controlled_ref['full_final_s_leak_abs'])}`",
                f"`{fmt(controlled_ref['full_final_mixed_ew2'])}`",
                f"`{fmt(controlled_ref['full_final_mixed_c2'])}`",
                f"`{controlled_ref['full_closure_status']}`",
                f"`{controlled_ref['full_transport_closure_status']}`",
            ]
        )
        + " |"
    )
    lines.append("")
    lines.append("## Full-Case Ablation Summary")
    lines.append("")
    lines.append(
        "| scenario | description | psi0_span | abs_jw_ew | S_leak_abs | mixed_ew2 | mixed_c2 | "
        "ratio_abs_jw_ew | ratio_S_leak | ratio_mixed_ew2 | ratio_mixed_c2 | closure | transport | "
        "correlation | activity |"
    )
    lines.append(
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
        ":---: | :---: | :---: | :---: |"
    )
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{row['scenario']}`",
                    row["description"],
                    f"`{fmt(row['full_final_psi0_span'])}`",
                    f"`{fmt(abs(row['full_final_jw_ew']))}`",
                    f"`{fmt(row['full_final_s_leak_abs'])}`",
                    f"`{fmt(row['full_final_mixed_ew2'])}`",
                    f"`{fmt(row['full_final_mixed_c2'])}`",
                    f"`{fmt(row['ratio_abs_final_jw_ew'])}`",
                    f"`{fmt(row['ratio_final_s_leak_abs'])}`",
                    f"`{fmt(row['ratio_final_mixed_ew2'])}`",
                    f"`{fmt(row['ratio_final_mixed_c2'])}`",
                    f"`{row['full_closure_status']}`",
                    f"`{row['full_transport_closure_status']}`",
                    f"`{row['full_correlation_status']}`",
                    f"`{row['full_activity_status']}`",
                ]
            )
            + " |"
        )
    lines.append("")
    lines.append("## Causality Checks")
    lines.append("")
    lines.append("| scenario | expected reduction checks | status | details |")
    lines.append("| :--- | :--- | :---: | :--- |")
    for item in ablation_checks:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{item['scenario']}`",
                    item["checks"],
                    f"`{item['status']}`",
                    item["details"],
                ]
            )
            + " |"
        )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_plots(plot_dir, rows):
    plot_dir.mkdir(parents=True, exist_ok=True)
    names = [row["scenario"] for row in rows]
    x = list(range(len(names)))

    abs_jw_ew = [abs(row["full_final_jw_ew"]) for row in rows]
    s_leak_abs = [row["full_final_s_leak_abs"] for row in rows]
    mixed_ew2 = [row["full_final_mixed_ew2"] for row in rows]
    mixed_c2 = [row["full_final_mixed_c2"] for row in rows]

    fig, axs = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)
    axs[0].plot(x, abs_jw_ew, marker="o", label="|JwEw|")
    axs[0].plot(x, s_leak_abs, marker="s", label="S_leak_abs")
    axs[0].plot(x, mixed_ew2, marker="^", label="mixed_ew2")
    axs[0].plot(x, mixed_c2, marker="d", label="mixed_c2")
    axs[0].set_yscale("log")
    axs[0].set_xticks(x, names, rotation=20, ha="right")
    axs[0].set_title("Ablation Channel Magnitudes")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend()

    ratio_jw_ew = [row["ratio_abs_final_jw_ew"] for row in rows]
    ratio_s_leak = [row["ratio_final_s_leak_abs"] for row in rows]
    ratio_mixed_ew2 = [row["ratio_final_mixed_ew2"] for row in rows]
    ratio_mixed_c2 = [row["ratio_final_mixed_c2"] for row in rows]
    axs[1].plot(x, ratio_jw_ew, marker="o", label="ratio |JwEw|")
    axs[1].plot(x, ratio_s_leak, marker="s", label="ratio S_leak_abs")
    axs[1].plot(x, ratio_mixed_ew2, marker="^", label="ratio mixed_ew2")
    axs[1].plot(x, ratio_mixed_c2, marker="d", label="ratio mixed_c2")
    axs[1].axhline(1.0, color="k", linewidth=1.0, alpha=0.5)
    axs[1].axhline(0.5, color="r", linewidth=1.0, alpha=0.4)
    axs[1].set_xticks(x, names, rotation=20, ha="right")
    axs[1].set_title("Ablation Ratios vs Baseline")
    axs[1].set_ylim(bottom=0.0)
    axs[1].grid(True, alpha=0.3)
    axs[1].legend()

    path = plot_dir / "ablation_channels.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument("--controlled-input", default="inputs/harris_4d_controlled.in")
    parser.add_argument("--full-input", default="inputs/harris_4d_full.in")
    parser.add_argument("--output-dir", default="harris_ablation_outputs")
    parser.add_argument("--tlim", type=float, default=0.1)
    parser.add_argument("--nlim", type=int, default=1000)
    parser.add_argument("--output-dt", type=float, default=0.01)
    parser.add_argument(
        "--athena-arg",
        action="append",
        default=[],
        help="Additional athenaPK runtime override passed to both controlled/full runs",
    )
    parser.add_argument(
        "--scan-arg",
        action="append",
        default=[],
        help="Additional harris_scan_matrix.py argument forwarded as-is",
    )
    parser.add_argument(
        "--ablation-strong-ratio-threshold",
        type=float,
        default=0.5,
        help="Expected strong-reduction threshold for selected channels",
    )
    parser.add_argument(
        "--check-transport-closure",
        action="store_true",
        help="Forward transport closure checks into scan runs",
    )
    parser.add_argument(
        "--check-em-bulk-ledger",
        action="store_true",
        help="Forward EM bulk-ledger checks into scan runs",
    )
    parser.add_argument(
        "--fail-on-ablation-check",
        action="store_true",
        help="Exit nonzero if causal ablation checks fail",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    scan_script = repo_root / "scripts" / "harris_scan_matrix.py"
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    output_dir = Path(args.output_dir)
    output_dir = output_dir.resolve() if output_dir.is_absolute() else (workdir / output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    def resolve_input(path_str):
        p = Path(path_str)
        if p.is_absolute():
            return p.resolve()
        for base in (Path.cwd(), workdir, repo_root):
            candidate = (base / p).resolve()
            if candidate.exists():
                return candidate
        return (repo_root / p).resolve()

    controlled_input = resolve_input(args.controlled_input)
    full_input = resolve_input(args.full_input)

    scenario_results = []
    controlled_reference = None
    for scenario in SCENARIOS:
        parsed = run_scenario(
            scenario,
            scan_script,
            binary,
            workdir,
            controlled_input,
            full_input,
            output_dir,
            args.tlim,
            args.nlim,
            args.output_dt,
            args.athena_arg,
            args.scan_arg,
            args.check_transport_closure,
            args.check_em_bulk_ledger,
        )
        controlled_row = parsed["controlled"]
        full_row = parsed["full"]
        if controlled_reference is None:
            controlled_reference = {
                "full_final_psi0_span": parse_float(controlled_row, "final_psi0_span"),
                "full_final_jw_ew": parse_float(controlled_row, "final_jw_ew"),
                "full_final_s_leak_abs": parse_float(controlled_row, "final_s_leak_abs"),
                "full_final_mixed_ew2": parse_float(controlled_row, "final_mixed_ew2"),
                "full_final_mixed_c2": parse_float(controlled_row, "final_mixed_c2"),
                "full_closure_status": controlled_row.get("closure_status", "N/A"),
                "full_transport_closure_status": controlled_row.get(
                    "transport_closure_status", "N/A"
                ),
            }

        scenario_results.append(
            {
                "scenario": scenario.name,
                "description": scenario.description,
                "full_final_psi0_span": parse_float(full_row, "final_psi0_span"),
                "full_final_jw_ew": parse_float(full_row, "final_jw_ew"),
                "full_final_s_leak_abs": parse_float(full_row, "final_s_leak_abs"),
                "full_final_mixed_ew2": parse_float(full_row, "final_mixed_ew2"),
                "full_final_mixed_c2": parse_float(full_row, "final_mixed_c2"),
                "full_closure_status": full_row.get("closure_status", "N/A"),
                "full_transport_closure_status": full_row.get(
                    "transport_closure_status", "N/A"
                ),
                "full_correlation_status": full_row.get("correlation_status", "N/A"),
                "full_activity_status": full_row.get("activity_status", "N/A"),
            }
        )

    baseline = next(row for row in scenario_results if row["scenario"] == "baseline_full4d")
    baseline_abs_jw_ew = abs(baseline["full_final_jw_ew"])
    baseline_s_leak_abs = baseline["full_final_s_leak_abs"]
    baseline_mixed_ew2 = baseline["full_final_mixed_ew2"]
    baseline_mixed_c2 = baseline["full_final_mixed_c2"]

    for row in scenario_results:
        row["ratio_abs_final_jw_ew"] = safe_ratio(abs(row["full_final_jw_ew"]), baseline_abs_jw_ew)
        row["ratio_final_s_leak_abs"] = safe_ratio(row["full_final_s_leak_abs"], baseline_s_leak_abs)
        row["ratio_final_mixed_ew2"] = safe_ratio(row["full_final_mixed_ew2"], baseline_mixed_ew2)
        row["ratio_final_mixed_c2"] = safe_ratio(row["full_final_mixed_c2"], baseline_mixed_c2)

    strict_check_specs = {
        "ablate_mixed_seed": ["ratio_final_mixed_ew2", "ratio_final_mixed_c2"],
        "ablate_transverse_force": ["ratio_abs_final_jw_ew", "ratio_final_s_leak_abs"],
    }

    ablation_checks = []
    checks_ok = True
    threshold = args.ablation_strong_ratio_threshold
    for row in scenario_results:
        if row["scenario"] == "baseline_full4d":
            continue
        keys = strict_check_specs.get(row["scenario"], [])
        if not keys:
            ablation_checks.append(
                {
                    "scenario": row["scenario"],
                    "checks": "none (informational)",
                    "status": "INFO",
                    "details": "not part of strict causal gate",
                }
            )
            continue
        failures = []
        for key in keys:
            value = row.get(key, math.nan)
            if not math.isfinite(value) or value > threshold:
                failures.append(f"{key}={fmt(value)} > {threshold:.2f}")
        status = "PASS" if not failures else "FAIL"
        checks_ok = checks_ok and (status == "PASS")
        ablation_checks.append(
            {
                "scenario": row["scenario"],
                "checks": ", ".join(keys),
                "status": status,
                "details": "none" if not failures else "; ".join(failures),
            }
        )

    overall_status = "PASS" if checks_ok else "FAIL"

    summary_csv = output_dir / "ablation_summary.csv"
    write_csv(summary_csv, scenario_results)
    plot_dir = output_dir / "plots"
    plot_path = make_plots(plot_dir, scenario_results)
    report_path = output_dir / "ablation_report.md"
    write_report(report_path, controlled_reference, scenario_results, ablation_checks, overall_status)

    print(f"summary_csv,{summary_csv}")
    print(f"plot,{plot_path}")
    print(f"report,{report_path}")
    print(f"ablation_status,{overall_status}")
    if args.fail_on_ablation_check and overall_status == "FAIL":
        raise SystemExit("ablation checks failed")


if __name__ == "__main__":
    main()
