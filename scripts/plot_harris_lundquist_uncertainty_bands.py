#!/usr/bin/env python3
"""Build uncertainty-band tables/plots from a phase-ensemble summary."""

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def parse_float(row, key):
    raw = row.get(key, "")
    if raw is None or raw == "":
        return math.nan
    try:
        return float(raw)
    except ValueError:
        return math.nan


def load_csv_rows(path):
    rows = list(csv.DictReader(path.open("r", encoding="utf-8")))
    if not rows:
        raise RuntimeError(f"CSV has no rows: {path}")
    return rows


def stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return {"count": 0, "mean": math.nan, "std": math.nan, "min": math.nan, "max": math.nan}
    return {
        "count": len(finite),
        "mean": statistics.fmean(finite),
        "std": statistics.pstdev(finite) if len(finite) > 1 else 0.0,
        "min": min(finite),
        "max": max(finite),
    }


def fmt(v):
    return "nan" if not math.isfinite(v) else f"{v:.6e}"


def write_stats_csv(path, stats_rows):
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["metric", "S", "count", "mean", "std", "min", "max"])
        for row in stats_rows:
            writer.writerow(row)


def write_seed_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(
            [
                "seed",
                "status",
                "phase_x",
                "phase_z",
                "min_ratio_psi_w0_span",
                "min_ratio_max_abs_dpsi_w0_dt",
                "slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt",
                "min_full_jw_ew_abs",
                "summary_csv",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.get("seed", ""),
                    row.get("status", ""),
                    row.get("phase_x", ""),
                    row.get("phase_z", ""),
                    row.get("min_ratio_psi_w0_span", ""),
                    row.get("min_ratio_max_abs_dpsi_w0_dt", ""),
                    row.get("slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt", ""),
                    row.get("min_full_jw_ew_abs", ""),
                    row.get("summary_csv", ""),
                ]
            )


def write_report(path, total_rows, pass_rows, metric_stats, scalar_stats):
    lines = []
    lines.append("# Harris Lundquist Uncertainty Bands")
    lines.append("")
    lines.append(f"- Total seeds: **{len(total_rows)}**")
    lines.append(f"- Passing seeds used: **{len(pass_rows)}**")
    lines.append("")
    lines.append("## Scalar Uncertainty (Across Passing Seeds)")
    lines.append("")
    for name, values in scalar_stats.items():
        s = stats(values)
        lines.append(
            f"- `{name}`: mean=`{fmt(s['mean'])}`, std=`{fmt(s['std'])}`, "
            f"min=`{fmt(s['min'])}`, max=`{fmt(s['max'])}`"
        )
    lines.append("")
    lines.append("## Per-S Bands")
    lines.append("")
    lines.append("| metric | S | count | mean | std | min | max |")
    lines.append("| :--- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for metric_name in sorted(metric_stats.keys()):
        for s in sorted(metric_stats[metric_name].keys()):
            st = metric_stats[metric_name][s]
            lines.append(
                "| "
                + " | ".join(
                    [
                        f"`{metric_name}`",
                        f"`{s:.6e}`",
                        str(st["count"]),
                        f"`{fmt(st['mean'])}`",
                        f"`{fmt(st['std'])}`",
                        f"`{fmt(st['min'])}`",
                        f"`{fmt(st['max'])}`",
                    ]
                )
                + " |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_metric_band(ax, s_vals, means, stds, title, ylabel, log_y=False):
    if not s_vals:
        ax.set_title(title)
        ax.set_xlabel("Lundquist S")
        ax.set_ylabel(ylabel)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return

    ax.plot(s_vals, means, marker="o", linewidth=1.8)
    lower = [m - d for m, d in zip(means, stds)]
    upper = [m + d for m, d in zip(means, stds)]
    if log_y:
        positive_floor = min([m for m in means if m > 0.0], default=1.0e-30) * 1.0e-3
        lower = [max(v, positive_floor) for v in lower]
    ax.fill_between(s_vals, lower, upper, alpha=0.25)
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel("Lundquist S")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-ensemble-summary-csv", required=True)
    parser.add_argument(
        "--output-dir",
        default="",
        help="Output directory (default: <phase-ensemble-summary parent>/uncertainty_bands)",
    )
    args = parser.parse_args()

    phase_summary = Path(args.phase_ensemble_summary_csv).resolve()
    if not phase_summary.exists():
        raise FileNotFoundError(f"Missing phase ensemble summary CSV: {phase_summary}")

    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else (phase_summary.parent / "uncertainty_bands").resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = load_csv_rows(phase_summary)
    pass_rows = [row for row in rows if (row.get("status", "").upper() == "PASS")]
    if not pass_rows:
        raise RuntimeError("No passing rows in phase ensemble summary.")

    metric_extractors = {
        "ratio_psi_w0_span": lambda row: parse_float(row, "ratio_full_over_controlled_psi_w0_span"),
        "ratio_max_abs_dpsi_w0_dt": lambda row: parse_float(
            row, "ratio_full_over_controlled_max_abs_dpsi_w0_dt"
        ),
        "full_jw_ew_abs": lambda row: abs(parse_float(row, "full_final_jw_ew")),
        "full_s_leak_abs": lambda row: parse_float(row, "full_final_s_leak_abs"),
        "full_mixed_ew2": lambda row: parse_float(row, "full_final_mixed_ew2"),
    }

    per_metric = {name: defaultdict(list) for name in metric_extractors}

    for phase_row in pass_rows:
        summary_csv_raw = (phase_row.get("summary_csv") or "").strip()
        if not summary_csv_raw:
            continue
        summary_csv = Path(summary_csv_raw)
        if not summary_csv.is_absolute():
            summary_csv = (phase_summary.parent / summary_csv).resolve()
        if not summary_csv.exists():
            raise FileNotFoundError(f"Missing per-seed summary CSV: {summary_csv}")
        for row in load_csv_rows(summary_csv):
            s = parse_float(row, "S")
            if not math.isfinite(s):
                continue
            for metric_name, metric_fn in metric_extractors.items():
                value = metric_fn(row)
                if math.isfinite(value):
                    per_metric[metric_name][s].append(value)

    metric_stats = {}
    stats_rows = []
    for metric_name, by_s in per_metric.items():
        metric_stats[metric_name] = {}
        for s in sorted(by_s.keys()):
            st = stats(by_s[s])
            metric_stats[metric_name][s] = st
            stats_rows.append(
                [metric_name, f"{s:.12g}", st["count"], fmt(st["mean"]), fmt(st["std"]), fmt(st["min"]), fmt(st["max"])]
            )

    scalar_stats = {
        "min_ratio_psi_w0_span": [parse_float(r, "min_ratio_psi_w0_span") for r in pass_rows],
        "min_ratio_max_abs_dpsi_w0_dt": [
            parse_float(r, "min_ratio_max_abs_dpsi_w0_dt") for r in pass_rows
        ],
        "slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt": [
            parse_float(r, "slope_delta_full_minus_controlled_max_abs_dpsi_w0_dt") for r in pass_rows
        ],
        "min_full_jw_ew_abs": [parse_float(r, "min_full_jw_ew_abs") for r in pass_rows],
    }

    figure_path = output_dir / "uncertainty_bands.png"
    stats_csv_path = output_dir / "uncertainty_bands_stats.csv"
    seeds_csv_path = output_dir / "uncertainty_bands_seed_scalars.csv"
    report_path = output_dir / "uncertainty_bands_report.md"

    fig, axs = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    panels = [
        ("ratio_psi_w0_span", "full/control psi_w0 span", False),
        ("ratio_max_abs_dpsi_w0_dt", "full/control max|dpsi_w0/dt|", False),
        ("full_jw_ew_abs", "|JwEw|", True),
        ("full_s_leak_abs", "|S_leak|", True),
    ]

    for ax, (metric_name, label, log_y) in zip(axs.flat, panels):
        by_s = metric_stats.get(metric_name, {})
        s_vals = sorted(by_s.keys())
        means = [by_s[s]["mean"] for s in s_vals]
        stds = [by_s[s]["std"] for s in s_vals]
        plot_metric_band(
            ax,
            s_vals,
            means,
            stds,
            f"{label} uncertainty band",
            label,
            log_y=log_y,
        )

    fig.savefig(figure_path, dpi=150)
    plt.close(fig)

    write_stats_csv(stats_csv_path, stats_rows)
    write_seed_csv(seeds_csv_path, pass_rows)
    write_report(report_path, rows, pass_rows, metric_stats, scalar_stats)

    print(f"uncertainty_bands_plot,{figure_path}")
    print(f"uncertainty_bands_stats_csv,{stats_csv_path}")
    print(f"uncertainty_bands_seed_csv,{seeds_csv_path}")
    print(f"uncertainty_bands_report,{report_path}")


if __name__ == "__main__":
    main()
