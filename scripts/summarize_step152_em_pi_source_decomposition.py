#!/usr/bin/env python3
"""Decompose EM pi mode-0 source integrals against existing mode-0 RHS diagnostics."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


CHANNEL_CONFIG = {
    "pi0": {
        "src_col": "m4d_int_srcpi0_mode0",
        "rhs_total_col": "m4d_int_rhs_a0_mode0_total",
        "components": {
            "lap": "m4d_int_rhs_a0_mode0_lap",
            "mass": "m4d_int_rhs_a0_mode0_mass",
            "current": "m4d_int_rhs_a0_mode0_current",
            "damping": "m4d_int_rhs_a0_mode0_damping",
            "timelike": "m4d_int_rhs_a0_mode0_timelike",
            "gauge": "m4d_int_rhs_a0_mode0_gauge",
        },
    },
    "piy": {
        "src_col": "m4d_int_srcpiy_mode0",
        "rhs_total_col": "m4d_int_rhs_ay_mode0_total",
        "components": {
            "lap": "m4d_int_rhs_ay_mode0_lap",
            "mass": "m4d_int_rhs_ay_mode0_mass",
            "current": "m4d_int_rhs_ay_mode0_current",
            "damping": "m4d_int_rhs_ay_mode0_damping",
            "spatial_mixed": "m4d_int_rhs_ay_mode0_spatial_mixed",
            "even_bridge": "m4d_int_rhs_ay_mode0_even_bridge",
            "geometry_shift": "m4d_int_rhs_ay_mode0_geometry_shift",
            "w0_response": "m4d_int_rhs_ay_mode0_w0_response",
        },
    },
    "piw": {
        "src_col": "m4d_int_srcpiw_mode0",
        "rhs_total_col": "m4d_int_rhs_aw_mode0_total",
        "components": {
            "lap": "m4d_int_rhs_aw_mode0_lap",
            "current": "m4d_int_rhs_aw_mode0_current",
            "damping": "m4d_int_rhs_aw_mode0_damping",
            "spatial_mixed": "m4d_int_rhs_aw_mode0_spatial_mixed",
            "timelike": "m4d_int_rhs_aw_mode0_timelike",
        },
    },
}


def load_harris_scan_matrix(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def derivative_series(times, values):
    if values is None or len(times) < 3 or len(values) != len(times):
        return [], []
    out_t = []
    out_v = []
    for i in range(2, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        out_t.append(times[i])
        out_v.append((values[i] - values[i - 1]) / dt)
    return out_t, out_v


def max_abs(values):
    if not values:
        return math.nan
    return max(abs(v) for v in values)


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


def parse_case_rows(run_log: Path) -> tuple[dict[str, str], dict[str, str]]:
    lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()
    header = None
    controlled = None
    full = None
    for i, line in enumerate(lines):
        if not line.startswith("case,final_psi0_span,"):
            continue
        header = line
        for j in range(i + 1, min(i + 24, len(lines))):
            if lines[j].startswith("controlled,"):
                controlled = lines[j]
            elif lines[j].startswith("full,"):
                full = lines[j]
        if header and controlled and full:
            break
    if header is None or controlled is None or full is None:
        raise RuntimeError(f"Could not find controlled/full rows in {run_log}")
    keys = next(csv.reader([header]))
    controlled_vals = next(csv.reader([controlled]))
    full_vals = next(csv.reader([full]))
    return dict(zip(keys, controlled_vals)), dict(zip(keys, full_vals))


def dominant_component(component_max):
    ranked = [(name, value) for name, value in component_max.items() if math.isfinite(value)]
    if not ranked:
        return "nan"
    return max(ranked, key=lambda item: item[1])[0]


def safe_ratio(baseline: float, value: float) -> float:
    if not math.isfinite(baseline) or not math.isfinite(value) or value == 0.0:
        return math.nan
    return baseline / value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Case glob")
    parser.add_argument("--baseline-case", default="baseline_pi_transport")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows = []
    summaries = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        hst = case_dir / "harris_full.out1.hst"
        if not run_log.exists() or not hst.exists():
            continue
        _, full = parse_case_rows(run_log)
        cols = hsm.parse_hst(hst)
        times = cols.get("time", [])
        case_rows = []
        for channel, config in CHANNEL_CONFIG.items():
            _, src_rate = derivative_series(times, cols.get(config["src_col"]))
            _, rhs_total_rate = derivative_series(times, cols.get(config["rhs_total_col"]))
            if not src_rate or len(src_rate) != len(rhs_total_rate):
                continue
            component_rates = {}
            component_max = {}
            for name, col in config["components"].items():
                _, series = derivative_series(times, cols.get(col))
                if series and len(series) == len(src_rate):
                    component_rates[name] = series
                    component_max[name] = max_abs(series)
                else:
                    component_max[name] = math.nan
            gap_rate = [src - rhs for src, rhs in zip(src_rate, rhs_total_rate)]
            max_src = max_abs(src_rate)
            max_rhs_total = max_abs(rhs_total_rate)
            max_gap = max_abs(gap_rate)
            row = {
                "case_name": case_dir.name,
                "channel": channel,
                "transport_status": full.get(f"{channel}_transport_status", ""),
                "ledger_status": full.get("em_bulk_ledger_status", ""),
                "max_abs_src_rate": max_src,
                "max_abs_rhs_total_rate": max_rhs_total,
                "max_abs_gap_rate": max_gap,
                "corr_src_vs_rhs_total": pearson(src_rate, rhs_total_rate),
                "corr_src_vs_gap": pearson(src_rate, gap_rate),
                "dominant_rhs_component": dominant_component(component_max),
            }
            for name, value in component_max.items():
                row[f"max_abs_{name}_rate"] = value
            rows.append(row)
            case_rows.append(row)

        if case_rows:
            top = max(case_rows, key=lambda row: row["max_abs_src_rate"])
            summaries.append(
                {
                    "case_name": case_dir.name,
                    "top_channel": top["channel"],
                    "top_max_abs_src_rate": top["max_abs_src_rate"],
                    "top_max_abs_rhs_total_rate": top["max_abs_rhs_total_rate"],
                    "top_max_abs_gap_rate": top["max_abs_gap_rate"],
                    "top_dominant_rhs_component": top["dominant_rhs_component"],
                    "top_corr_src_vs_rhs_total": top["corr_src_vs_rhs_total"],
                    "top_corr_src_vs_gap": top["corr_src_vs_gap"],
                }
            )

    baseline_by_channel = {
        row["channel"]: row["max_abs_gap_rate"]
        for row in rows
        if row["case_name"] == args.baseline_case
    }
    for row in rows:
        row["gap_reduction_vs_baseline"] = safe_ratio(
            baseline_by_channel.get(row["channel"], math.nan), row["max_abs_gap_rate"]
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_name",
        "channel",
        "transport_status",
        "ledger_status",
        "max_abs_src_rate",
        "max_abs_rhs_total_rate",
        "max_abs_gap_rate",
        "corr_src_vs_rhs_total",
        "corr_src_vs_gap",
        "dominant_rhs_component",
        "max_abs_lap_rate",
        "max_abs_mass_rate",
        "max_abs_current_rate",
        "max_abs_damping_rate",
        "max_abs_timelike_rate",
        "max_abs_gauge_rate",
        "max_abs_spatial_mixed_rate",
        "max_abs_even_bridge_rate",
        "max_abs_geometry_shift_rate",
        "max_abs_w0_response_rate",
        "gap_reduction_vs_baseline",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary_csv = Path(args.summary_csv)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_fields = [
        "case_name",
        "top_channel",
        "top_max_abs_src_rate",
        "top_max_abs_rhs_total_rate",
        "top_max_abs_gap_rate",
        "top_dominant_rhs_component",
        "top_corr_src_vs_rhs_total",
        "top_corr_src_vs_gap",
    ]
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summaries)

    out_md = Path(args.out_md)
    lines = [
        "# Step-152 EM Pi Source Decomposition",
        "",
        f"Parsed `{len(rows)}` channel rows from `{args.glob}`.",
        "",
        "Definitions:",
        "",
        "- `src_rate = d/dt(m4d_int_srcpi*_mode0)`",
        "- `rhs_total_rate = d/dt(m4d_int_rhs_*_mode0_total)`",
        "- `gap_rate = src_rate - rhs_total_rate`",
        "",
        "The gap isolates source-integral contribution not explained by the signed mode-0 RHS diagnostics.",
        "",
        f"Baseline case: `{args.baseline_case}`.",
        "",
        "| case | channel | transport | ledger | max|src| | max|rhs_total| | max|gap| | gap reduction vs baseline | corr(src,rhs_total) | corr(src,gap) | dominant rhs component |",
        "|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['case_name']} | {row['channel']} | {row['transport_status']} | {row['ledger_status']} | "
            f"{fmt(row['max_abs_src_rate'])} | {fmt(row['max_abs_rhs_total_rate'])} | "
            f"{fmt(row['max_abs_gap_rate'])} | {fmt(row['gap_reduction_vs_baseline'])} | "
            f"{fmt(row['corr_src_vs_rhs_total'])} | {fmt(row['corr_src_vs_gap'])} | "
            f"{row['dominant_rhs_component']} |"
        )

    baseline_rows = [row for row in rows if row["case_name"] == args.baseline_case]
    if baseline_rows:
        lines.extend(
            [
                "",
                "## Baseline Component Ranking",
                "",
                "| channel | dominant rhs component | lap | mass | current | damping | timelike | gauge | spatial_mixed | even_bridge | geometry_shift | w0_response | gap |",
                "|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in baseline_rows:
            lines.append(
                f"| {row['channel']} | {row['dominant_rhs_component']} | "
                f"{fmt(row.get('max_abs_lap_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_mass_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_current_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_damping_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_timelike_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_gauge_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_spatial_mixed_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_even_bridge_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_geometry_shift_rate', math.nan))} | "
                f"{fmt(row.get('max_abs_w0_response_rate', math.nan))} | "
                f"{fmt(row['max_abs_gap_rate'])} |"
            )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{summary_csv}")
    print(f"channels_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
