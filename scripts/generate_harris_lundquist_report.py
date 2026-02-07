#!/usr/bin/env python3
"""Generate a concise markdown report from a Harris Lundquist summary CSV."""

import argparse
import csv
import math
from datetime import datetime, timezone
from pathlib import Path


def parse_float(row, key):
    raw = row.get(key, "")
    if raw is None or raw == "":
        return math.nan
    try:
        return float(raw)
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


def fmt_val(value):
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def gate_abs(value, threshold):
    if not math.isfinite(value):
        return "FAIL"
    return "PASS" if abs(value) >= threshold else "FAIL"


def gate_max(value, threshold):
    if not math.isfinite(value):
        return "FAIL"
    return "PASS" if value <= threshold else "FAIL"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-csv", required=True, help="Path to lundquist_scan_summary.csv")
    parser.add_argument(
        "--plots-dir",
        default="",
        help="Optional directory containing generated plots/status CSV",
    )
    parser.add_argument(
        "--report-path",
        default="",
        help="Output markdown report path (default: <summary parent>/production_report.md)",
    )
    args = parser.parse_args()

    summary_csv = Path(args.summary_csv).resolve()
    if not summary_csv.exists():
        raise FileNotFoundError(f"Summary CSV not found: {summary_csv}")

    rows = list(csv.DictReader(summary_csv.open("r", encoding="utf-8")))
    if not rows:
        raise RuntimeError("Summary CSV is empty.")

    rows.sort(key=lambda r: parse_float(r, "S"))
    s_vals = [parse_float(r, "S") for r in rows]
    log_s = [math.log10(s) for s in s_vals]

    metrics = {
        "corr_logS_full_final_psi0_span": pearson(
            log_s, [parse_float(r, "full_final_psi0_span") for r in rows]
        ),
        "corr_logS_full_final_s_leak_abs": pearson(
            log_s, [parse_float(r, "full_final_s_leak_abs") for r in rows]
        ),
        "corr_logS_full_final_jw_ew_abs": pearson(
            log_s, [abs(parse_float(r, "full_final_jw_ew")) for r in rows]
        ),
        "corr_logS_full_final_mixed_ew2": pearson(
            log_s, [parse_float(r, "full_final_mixed_ew2") for r in rows]
        ),
        "corr_logS_full_final_em_leak_w_abs": pearson(
            log_s, [abs(parse_float(r, "full_final_em_leak_w")) for r in rows]
        ),
        "corr_logS_full_final_helicity_sub_abs": pearson(
            log_s, [abs(parse_float(r, "full_final_helicity_sub")) for r in rows]
        ),
        "corr_logS_full_final_edotb_sub_abs": pearson(
            log_s, [abs(parse_float(r, "full_final_edotb_sub")) for r in rows]
        ),
        "corr_logS_full_final_jw_mode_activity_proxy": pearson(
            log_s,
            [
                max(
                    parse_float(r, "full_final_jw_mode_l2_1"),
                    abs(parse_float(r, "full_final_jw_ew")),
                )
                for r in rows
            ],
        ),
    }

    gate_specs = [
        ("corr_logS_full_final_psi0_span", "abs_min", 7.0e-1),
        ("corr_logS_full_final_s_leak_abs", "abs_min", 5.0e-1),
        ("corr_logS_full_final_jw_ew_abs", "abs_min", 5.0e-1),
        ("corr_logS_full_final_mixed_ew2", "abs_min", 7.0e-1),
        ("corr_logS_full_final_em_leak_w_abs", "abs_min", 2.0e-1),
        ("corr_logS_full_final_helicity_sub_abs", "abs_min", 5.0e-1),
        ("corr_logS_full_final_edotb_sub_abs", "abs_min", 5.0e-1),
        ("corr_logS_full_final_jw_mode_activity_proxy", "abs_min", 5.0e-1),
        ("corr_logS_full_final_helicity_sub_abs", "max", -2.0e-1),
        ("corr_logS_full_final_edotb_sub_abs", "max", -2.0e-1),
    ]

    gate_rows = []
    gate_ok = True
    for key, mode, threshold in gate_specs:
        value = metrics[key]
        if mode == "abs_min":
            status = gate_abs(value, threshold)
            condition = f"|{key}| >= {threshold:.1e}"
        else:
            status = gate_max(value, threshold)
            condition = f"{key} <= {threshold:.1e}"
        gate_ok = gate_ok and (status == "PASS")
        gate_rows.append((condition, value, status))

    per_case_ok = True
    for row in rows:
        for status_key in (
            "full_closure_status",
            "full_transport_closure_status",
            "full_correlation_status",
            "full_activity_status",
        ):
            if row.get(status_key, "FAIL") != "PASS":
                per_case_ok = False
                break

    overall = "PASS" if (per_case_ok and gate_ok) else "FAIL"

    report_path = (
        Path(args.report_path).resolve()
        if args.report_path
        else (summary_csv.parent / "production_report.md").resolve()
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)

    plots_dir = Path(args.plots_dir).resolve() if args.plots_dir else summary_csv.parent / "plots"
    primary_plot = plots_dir / "lundquist_primary_channels.png"
    subscale_plot = plots_dir / "lundquist_subscale_channels.png"
    status_csv = plots_dir / "lundquist_status.csv"

    lines = []
    lines.append("# Harris Lundquist Production Report")
    lines.append("")
    lines.append(f"- Generated (UTC): {datetime.now(timezone.utc).isoformat()}")
    lines.append(f"- Overall status: **{overall}**")
    lines.append(f"- Summary CSV: `{summary_csv}`")
    lines.append(f"- Plot (primary): `{primary_plot}`")
    lines.append(f"- Plot (subscale): `{subscale_plot}`")
    lines.append(f"- Status CSV: `{status_csv}`")
    lines.append("")
    lines.append("## Scan-Level Correlation Gates")
    lines.append("")
    lines.append("| Gate | Value | Status |")
    lines.append("| --- | ---: | :---: |")
    for condition, value, status in gate_rows:
        lines.append(f"| `{condition}` | `{fmt_val(value)}` | `{status}` |")
    lines.append("")
    lines.append("## Per-S Summary")
    lines.append("")
    lines.append(
        "| S | psi0_span | |JwEw| | S_leak_abs | mixed_ew2 | mixed_c2 | "
        "jw_mode_activity_proxy | closure | transport | correlation | activity |"
    )
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: | :---: | :---: |")
    for row in rows:
        jw_mode_activity_proxy = max(
            parse_float(row, "full_final_jw_mode_l2_1"),
            abs(parse_float(row, "full_final_jw_ew")),
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{parse_float(row, 'S'):.0f}`",
                    f"`{fmt_val(parse_float(row, 'full_final_psi0_span'))}`",
                    f"`{fmt_val(abs(parse_float(row, 'full_final_jw_ew')))}`",
                    f"`{fmt_val(parse_float(row, 'full_final_s_leak_abs'))}`",
                    f"`{fmt_val(parse_float(row, 'full_final_mixed_ew2'))}`",
                    f"`{fmt_val(parse_float(row, 'full_final_mixed_c2'))}`",
                    f"`{fmt_val(jw_mode_activity_proxy)}`",
                    f"`{row.get('full_closure_status', 'N/A')}`",
                    f"`{row.get('full_transport_closure_status', 'N/A')}`",
                    f"`{row.get('full_correlation_status', 'N/A')}`",
                    f"`{row.get('full_activity_status', 'N/A')}`",
                ]
            )
            + " |"
        )

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"production_report,{report_path}")
    print(f"overall_status,{overall}")


if __name__ == "__main__":
    main()
