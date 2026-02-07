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


def mode_activity_proxy(jw_mode_l2, jw_ew):
    candidates = []
    if math.isfinite(jw_mode_l2):
        candidates.append(jw_mode_l2)
    if math.isfinite(jw_ew):
        candidates.append(abs(jw_ew))
    if not candidates:
        return math.nan
    return max(candidates)


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
    parser.add_argument(
        "--psi0-enhancement-min",
        type=float,
        default=1.0,
        help="Minimum min(full_psi0_span/controlled_psi0_span) for PASS",
    )
    parser.add_argument(
        "--max-abs-dpsi0-enhancement-min",
        type=float,
        default=1.0,
        help="Minimum min(full max_abs_dpsi0_dt / controlled max_abs_dpsi0_dt) for PASS",
    )
    parser.add_argument(
        "--controlled-max-jw-ew-abs",
        type=float,
        default=1.0e-12,
        help="Maximum max(|controlled JwEw|) allowed for controlled inactivity PASS",
    )
    parser.add_argument(
        "--controlled-max-s-leak-abs",
        type=float,
        default=1.0e-12,
        help="Maximum max(controlled S_leak_abs) allowed for controlled inactivity PASS",
    )
    parser.add_argument(
        "--controlled-max-jw-mode-activity-proxy",
        type=float,
        default=1.0e-12,
        help="Maximum max(controlled jw_mode_activity_proxy) for controlled inactivity PASS",
    )
    parser.add_argument(
        "--full-min-jw-ew-abs",
        type=float,
        default=1.0e-14,
        help="Minimum min(|full JwEw|) required for full activation PASS",
    )
    parser.add_argument(
        "--full-min-s-leak-abs",
        type=float,
        default=1.0e-13,
        help="Minimum min(full S_leak_abs) required for full activation PASS",
    )
    parser.add_argument(
        "--full-min-jw-mode-activity-proxy",
        type=float,
        default=1.0e-14,
        help="Minimum min(full jw_mode_activity_proxy) required for full activation PASS",
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
        "corr_logS_full_final_psi_proj_span": pearson(
            log_s, [parse_float(r, "full_final_psi_proj_span") for r in rows]
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
                mode_activity_proxy(
                    parse_float(r, "full_final_jw_mode_l2_1"),
                    parse_float(r, "full_final_jw_ew"),
                )
                for r in rows
            ],
        ),
        "corr_logS_full_max_abs_dpsi0_dt": pearson(
            log_s, [parse_float(r, "full_max_abs_dpsi0_dt") for r in rows]
        ),
        "corr_logS_full_max_dpsi0_dt": pearson(
            log_s, [parse_float(r, "full_max_dpsi0_dt") for r in rows]
        ),
        "corr_logS_full_max_abs_dpsi_proj_dt": pearson(
            log_s, [parse_float(r, "full_max_abs_dpsi_proj_dt") for r in rows]
        ),
        "corr_logS_full_max_dpsi_proj_dt": pearson(
            log_s, [parse_float(r, "full_max_dpsi_proj_dt") for r in rows]
        ),
    }
    corr_signal_spans = {
        "corr_logS_full_final_psi0_span": series_span(
            [parse_float(r, "full_final_psi0_span") for r in rows]
        ),
        "corr_logS_full_final_psi_proj_span": series_span(
            [parse_float(r, "full_final_psi_proj_span") for r in rows]
        ),
        "corr_logS_full_final_s_leak_abs": series_span(
            [parse_float(r, "full_final_s_leak_abs") for r in rows]
        ),
        "corr_logS_full_final_jw_ew_abs": series_span(
            [abs(parse_float(r, "full_final_jw_ew")) for r in rows]
        ),
        "corr_logS_full_final_jw_mode_activity_proxy": series_span(
            [
                mode_activity_proxy(
                    parse_float(r, "full_final_jw_mode_l2_1"),
                    parse_float(r, "full_final_jw_ew"),
                )
                for r in rows
            ]
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

    psi_ratio = []
    psi_proj_ratio = []
    rate_ratio = []
    rate_proj_ratio = []
    controlled_jw_ew_abs = []
    controlled_s_leak_abs = []
    controlled_jw_mode_proxy = []
    full_jw_ew_abs = []
    full_s_leak_abs = []
    full_jw_mode_proxy = []
    for row in rows:
        psi_c = parse_float(row, "controlled_final_psi0_span")
        psi_f = parse_float(row, "full_final_psi0_span")
        if math.isfinite(psi_c) and psi_c > 0.0 and math.isfinite(psi_f):
            psi_ratio.append(psi_f / psi_c)
        else:
            psi_ratio.append(math.nan)
        psi_proj_c = parse_float(row, "controlled_final_psi_proj_span")
        psi_proj_f = parse_float(row, "full_final_psi_proj_span")
        if (
            math.isfinite(psi_proj_c)
            and psi_proj_c > 0.0
            and math.isfinite(psi_proj_f)
        ):
            psi_proj_ratio.append(psi_proj_f / psi_proj_c)
        else:
            psi_proj_ratio.append(math.nan)
        rate_c = parse_float(row, "controlled_max_abs_dpsi0_dt")
        rate_f = parse_float(row, "full_max_abs_dpsi0_dt")
        if math.isfinite(rate_c) and rate_c > 0.0 and math.isfinite(rate_f):
            rate_ratio.append(rate_f / rate_c)
        else:
            rate_ratio.append(math.nan)
        rate_proj_c = parse_float(row, "controlled_max_abs_dpsi_proj_dt")
        rate_proj_f = parse_float(row, "full_max_abs_dpsi_proj_dt")
        if (
            math.isfinite(rate_proj_c)
            and rate_proj_c > 0.0
            and math.isfinite(rate_proj_f)
        ):
            rate_proj_ratio.append(rate_proj_f / rate_proj_c)
        else:
            rate_proj_ratio.append(math.nan)

        ctrl_jw = abs(parse_float(row, "controlled_final_jw_ew"))
        full_jw = abs(parse_float(row, "full_final_jw_ew"))
        ctrl_sl = parse_float(row, "controlled_final_s_leak_abs")
        full_sl = parse_float(row, "full_final_s_leak_abs")
        ctrl_proxy = mode_activity_proxy(parse_float(row, "controlled_final_jw_mode_l2_1"), ctrl_jw)
        full_proxy = mode_activity_proxy(parse_float(row, "full_final_jw_mode_l2_1"), full_jw)
        controlled_jw_ew_abs.append(ctrl_jw)
        full_jw_ew_abs.append(full_jw)
        controlled_s_leak_abs.append(ctrl_sl)
        full_s_leak_abs.append(full_sl)
        controlled_jw_mode_proxy.append(ctrl_proxy)
        full_jw_mode_proxy.append(full_proxy)

    compare_gate_rows = []
    compare_specs = [
        (
            "min(full_psi0_span/controlled_psi0_span)",
            series_min(psi_ratio),
            "min",
            args.psi0_enhancement_min,
        ),
        (
            "max(|controlled JwEw|)",
            series_max(controlled_jw_ew_abs),
            "max",
            args.controlled_max_jw_ew_abs,
        ),
        (
            "max(controlled S_leak_abs)",
            series_max(controlled_s_leak_abs),
            "max",
            args.controlled_max_s_leak_abs,
        ),
        (
            "max(controlled jw_mode_activity_proxy)",
            series_max(controlled_jw_mode_proxy),
            "max",
            args.controlled_max_jw_mode_activity_proxy,
        ),
        (
            "min(full_max_abs_dpsi0_dt/controlled_max_abs_dpsi0_dt)",
            series_min(rate_ratio),
            "min",
            args.max_abs_dpsi0_enhancement_min,
        ),
        (
            "min(|full JwEw|)",
            series_min(full_jw_ew_abs),
            "min",
            args.full_min_jw_ew_abs,
        ),
        (
            "min(full S_leak_abs)",
            series_min(full_s_leak_abs),
            "min",
            args.full_min_s_leak_abs,
        ),
        (
            "min(full jw_mode_activity_proxy)",
            series_min(full_jw_mode_proxy),
            "min",
            args.full_min_jw_mode_activity_proxy,
        ),
    ]
    for name, value, mode, threshold in compare_specs:
        if mode == "min":
            status = "FAIL" if (not math.isfinite(value) or value < threshold) else "PASS"
            cond = f"{name} >= {threshold:.1e}"
        elif mode == "max":
            status = gate_max(value, threshold)
            cond = f"{name} <= {threshold:.1e}"
        else:
            status = gate_abs(value, threshold)
            cond = f"|{name}| >= {threshold:.1e}"
        gate_ok = gate_ok and (status == "PASS")
        compare_gate_rows.append((cond, value, status))

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
        "| S | psi0_ctrl | psi0_full | psi0_ratio | "
        "psi_proj_ctrl | psi_proj_full | psi_proj_ratio | "
        "max|dpsi/dt|_ctrl | max|dpsi/dt|_full | rate_ratio | "
        "max|dpsi_proj/dt|_ctrl | max|dpsi_proj/dt|_full | rate_proj_ratio | "
        "|JwEw|_ctrl | |JwEw|_full | JwEw_ratio | "
        "S_leak_ctrl | S_leak_full | S_leak_ratio | "
        "mixed_ew2_ctrl | mixed_ew2_full | mixed_ew2_ratio | "
        "jw_mode_proxy_ctrl | jw_mode_proxy_full | jw_mode_proxy_ratio | "
        "closure | transport | correlation | activity |"
    )
    lines.append(
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
        "---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
        "---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: | :---: | :---: |"
    )
    for row in rows:
        psi_ctrl = parse_float(row, "controlled_final_psi0_span")
        psi_full = parse_float(row, "full_final_psi0_span")
        psi_ratio_val = (
            psi_full / psi_ctrl if (math.isfinite(psi_ctrl) and psi_ctrl > 0.0) else math.nan
        )
        psi_proj_ctrl = parse_float(row, "controlled_final_psi_proj_span")
        psi_proj_full = parse_float(row, "full_final_psi_proj_span")
        psi_proj_ratio_val = (
            psi_proj_full / psi_proj_ctrl
            if (math.isfinite(psi_proj_ctrl) and psi_proj_ctrl > 0.0)
            else math.nan
        )
        jw_ctrl_abs = abs(parse_float(row, "controlled_final_jw_ew"))
        jw_full_abs = abs(parse_float(row, "full_final_jw_ew"))
        jw_ratio = (
            jw_full_abs / jw_ctrl_abs if (math.isfinite(jw_ctrl_abs) and jw_ctrl_abs > 0.0) else math.nan
        )
        s_leak_ctrl = parse_float(row, "controlled_final_s_leak_abs")
        s_leak_full = parse_float(row, "full_final_s_leak_abs")
        s_leak_ratio = (
            s_leak_full / s_leak_ctrl
            if (math.isfinite(s_leak_ctrl) and s_leak_ctrl > 0.0)
            else math.nan
        )
        mixed_ew2_ctrl = parse_float(row, "controlled_final_mixed_ew2")
        mixed_ew2_full = parse_float(row, "full_final_mixed_ew2")
        mixed_ew2_ratio = (
            mixed_ew2_full / mixed_ew2_ctrl
            if (math.isfinite(mixed_ew2_ctrl) and mixed_ew2_ctrl > 0.0)
            else math.nan
        )
        dpsi_rate_ctrl = parse_float(row, "controlled_max_abs_dpsi0_dt")
        dpsi_rate_full = parse_float(row, "full_max_abs_dpsi0_dt")
        dpsi_rate_ratio = (
            dpsi_rate_full / dpsi_rate_ctrl
            if (math.isfinite(dpsi_rate_ctrl) and dpsi_rate_ctrl > 0.0)
            else math.nan
        )
        dpsi_proj_rate_ctrl = parse_float(row, "controlled_max_abs_dpsi_proj_dt")
        dpsi_proj_rate_full = parse_float(row, "full_max_abs_dpsi_proj_dt")
        dpsi_proj_rate_ratio = (
            dpsi_proj_rate_full / dpsi_proj_rate_ctrl
            if (math.isfinite(dpsi_proj_rate_ctrl) and dpsi_proj_rate_ctrl > 0.0)
            else math.nan
        )
        jw_mode_activity_proxy = mode_activity_proxy(
            parse_float(row, "full_final_jw_mode_l2_1"),
            parse_float(row, "full_final_jw_ew"),
        )
        jw_mode_activity_proxy_ctrl = mode_activity_proxy(
            parse_float(row, "controlled_final_jw_mode_l2_1"),
            parse_float(row, "controlled_final_jw_ew"),
        )
        jw_mode_proxy_ratio = (
            jw_mode_activity_proxy / jw_mode_activity_proxy_ctrl
            if (
                math.isfinite(jw_mode_activity_proxy_ctrl)
                and jw_mode_activity_proxy_ctrl > 0.0
            )
            else math.nan
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{parse_float(row, 'S'):.0f}`",
                    f"`{fmt_val(psi_ctrl)}`",
                    f"`{fmt_val(psi_full)}`",
                    f"`{fmt_val(psi_ratio_val)}`",
                    f"`{fmt_val(psi_proj_ctrl)}`",
                    f"`{fmt_val(psi_proj_full)}`",
                    f"`{fmt_val(psi_proj_ratio_val)}`",
                    f"`{fmt_val(dpsi_rate_ctrl)}`",
                    f"`{fmt_val(dpsi_rate_full)}`",
                    f"`{fmt_val(dpsi_rate_ratio)}`",
                    f"`{fmt_val(dpsi_proj_rate_ctrl)}`",
                    f"`{fmt_val(dpsi_proj_rate_full)}`",
                    f"`{fmt_val(dpsi_proj_rate_ratio)}`",
                    f"`{fmt_val(jw_ctrl_abs)}`",
                    f"`{fmt_val(jw_full_abs)}`",
                    f"`{fmt_val(jw_ratio)}`",
                    f"`{fmt_val(s_leak_ctrl)}`",
                    f"`{fmt_val(s_leak_full)}`",
                    f"`{fmt_val(s_leak_ratio)}`",
                    f"`{fmt_val(mixed_ew2_ctrl)}`",
                    f"`{fmt_val(mixed_ew2_full)}`",
                    f"`{fmt_val(mixed_ew2_ratio)}`",
                    f"`{fmt_val(jw_mode_activity_proxy_ctrl)}`",
                    f"`{fmt_val(jw_mode_activity_proxy)}`",
                    f"`{fmt_val(jw_mode_proxy_ratio)}`",
                    f"`{row.get('full_closure_status', 'N/A')}`",
                    f"`{row.get('full_transport_closure_status', 'N/A')}`",
                    f"`{row.get('full_correlation_status', 'N/A')}`",
                    f"`{row.get('full_activity_status', 'N/A')}`",
                ]
            )
            + " |"
        )

    lines.append("")
    lines.append("## Controlled-vs-Full Activation Gates")
    lines.append("")
    lines.append("| Gate | Value | Status |")
    lines.append("| --- | ---: | :---: |")
    for condition, value, status in compare_gate_rows:
        lines.append(f"| `{condition}` | `{fmt_val(value)}` | `{status}` |")

    lines.append("")
    lines.append("## Correlation Signal Variation")
    lines.append("")
    lines.append("| Corr Metric | Span Across S |")
    lines.append("| --- | ---: |")
    for key, value in corr_signal_spans.items():
        lines.append(f"| `{key}` | `{fmt_val(value)}` |")

    lines.append("")
    lines.append("## Rate Proxy Correlations")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("| --- | ---: |")
    lines.append(
        f"| `corr_logS_full_max_abs_dpsi0_dt` | `{fmt_val(metrics['corr_logS_full_max_abs_dpsi0_dt'])}` |"
    )
    lines.append(
        f"| `corr_logS_full_max_dpsi0_dt` | `{fmt_val(metrics['corr_logS_full_max_dpsi0_dt'])}` |"
    )
    lines.append(
        f"| `corr_logS_full_final_psi_proj_span` | `{fmt_val(metrics['corr_logS_full_final_psi_proj_span'])}` |"
    )
    lines.append(
        f"| `corr_logS_full_max_abs_dpsi_proj_dt` | `{fmt_val(metrics['corr_logS_full_max_abs_dpsi_proj_dt'])}` |"
    )
    lines.append(
        f"| `corr_logS_full_max_dpsi_proj_dt` | `{fmt_val(metrics['corr_logS_full_max_dpsi_proj_dt'])}` |"
    )

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"production_report,{report_path}")
    print(f"overall_status,{overall}")


if __name__ == "__main__":
    main()
