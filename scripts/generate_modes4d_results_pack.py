#!/usr/bin/env python3
"""Generate a combined results pack from production and ablation outputs."""

import argparse
import csv
import math
from datetime import datetime, timezone
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


def parse_scalar_float(raw):
    if raw is None or raw == "":
        return math.nan
    try:
        return float(raw)
    except ValueError:
        return math.nan


def fmt(value):
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


def load_csv_rows(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing CSV: {path}")
    rows = list(csv.DictReader(path.open("r", encoding="utf-8")))
    if not rows:
        raise RuntimeError(f"CSV has no rows: {path}")
    return rows


def load_status_map(path):
    rows = load_csv_rows(path)
    status_map = {}
    for row in rows:
        component = (row.get("component") or "").strip()
        status = (row.get("status") or "").strip()
        if component:
            status_map[component] = status
    if not status_map:
        raise RuntimeError(f"Status CSV has no component rows: {path}")
    return status_map


def load_topology_gate_rows(path):
    rows = load_csv_rows(path)
    gate_rows = []
    for row in rows:
        gate_rows.append(
            (
                row.get("gate", ""),
                parse_scalar_float(row.get("value", "")),
                parse_scalar_float(row.get("threshold", "")),
                (row.get("status", "") or "FAIL").strip(),
                (row.get("gate_class", "") or "blocking").strip(),
            )
        )
    return gate_rows


def build_lundquist_gate_table(rows):
    rows_sorted = sorted(rows, key=lambda r: parse_float(r, "S"))
    log_s = [math.log10(parse_float(r, "S")) for r in rows_sorted]

    def corr_from_col(col, abs_value=False):
        ys = []
        for row in rows_sorted:
            value = parse_float(row, col)
            if abs_value:
                value = abs(value)
            ys.append(value)
        mx = sum(log_s) / len(log_s)
        my = sum(ys) / len(ys)
        vx = sum((x - mx) * (x - mx) for x in log_s)
        vy = sum((y - my) * (y - my) for y in ys)
        if vx <= 0.0 or vy <= 0.0:
            return math.nan
        cov = sum((x - mx) * (y - my) for x, y in zip(log_s, ys))
        return cov / math.sqrt(vx * vy)

    metrics = {
        "corr_logS_full_final_psi0_span": corr_from_col("full_final_psi0_span"),
        "corr_logS_full_final_s_leak_abs": corr_from_col("full_final_s_leak_abs"),
        "corr_logS_full_final_jw_ew_abs": corr_from_col("full_final_jw_ew", abs_value=True),
        "corr_logS_full_final_mixed_ew2": corr_from_col("full_final_mixed_ew2"),
        "corr_logS_full_final_em_leak_w_abs": corr_from_col("full_final_em_leak_w", abs_value=True),
        "corr_logS_full_final_helicity_sub_abs": corr_from_col("full_final_helicity_sub", abs_value=True),
        "corr_logS_full_final_edotb_sub_abs": corr_from_col("full_final_edotb_sub", abs_value=True),
        "corr_logS_full_final_jw_mode_activity_proxy": corr_from_col(
            "full_final_jw_mode_activity_proxy"
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
    overall = "PASS"
    for key, mode, threshold in gate_specs:
        value = metrics.get(key, math.nan)
        if mode == "abs_min":
            status = gate_abs(value, threshold)
            condition = f"|{key}| >= {threshold:.1e}"
        else:
            status = gate_max(value, threshold)
            condition = f"{key} <= {threshold:.1e}"
        if status != "PASS":
            overall = "FAIL"
        gate_rows.append((condition, value, status))

    def min_finite_from_col(col):
        vals = [parse_float(r, col) for r in rows_sorted]
        vals = [v for v in vals if math.isfinite(v)]
        if not vals:
            return math.nan
        return min(vals)

    enhancement_specs = [
        ("min(full/controlled psi_w0_span) >= 1.2e+00", min_finite_from_col("ratio_full_over_controlled_psi_w0_span"), 1.2),
        (
            "min(full/controlled max_abs_dpsi_w0_dt) >= 1.2e+00",
            min_finite_from_col("ratio_full_over_controlled_max_abs_dpsi_w0_dt"),
            1.2,
        ),
    ]
    for condition, value, threshold in enhancement_specs:
        status = "PASS" if (math.isfinite(value) and value >= threshold) else "FAIL"
        if status != "PASS":
            overall = "FAIL"
        gate_rows.append((condition, value, status))
    return gate_rows, overall


def build_ablation_table(rows):
    rows_sorted = sorted(rows, key=lambda r: r["scenario"])
    overall = "PASS"
    table_rows = []
    for row in rows_sorted:
        if row["scenario"] == "baseline_full4d":
            continue
        r_jw = parse_float(row, "ratio_abs_final_jw_ew")
        r_sl = parse_float(row, "ratio_final_s_leak_abs")
        r_me = parse_float(row, "ratio_final_mixed_ew2")
        r_mc = parse_float(row, "ratio_final_mixed_c2")
        status = "INFO"
        details = "informational"
        if row["scenario"] == "ablate_mixed_seed":
            status = "PASS" if (r_me <= 0.5 and r_mc <= 0.5) else "FAIL"
            details = f"ratio_mixed_ew2={fmt(r_me)}, ratio_mixed_c2={fmt(r_mc)}"
        if row["scenario"] == "ablate_transverse_force":
            status = "PASS" if (r_jw <= 0.5 and r_sl <= 0.5) else "FAIL"
            details = f"ratio_abs_jw_ew={fmt(r_jw)}, ratio_s_leak={fmt(r_sl)}"
        if status == "FAIL":
            overall = "FAIL"
        table_rows.append((row["scenario"], row["description"], r_jw, r_sl, r_me, r_mc, status, details))
    return table_rows, overall


def write_summary_csv(path, lundquist_status, ablation_status, topology_status, overall_status):
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["component", "status"])
        writer.writerow(["lundquist_gates", lundquist_status])
        writer.writerow(["ablation_gates", ablation_status])
        writer.writerow(["topology_gates", topology_status])
        writer.writerow(["overall", overall_status])


def write_markdown_report(
    path,
    lundquist_summary_csv,
    lundquist_report,
    topology_report,
    ablation_summary_csv,
    ablation_report,
    panel_png,
    lundquist_gate_rows,
    topology_gate_rows,
    ablation_rows,
    lundquist_status,
    topology_status,
    ablation_status,
    overall_status,
):
    lines = []
    lines.append("# Modes4D Results Pack")
    lines.append("")
    lines.append(f"- Generated (UTC): {datetime.now(timezone.utc).isoformat()}")
    lines.append(f"- Overall status: **{overall_status}**")
    lines.append(f"- Lundquist gate status: **{lundquist_status}**")
    lines.append(f"- Topology gate status: **{topology_status}**")
    lines.append(f"- Ablation gate status: **{ablation_status}**")
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append(f"- Lundquist summary CSV: `{lundquist_summary_csv}`")
    lines.append(f"- Lundquist report: `{lundquist_report}`")
    lines.append(f"- Topology report: `{topology_report}`")
    lines.append(f"- Ablation summary CSV: `{ablation_summary_csv}`")
    lines.append(f"- Ablation report: `{ablation_report}`")
    lines.append(f"- Combined panel: `{panel_png}`")
    lines.append("")
    lines.append("## Lundquist Gates")
    lines.append("")
    lines.append("| Gate | Value | Status |")
    lines.append("| :--- | ---: | :---: |")
    for cond, value, status in lundquist_gate_rows:
        lines.append(f"| `{cond}` | `{fmt(value)}` | `{status}` |")
    lines.append("")
    lines.append("## Topology Gates")
    lines.append("")
    lines.append("| Gate | Value | Threshold | Status | Class |")
    lines.append("| :--- | ---: | ---: | :---: | :--- |")
    for gate, value, threshold, status, gate_class in topology_gate_rows:
        lines.append(
            f"| `{gate}` | `{fmt(value)}` | `{fmt(threshold)}` | `{status}` | `{gate_class}` |"
        )
    lines.append("")
    lines.append("## Ablation Causality Checks")
    lines.append("")
    lines.append("| scenario | ratio_abs_jw_ew | ratio_s_leak | ratio_mixed_ew2 | ratio_mixed_c2 | status | details |")
    lines.append("| :--- | ---: | ---: | ---: | ---: | :---: | :--- |")
    for scenario, desc, r_jw, r_sl, r_me, r_mc, status, details in ablation_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{scenario}`",
                    f"`{fmt(r_jw)}`",
                    f"`{fmt(r_sl)}`",
                    f"`{fmt(r_me)}`",
                    f"`{fmt(r_mc)}`",
                    f"`{status}`",
                    f"{desc}; {details}",
                ]
            )
            + " |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_panel(
    panel_path,
    lundquist_primary_png,
    lundquist_subscale_png,
    topology_png,
    ablation_png,
    lundquist_status,
    topology_status,
    ablation_status,
    overall_status,
):
    fig, axs = plt.subplots(2, 3, figsize=(20, 10), constrained_layout=True)

    axs[0, 0].imshow(plt.imread(lundquist_primary_png))
    axs[0, 0].set_title("Lundquist Primary Channels")
    axs[0, 0].axis("off")

    axs[0, 1].imshow(plt.imread(lundquist_subscale_png))
    axs[0, 1].set_title("Lundquist Subscale Channels")
    axs[0, 1].axis("off")

    axs[0, 2].imshow(plt.imread(topology_png))
    axs[0, 2].set_title("Controlled vs Full Topology")
    axs[0, 2].axis("off")

    axs[1, 0].imshow(plt.imread(ablation_png))
    axs[1, 0].set_title("Ablation Channels")
    axs[1, 0].axis("off")

    axs[1, 1].axis("off")
    summary_text = "\n".join(
        [
            "Modes4D Results Pack",
            "",
            f"Lundquist gates: {lundquist_status}",
            f"Topology gates: {topology_status}",
            f"Ablation gates: {ablation_status}",
            f"Overall: {overall_status}",
            "",
            "Novelty evidence bundle:",
            "- scan-level S trends with ledger closure",
            "- controlled-vs-full topology comparison",
            "- channel-causality ablations",
            "- combined report/plots for PR review",
        ]
    )
    axs[1, 1].text(0.0, 0.95, summary_text, va="top", fontsize=12)
    axs[1, 2].axis("off")

    fig.savefig(panel_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lundquist-summary-csv", required=True)
    parser.add_argument("--lundquist-report", required=True)
    parser.add_argument("--lundquist-primary-plot", required=True)
    parser.add_argument("--lundquist-subscale-plot", required=True)
    parser.add_argument("--topology-report", required=True)
    parser.add_argument("--topology-plot", required=True)
    parser.add_argument("--topology-status-csv", required=True)
    parser.add_argument("--topology-gates-csv", required=True)
    parser.add_argument("--ablation-summary-csv", required=True)
    parser.add_argument("--ablation-report", required=True)
    parser.add_argument("--ablation-plot", required=True)
    parser.add_argument("--output-dir", default="modes4d_results_pack")
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    lundquist_summary_csv = Path(args.lundquist_summary_csv).resolve()
    lundquist_report = Path(args.lundquist_report).resolve()
    lundquist_primary_png = Path(args.lundquist_primary_plot).resolve()
    lundquist_subscale_png = Path(args.lundquist_subscale_plot).resolve()
    topology_report = Path(args.topology_report).resolve()
    topology_png = Path(args.topology_plot).resolve()
    topology_status_csv = Path(args.topology_status_csv).resolve()
    topology_gates_csv = Path(args.topology_gates_csv).resolve()
    ablation_summary_csv = Path(args.ablation_summary_csv).resolve()
    ablation_report = Path(args.ablation_report).resolve()
    ablation_png = Path(args.ablation_plot).resolve()

    for path in [
        lundquist_summary_csv,
        lundquist_report,
        lundquist_primary_png,
        lundquist_subscale_png,
        topology_report,
        topology_png,
        topology_status_csv,
        topology_gates_csv,
        ablation_summary_csv,
        ablation_report,
        ablation_png,
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")

    lundquist_rows = load_csv_rows(lundquist_summary_csv)
    topology_status_rows = load_status_map(topology_status_csv)
    topology_gate_rows = load_topology_gate_rows(topology_gates_csv)
    ablation_rows_src = load_csv_rows(ablation_summary_csv)

    lundquist_gate_rows, lundquist_status = build_lundquist_gate_table(lundquist_rows)
    topology_status = topology_status_rows.get("overall", "FAIL")
    ablation_rows, ablation_status = build_ablation_table(ablation_rows_src)
    overall_status = (
        "PASS"
        if (lundquist_status == "PASS" and topology_status == "PASS" and ablation_status == "PASS")
        else "FAIL"
    )

    panel_png = output_dir / "modes4d_results_panel.png"
    summary_csv = output_dir / "modes4d_results_status.csv"
    report_md = output_dir / "modes4d_results_pack.md"

    render_panel(
        panel_png,
        lundquist_primary_png,
        lundquist_subscale_png,
        topology_png,
        ablation_png,
        lundquist_status,
        topology_status,
        ablation_status,
        overall_status,
    )
    write_summary_csv(summary_csv, lundquist_status, ablation_status, topology_status, overall_status)
    write_markdown_report(
        report_md,
        lundquist_summary_csv,
        lundquist_report,
        topology_report,
        ablation_summary_csv,
        ablation_report,
        panel_png,
        lundquist_gate_rows,
        topology_gate_rows,
        ablation_rows,
        lundquist_status,
        topology_status,
        ablation_status,
        overall_status,
    )

    print(f"results_panel,{panel_png}")
    print(f"results_status_csv,{summary_csv}")
    print(f"results_report,{report_md}")
    print(f"results_status,{overall_status}")


if __name__ == "__main__":
    main()
