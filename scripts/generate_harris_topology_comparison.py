#!/usr/bin/env python3
"""Generate controlled-vs-full topology comparison artifacts for a selected S."""

import argparse
import csv
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def fmt(value):
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def finite(values):
    return [v for v in values if math.isfinite(v)]


def series_max(values):
    vals = finite(values)
    if not vals:
        return math.nan
    return max(vals)


def safe_ratio(num, den):
    if not math.isfinite(num) or not math.isfinite(den) or den <= 0.0:
        return math.nan
    return num / den


def parse_hst(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    header = None
    for line in lines:
        if line.startswith("# [1]="):
            header = line
            break
    if header is None:
        raise RuntimeError(f"Could not find history header in {path}")

    labels = []
    for token in header.replace("#", "").split():
        if token.startswith("[") and "]=" in token:
            labels.append(token.split("]=", 1)[1])

    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])

    if not rows:
        raise RuntimeError(f"No history data rows in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


def maybe_col(cols, key):
    return cols.get(key)


def abs_series(values):
    if values is None:
        return None
    return [abs(v) if math.isfinite(v) else math.nan for v in values]


def derive_dpsi_dt(times, psi):
    if times is None or psi is None or len(times) != len(psi) or len(times) < 2:
        return ([], [])
    rates_t = []
    rates = []
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        rates_t.append(times[i])
        rates.append((psi[i] - psi[i - 1]) / dt)
    return rates_t, rates


def choose_scan_row(rows, s_value):
    rows_sorted = sorted(rows, key=lambda r: parse_float(r.get("S")))
    if not rows_sorted:
        raise RuntimeError("Summary CSV has no rows.")
    if s_value is None:
        return rows_sorted[0]
    best = min(rows_sorted, key=lambda r: abs(parse_float(r.get("S")) - s_value))
    return best


def gate_min(value, threshold):
    return "PASS" if (math.isfinite(value) and value >= threshold) else "FAIL"


def gate_max(value, threshold):
    return "PASS" if (math.isfinite(value) and value <= threshold) else "FAIL"


def write_status_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["component", "status"])
        for component, status in rows:
            writer.writerow([component, status])


def write_gate_details_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["gate", "value", "threshold", "status"])
        for gate, value, threshold, status in rows:
            writer.writerow([gate, fmt(value), f"{threshold:.6e}", status])


def render_plot(path, ctrl, full):
    times_c = ctrl["time"]
    times_f = full["time"]
    psi_c = ctrl["psi0"]
    psi_f = full["psi0"]
    dpsi_t_c, dpsi_c = derive_dpsi_dt(times_c, psi_c)
    dpsi_t_f, dpsi_f = derive_dpsi_dt(times_f, psi_f)

    fig, axs = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    axs[0, 0].plot(times_c, psi_c, label="controlled", lw=2)
    axs[0, 0].plot(times_f, psi_f, label="full", lw=2)
    axs[0, 0].set_title("psi0_span")
    axs[0, 0].set_xlabel("time")
    axs[0, 0].set_ylabel("psi0_span")
    axs[0, 0].legend()

    axs[0, 1].plot(dpsi_t_c, abs_series(dpsi_c), label="controlled", lw=2)
    axs[0, 1].plot(dpsi_t_f, abs_series(dpsi_f), label="full", lw=2)
    axs[0, 1].set_title("|d(psi0_span)/dt|")
    axs[0, 1].set_xlabel("time")
    axs[0, 1].set_ylabel("rate")
    axs[0, 1].legend()

    eps = 1.0e-30
    epar_c = [max(v, eps) for v in abs_series(ctrl["epar2"])]
    epar_f = [max(v, eps) for v in abs_series(full["epar2"])]
    mix_c = [max(v, eps) for v in abs_series(ctrl["mixed_ew2"])]
    mix_f = [max(v, eps) for v in abs_series(full["mixed_ew2"])]
    axs[1, 0].semilogy(times_c, epar_c, label="|E_parallel|^2 controlled", lw=2)
    axs[1, 0].semilogy(times_f, epar_f, label="|E_parallel|^2 full", lw=2)
    axs[1, 0].semilogy(times_c, mix_c, "--", label="mixed_ew2 controlled", lw=2)
    axs[1, 0].semilogy(times_f, mix_f, "--", label="mixed_ew2 full", lw=2)
    axs[1, 0].set_title("Topology Proxies")
    axs[1, 0].set_xlabel("time")
    axs[1, 0].set_ylabel("magnitude")
    axs[1, 0].legend()

    jw_c = [max(v, eps) for v in abs_series(ctrl["jw_ew"])]
    jw_f = [max(v, eps) for v in abs_series(full["jw_ew"])]
    sl_c = [max(v, eps) for v in ctrl["s_leak_abs"]]
    sl_f = [max(v, eps) for v in full["s_leak_abs"]]
    axs[1, 1].semilogy(times_c, jw_c, label="|JwEw| controlled", lw=2)
    axs[1, 1].semilogy(times_f, jw_f, label="|JwEw| full", lw=2)
    axs[1, 1].semilogy(times_c, sl_c, "--", label="S_leak_abs controlled", lw=2)
    axs[1, 1].semilogy(times_f, sl_f, "--", label="S_leak_abs full", lw=2)
    axs[1, 1].set_title("Leakage/Work Channels")
    axs[1, 1].set_xlabel("time")
    axs[1, 1].set_ylabel("magnitude")
    axs[1, 1].legend()

    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--scan-output-dir", required=True)
    parser.add_argument("--s-value", type=float, default=None)
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--full-min-jw-ew-abs", type=float, default=1.0e-14)
    parser.add_argument("--full-min-s-leak-abs", type=float, default=1.0e-13)
    parser.add_argument("--full-min-mixed-ew2", type=float, default=1.0e-2)
    parser.add_argument("--controlled-max-jw-ew-abs", type=float, default=1.0e-12)
    parser.add_argument("--controlled-max-s-leak-abs", type=float, default=1.0e-12)
    parser.add_argument("--max-abs-dpsi0-enhancement-min", type=float, default=1.0)
    args = parser.parse_args()

    summary_csv = Path(args.summary_csv).resolve()
    scan_output_dir = Path(args.scan_output_dir).resolve()
    if not summary_csv.exists():
        raise FileNotFoundError(f"Missing summary CSV: {summary_csv}")
    if not scan_output_dir.exists():
        raise FileNotFoundError(f"Missing scan output dir: {scan_output_dir}")

    rows = list(csv.DictReader(summary_csv.open("r", encoding="utf-8")))
    row = choose_scan_row(rows, args.s_value)
    selected_s = parse_float(row.get("S"))
    tag = str(selected_s).replace(".", "p")
    case_dir = scan_output_dir / f"S_{tag}"
    if not case_dir.exists():
        raise FileNotFoundError(f"Missing selected S case directory: {case_dir}")

    controlled_hst = case_dir / "harris_controlled.out1.hst"
    full_hst = case_dir / "harris_full.out1.hst"
    if not controlled_hst.exists() or not full_hst.exists():
        raise FileNotFoundError(
            f"Missing hst file(s): {controlled_hst} or {full_hst}"
        )

    outdir = (
        (scan_output_dir / "topology")
        if not args.output_dir
        else Path(args.output_dir).resolve()
    )
    outdir.mkdir(parents=True, exist_ok=True)

    cols_c = parse_hst(controlled_hst)
    cols_f = parse_hst(full_hst)
    ctrl = {
        "time": maybe_col(cols_c, "time"),
        "psi0": maybe_col(cols_c, "m4d_psi0_span"),
        "jw_ew": maybe_col(cols_c, "m4d_int_jw_ew"),
        "s_leak_abs": maybe_col(cols_c, "m4d_int_s_leak_abs") or abs_series(maybe_col(cols_c, "m4d_int_s_leak")),
        "mixed_ew2": maybe_col(cols_c, "m4d_mixed_ew2"),
        "epar2": maybe_col(cols_c, "m4d_brane_epar2"),
    }
    full = {
        "time": maybe_col(cols_f, "time"),
        "psi0": maybe_col(cols_f, "m4d_psi0_span"),
        "jw_ew": maybe_col(cols_f, "m4d_int_jw_ew"),
        "s_leak_abs": maybe_col(cols_f, "m4d_int_s_leak_abs") or abs_series(maybe_col(cols_f, "m4d_int_s_leak")),
        "mixed_ew2": maybe_col(cols_f, "m4d_mixed_ew2"),
        "epar2": maybe_col(cols_f, "m4d_brane_epar2"),
    }

    dpsi_t_c, dpsi_c = derive_dpsi_dt(ctrl["time"], ctrl["psi0"])
    dpsi_t_f, dpsi_f = derive_dpsi_dt(full["time"], full["psi0"])
    del dpsi_t_c, dpsi_t_f

    ctrl_max_abs_jw_ew = series_max(abs_series(ctrl["jw_ew"]))
    ctrl_max_s_leak_abs = series_max(ctrl["s_leak_abs"])
    full_max_abs_jw_ew = series_max(abs_series(full["jw_ew"]))
    full_max_s_leak_abs = series_max(full["s_leak_abs"])
    full_max_mixed_ew2 = series_max(full["mixed_ew2"])
    ctrl_max_abs_dpsi = series_max(abs_series(dpsi_c))
    full_max_abs_dpsi = series_max(abs_series(dpsi_f))
    ratio_max_abs_dpsi = safe_ratio(full_max_abs_dpsi, ctrl_max_abs_dpsi)

    gate_rows = []
    gate_rows.append(
        (
            "max|JwEw|_controlled <= tol",
            ctrl_max_abs_jw_ew,
            args.controlled_max_jw_ew_abs,
            gate_max(ctrl_max_abs_jw_ew, args.controlled_max_jw_ew_abs),
        )
    )
    gate_rows.append(
        (
            "max(S_leak_abs)_controlled <= tol",
            ctrl_max_s_leak_abs,
            args.controlled_max_s_leak_abs,
            gate_max(ctrl_max_s_leak_abs, args.controlled_max_s_leak_abs),
        )
    )
    gate_rows.append(
        (
            "max|JwEw|_full >= floor",
            full_max_abs_jw_ew,
            args.full_min_jw_ew_abs,
            gate_min(full_max_abs_jw_ew, args.full_min_jw_ew_abs),
        )
    )
    gate_rows.append(
        (
            "max(S_leak_abs)_full >= floor",
            full_max_s_leak_abs,
            args.full_min_s_leak_abs,
            gate_min(full_max_s_leak_abs, args.full_min_s_leak_abs),
        )
    )
    gate_rows.append(
        (
            "max(mixed_ew2)_full >= floor",
            full_max_mixed_ew2,
            args.full_min_mixed_ew2,
            gate_min(full_max_mixed_ew2, args.full_min_mixed_ew2),
        )
    )
    gate_rows.append(
        (
            "max|dpsi/dt|_full / controlled >= floor",
            ratio_max_abs_dpsi,
            args.max_abs_dpsi0_enhancement_min,
            gate_min(ratio_max_abs_dpsi, args.max_abs_dpsi0_enhancement_min),
        )
    )

    controlled_inactivity = (
        "PASS"
        if all(
            status == "PASS"
            for _, _, _, status in gate_rows[:2]
        )
        else "FAIL"
    )
    full_activation = (
        "PASS"
        if all(
            status == "PASS"
            for _, _, _, status in gate_rows[2:5]
        )
        else "FAIL"
    )
    rate_enhancement = gate_rows[5][3]

    ledger_fields = [
        "full_closure_status",
        "full_transport_closure_status",
        "full_correlation_status",
        "full_activity_status",
    ]
    ledger_status = (
        "PASS" if all(row.get(field, "FAIL") == "PASS" for field in ledger_fields) else "FAIL"
    )
    overall = (
        "PASS"
        if all(
            status == "PASS"
            for status in [controlled_inactivity, full_activation, rate_enhancement, ledger_status]
        )
        else "FAIL"
    )

    plot_path = outdir / "topology_timeseries.png"
    status_csv = outdir / "topology_status.csv"
    gates_csv = outdir / "topology_gate_details.csv"
    report_md = outdir / "topology_comparison_report.md"

    render_plot(plot_path, ctrl, full)
    write_status_csv(
        status_csv,
        [
            ("controlled_inactivity", controlled_inactivity),
            ("full_activation", full_activation),
            ("rate_enhancement", rate_enhancement),
            ("ledger_status", ledger_status),
            ("overall", overall),
        ],
    )
    write_gate_details_csv(gates_csv, gate_rows)

    lines = [
        "# Harris Topology Comparison",
        "",
        f"- Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"- Selected S: `{selected_s:.0f}`",
        f"- Overall status: **{overall}**",
        f"- Controlled inactivity: **{controlled_inactivity}**",
        f"- Full activation: **{full_activation}**",
        f"- Rate enhancement: **{rate_enhancement}**",
        f"- Ledger status (from scan row): **{ledger_status}**",
        "",
        "## Artifacts",
        "",
        f"- Controlled history: `{controlled_hst}`",
        f"- Full history: `{full_hst}`",
        f"- Topology plot: `{plot_path}`",
        f"- Topology status CSV: `{status_csv}`",
        f"- Topology gate details CSV: `{gates_csv}`",
        "",
        "## Gates",
        "",
        "| Gate | Value | Threshold | Status |",
        "| :--- | ---: | ---: | :---: |",
    ]
    for gate, value, threshold, status in gate_rows:
        lines.append(f"| `{gate}` | `{fmt(value)}` | `{threshold:.6e}` | `{status}` |")
    lines.append("")
    lines.append("## Ledger Row Statuses")
    lines.append("")
    lines.append("| Field | Value |")
    lines.append("| :--- | :---: |")
    for field in ledger_fields:
        lines.append(f"| `{field}` | `{row.get(field, 'N/A')}` |")
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"topology_plot,{plot_path}")
    print(f"topology_status_csv,{status_csv}")
    print(f"topology_gate_details_csv,{gates_csv}")
    print(f"topology_report,{report_md}")
    print(f"topology_status,{overall}")


if __name__ == "__main__":
    main()
