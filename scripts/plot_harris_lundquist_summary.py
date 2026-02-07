#!/usr/bin/env python3
"""Create standard production plots from a Lundquist scan summary CSV."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def parse_float(row, key):
    value = row.get(key, "")
    if value is None or value == "":
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def read_rows(summary_csv):
    rows = list(csv.DictReader(summary_csv.open("r", encoding="utf-8")))
    parsed = []
    for row in rows:
        parsed.append(
            {
                "S": parse_float(row, "S"),
                "psi0": parse_float(row, "full_final_psi0_span"),
                "psi_proj": parse_float(row, "full_final_psi_proj_span"),
                "psi_w0": parse_float(row, "full_final_psi_w0_span"),
                "jw_ew_abs": abs(parse_float(row, "full_final_jw_ew")),
                "s_leak_abs": parse_float(row, "full_final_s_leak_abs"),
                "mixed_ew2": parse_float(row, "full_final_mixed_ew2"),
                "mixed_c2": parse_float(row, "full_final_mixed_c2"),
                "em_leak_w": parse_float(row, "full_final_em_leak_w"),
                "helicity_sub": parse_float(row, "full_final_helicity_sub"),
                "edotb_sub": parse_float(row, "full_final_edotb_sub"),
                "em_a2_mode_1": parse_float(row, "full_final_em_a2_mode_1"),
                "em_pi2_mode_1": parse_float(row, "full_final_em_pi2_mode_1"),
                "jw_mode_l2_1": parse_float(row, "full_final_jw_mode_l2_1"),
                "jw_mode_activity_proxy": parse_float(
                    row, "full_final_jw_mode_activity_proxy"
                ),
                "closure_status": row.get("full_closure_status", "N/A"),
                "transport_status": row.get("full_transport_closure_status", "N/A"),
                "correlation_status": row.get("full_correlation_status", "N/A"),
                "activity_status": row.get("full_activity_status", "N/A"),
            }
        )
    parsed.sort(key=lambda r: r["S"])
    return parsed


def finite_pair(xs, ys):
    out_x, out_y = [], []
    for x, y in zip(xs, ys):
        if math.isfinite(x) and math.isfinite(y):
            out_x.append(x)
            out_y.append(y)
    return out_x, out_y


def apply_log_y_if_positive(ax, ys):
    finite = [y for y in ys if math.isfinite(y)]
    if finite and min(finite) > 0.0:
        ax.set_yscale("log")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-csv", required=True, help="Path to lundquist_scan_summary.csv")
    parser.add_argument(
        "--output-dir",
        default="",
        help="Directory for output plots (default: <summary-csv parent>/plots)",
    )
    args = parser.parse_args()

    summary_csv = Path(args.summary_csv).resolve()
    if not summary_csv.exists():
        raise FileNotFoundError(f"Summary CSV not found: {summary_csv}")

    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else (summary_csv.parent / "plots").resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = read_rows(summary_csv)
    if not rows:
        raise RuntimeError("Summary CSV is empty.")

    s_vals = [r["S"] for r in rows]
    psi0 = [r["psi0"] for r in rows]
    psi_proj = [r["psi_proj"] for r in rows]
    psi_w0 = [r["psi_w0"] for r in rows]
    jw_ew_abs = [r["jw_ew_abs"] for r in rows]
    s_leak_abs = [r["s_leak_abs"] for r in rows]
    mixed_ew2 = [r["mixed_ew2"] for r in rows]
    mixed_c2 = [r["mixed_c2"] for r in rows]
    em_leak_w = [r["em_leak_w"] for r in rows]
    helicity_sub = [r["helicity_sub"] for r in rows]
    edotb_sub = [r["edotb_sub"] for r in rows]
    em_a2_mode_1 = [r["em_a2_mode_1"] for r in rows]
    em_pi2_mode_1 = [r["em_pi2_mode_1"] for r in rows]
    jw_mode_l2_1 = [r["jw_mode_l2_1"] for r in rows]
    jw_mode_activity_proxy = [r["jw_mode_activity_proxy"] for r in rows]

    fig1, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    x, y = finite_pair(s_vals, psi0)
    axs[0, 0].plot(x, y, marker="o", label="full_final_psi0_span")
    xp, yp = finite_pair(s_vals, psi_proj)
    if yp:
        axs[0, 0].plot(xp, yp, marker="s", label="full_final_psi_proj_span")
    xw, yw = finite_pair(s_vals, psi_w0)
    if yw:
        axs[0, 0].plot(xw, yw, marker="^", label="full_final_psi_w0_span")
    axs[0, 0].set_xscale("log")
    axs[0, 0].set_title("Reconnection Proxy")
    axs[0, 0].set_xlabel("Lundquist S")
    axs[0, 0].set_ylabel("psi0 span")
    axs[0, 0].grid(True, alpha=0.3)
    axs[0, 0].legend()

    x1, y1 = finite_pair(s_vals, jw_ew_abs)
    x2, y2 = finite_pair(s_vals, s_leak_abs)
    axs[0, 1].plot(x1, y1, marker="o", label="|JwEw|")
    axs[0, 1].plot(x2, y2, marker="s", label="S_leak_abs")
    axs[0, 1].set_xscale("log")
    apply_log_y_if_positive(axs[0, 1], y1 + y2)
    axs[0, 1].set_title("Leakage Work Channels")
    axs[0, 1].set_xlabel("Lundquist S")
    axs[0, 1].set_ylabel("magnitude")
    axs[0, 1].grid(True, alpha=0.3)
    axs[0, 1].legend()

    x1, y1 = finite_pair(s_vals, mixed_ew2)
    x2, y2 = finite_pair(s_vals, mixed_c2)
    axs[1, 0].plot(x1, y1, marker="o", label="mixed_ew2")
    axs[1, 0].plot(x2, y2, marker="s", label="mixed_c2")
    axs[1, 0].set_xscale("log")
    apply_log_y_if_positive(axs[1, 0], y1 + y2)
    axs[1, 0].set_title("Mixed-Sector Intensity")
    axs[1, 0].set_xlabel("Lundquist S")
    axs[1, 0].set_ylabel("magnitude")
    axs[1, 0].grid(True, alpha=0.3)
    axs[1, 0].legend()

    x1, y1 = finite_pair(s_vals, em_a2_mode_1)
    x2, y2 = finite_pair(s_vals, em_pi2_mode_1)
    x3, y3 = finite_pair(s_vals, jw_mode_l2_1)
    x4, y4 = finite_pair(s_vals, jw_mode_activity_proxy)
    axs[1, 1].plot(x1, y1, marker="o", label="em_a2_mode_1")
    axs[1, 1].plot(x2, y2, marker="s", label="em_pi2_mode_1")
    if y3:
        axs[1, 1].plot(x3, y3, marker="^", label="jw_mode_l2_1")
    if y4:
        axs[1, 1].plot(x4, y4, marker="d", label="jw_mode_activity_proxy")
    axs[1, 1].set_xscale("log")
    apply_log_y_if_positive(axs[1, 1], y1 + y2 + y3 + y4)
    axs[1, 1].set_title("Mode-1 Transfer Channels")
    axs[1, 1].set_xlabel("Lundquist S")
    axs[1, 1].set_ylabel("magnitude")
    axs[1, 1].grid(True, alpha=0.3)
    axs[1, 1].legend()

    primary_path = output_dir / "lundquist_primary_channels.png"
    fig1.savefig(primary_path, dpi=150)
    plt.close(fig1)

    fig2, axs2 = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)

    x, y = finite_pair(s_vals, em_leak_w)
    axs2[0].plot(x, y, marker="o", label="em_leak_w")
    axs2[0].axhline(0.0, color="k", linewidth=1.0, alpha=0.5)
    axs2[0].set_xscale("log")
    axs2[0].set_title("Projected EM Leakage")
    axs2[0].set_xlabel("Lundquist S")
    axs2[0].set_ylabel("signed value")
    axs2[0].grid(True, alpha=0.3)

    x, y = finite_pair(s_vals, helicity_sub)
    axs2[1].plot(x, y, marker="o", label="helicity_sub")
    axs2[1].axhline(0.0, color="k", linewidth=1.0, alpha=0.5)
    axs2[1].set_xscale("log")
    axs2[1].set_title("Subscale Helicity")
    axs2[1].set_xlabel("Lundquist S")
    axs2[1].set_ylabel("signed value")
    axs2[1].grid(True, alpha=0.3)

    x, y = finite_pair(s_vals, edotb_sub)
    axs2[2].plot(x, y, marker="o", label="edotb_sub")
    axs2[2].axhline(0.0, color="k", linewidth=1.0, alpha=0.5)
    axs2[2].set_xscale("log")
    axs2[2].set_title("Subscale E.B")
    axs2[2].set_xlabel("Lundquist S")
    axs2[2].set_ylabel("signed value")
    axs2[2].grid(True, alpha=0.3)

    subscale_path = output_dir / "lundquist_subscale_channels.png"
    fig2.savefig(subscale_path, dpi=150)
    plt.close(fig2)

    status_path = output_dir / "lundquist_status.csv"
    with status_path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(
            [
                "S",
                "closure_status",
                "transport_status",
                "correlation_status",
                "activity_status",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["S"],
                    row["closure_status"],
                    row["transport_status"],
                    row["correlation_status"],
                    row["activity_status"],
                ]
            )

    print(f"summary_csv,{summary_csv}")
    print(f"plot,{primary_path}")
    print(f"plot,{subscale_path}")
    print(f"status_csv,{status_path}")


if __name__ == "__main__":
    main()
