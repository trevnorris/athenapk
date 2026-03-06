#!/usr/bin/env python3
"""Decompose EM mode-0 pi transport residuals into dQ/dt, div, and src pieces."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


CHANNELS = (
    ("pi0", "m4d_pi0_mode_0", "m4d_int_divpi0_mode0", "m4d_int_srcpi0_mode0"),
    ("piy", "m4d_piy_mode_0", "m4d_int_divpiy_mode0", "m4d_int_srcpiy_mode0"),
    ("piw", "m4d_piw_mode_0", "m4d_int_divpiw_mode0", "m4d_int_srcpiw_mode0"),
)


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


def safe_float(value: str | None) -> float:
    if value is None:
        return math.nan
    text = value.strip()
    if not text:
        return math.nan
    try:
        return float(text)
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


def rms_abs(values):
    if not values:
        return math.nan
    return math.sqrt(sum(v * v for v in values) / len(values))


def max_abs(values):
    if not values:
        return math.nan
    return max(abs(v) for v in values)


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


def compute_balance_series(times, quantity, int_div, int_src):
    t_out = []
    dqdt = []
    ddivdt = []
    dsrcdt = []
    residual = []
    if (
        len(times) < 3
        or quantity is None
        or int_div is None
        or len(quantity) != len(times)
        or len(int_div) != len(times)
    ):
        return (t_out, dqdt, ddivdt, dsrcdt, residual)
    if int_src is not None and len(int_src) != len(times):
        return (t_out, dqdt, ddivdt, dsrcdt, residual)

    for i in range(2, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        q = (quantity[i] - quantity[i - 1]) / dt
        div = (int_div[i] - int_div[i - 1]) / dt
        src = 0.0
        if int_src is not None:
            src = (int_src[i] - int_src[i - 1]) / dt
        t_out.append(times[i])
        dqdt.append(q)
        ddivdt.append(div)
        dsrcdt.append(src)
        residual.append(q + div - src)
    return (t_out, dqdt, ddivdt, dsrcdt, residual)


def dominant_term(max_dqdt: float, max_ddivdt: float, max_dsrcdt: float) -> str:
    items = (
        ("dqdt", max_dqdt),
        ("ddivdt", max_ddivdt),
        ("dsrcdt", max_dsrcdt),
    )
    ranked = [item for item in items if math.isfinite(item[1])]
    if not ranked:
        return "nan"
    return max(ranked, key=lambda item: item[1])[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Case glob (absolute or relative)")
    parser.add_argument("--baseline-case", default="baseline_pi_transport")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--summary-csv", required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows = []
    summary_rows = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        hst = case_dir / "harris_full.out1.hst"
        if not run_log.exists() or not hst.exists():
            continue

        try:
            _, full = parse_case_rows(run_log)
        except Exception:
            full = {}
        cols = hsm.parse_hst(hst)
        times = cols.get("time", [])

        case_channels = []
        for channel, q_key, div_key, src_key in CHANNELS:
            series = compute_balance_series(
                times,
                cols.get(q_key),
                cols.get(div_key),
                cols.get(src_key),
            )
            t_out, dqdt, ddivdt, dsrcdt, residual = series
            if not residual:
                continue

            peak_index = max(range(len(residual)), key=lambda idx: abs(residual[idx]))
            max_dqdt = max_abs(dqdt)
            max_ddivdt = max_abs(ddivdt)
            max_dsrcdt = max_abs(dsrcdt)
            max_residual = max_abs(residual)
            row = {
                "case_name": case_dir.name,
                "channel": channel,
                "transport_status": full.get(f"{channel}_transport_status", ""),
                "ledger_status": full.get("em_bulk_ledger_status", ""),
                "max_abs_dqdt": max_dqdt,
                "rms_abs_dqdt": rms_abs(dqdt),
                "max_abs_ddivdt": max_ddivdt,
                "rms_abs_ddivdt": rms_abs(ddivdt),
                "max_abs_dsrcdt": max_dsrcdt,
                "rms_abs_dsrcdt": rms_abs(dsrcdt),
                "max_abs_residual": max_residual,
                "rms_abs_residual": rms_abs(residual),
                "final_residual": residual[-1],
                "peak_time": t_out[peak_index],
                "peak_dqdt": dqdt[peak_index],
                "peak_ddivdt": ddivdt[peak_index],
                "peak_dsrcdt": dsrcdt[peak_index],
                "peak_residual": residual[peak_index],
                "corr_residual_dqdt": pearson(residual, dqdt),
                "corr_residual_ddivdt": pearson(residual, ddivdt),
                "corr_residual_negdsrcdt": pearson(residual, [-x for x in dsrcdt]),
                "dominant_term": dominant_term(max_dqdt, max_ddivdt, max_dsrcdt),
            }
            rows.append(row)
            case_channels.append(row)

        if case_channels:
            top = max(case_channels, key=lambda row: row["max_abs_residual"])
            summary_rows.append(
                {
                    "case_name": case_dir.name,
                    "top_channel": top["channel"],
                    "top_max_abs_residual": top["max_abs_residual"],
                    "top_transport_status": top["transport_status"],
                    "top_ledger_status": top["ledger_status"],
                    "top_dominant_term": top["dominant_term"],
                    "top_peak_time": top["peak_time"],
                    "top_peak_dqdt": top["peak_dqdt"],
                    "top_peak_ddivdt": top["peak_ddivdt"],
                    "top_peak_dsrcdt": top["peak_dsrcdt"],
                    "top_peak_residual": top["peak_residual"],
                }
            )

    baseline = next((row for row in rows if row["case_name"] == args.baseline_case), None)
    baseline_by_channel = {
        row["channel"]: row["max_abs_residual"]
        for row in rows
        if row["case_name"] == args.baseline_case
    }
    for row in rows:
        row["residual_reduction_vs_baseline"] = (
            baseline_by_channel.get(row["channel"], math.nan) / row["max_abs_residual"]
            if math.isfinite(baseline_by_channel.get(row["channel"], math.nan))
            and math.isfinite(row["max_abs_residual"])
            and row["max_abs_residual"] != 0.0
            else math.nan
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_name",
        "channel",
        "transport_status",
        "ledger_status",
        "max_abs_dqdt",
        "rms_abs_dqdt",
        "max_abs_ddivdt",
        "rms_abs_ddivdt",
        "max_abs_dsrcdt",
        "rms_abs_dsrcdt",
        "max_abs_residual",
        "rms_abs_residual",
        "final_residual",
        "peak_time",
        "peak_dqdt",
        "peak_ddivdt",
        "peak_dsrcdt",
        "peak_residual",
        "corr_residual_dqdt",
        "corr_residual_ddivdt",
        "corr_residual_negdsrcdt",
        "dominant_term",
        "residual_reduction_vs_baseline",
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
        "top_max_abs_residual",
        "top_transport_status",
        "top_ledger_status",
        "top_dominant_term",
        "top_peak_time",
        "top_peak_dqdt",
        "top_peak_ddivdt",
        "top_peak_dsrcdt",
        "top_peak_residual",
    ]
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    out_md = Path(args.out_md)
    lines = [
        "# Step-151 EM Pi Term Isolation",
        "",
        f"Parsed `{len(rows)}` channel rows from `{args.glob}`.",
        "",
        "Residual definition per channel:",
        "",
        "`residual = dQ/dt + d(div-int)/dt - d(src-int)/dt`",
        "",
    ]
    if baseline:
        lines.append(f"Baseline case: `{args.baseline_case}`.")
        lines.append("")
    lines.append(
        "| case | channel | transport | ledger | max|dq/dt| | max|ddiv/dt| | max|dsrc/dt| | max|residual| | reduction vs baseline | dominant term | peak t | peak dq/dt | peak ddiv/dt | peak dsrc/dt | peak residual |"
    )
    lines.append(
        "|:---|:---|:---|:---|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|"
    )
    for row in rows:
        lines.append(
            f"| {row['case_name']} | {row['channel']} | {row['transport_status']} | "
            f"{row['ledger_status']} | {fmt(row['max_abs_dqdt'])} | "
            f"{fmt(row['max_abs_ddivdt'])} | {fmt(row['max_abs_dsrcdt'])} | "
            f"{fmt(row['max_abs_residual'])} | {fmt(row['residual_reduction_vs_baseline'])} | "
            f"{row['dominant_term']} | {fmt(row['peak_time'])} | "
            f"{fmt(row['peak_dqdt'])} | {fmt(row['peak_ddivdt'])} | "
            f"{fmt(row['peak_dsrcdt'])} | {fmt(row['peak_residual'])} |"
        )

    ranked_baseline = [
        row for row in rows if row["case_name"] == args.baseline_case
    ]
    if ranked_baseline:
        ranked_baseline.sort(key=lambda row: row["max_abs_residual"], reverse=True)
        lines.extend(
            [
                "",
                "## Baseline Ranking",
                "",
                "| channel | max|residual| | dominant term | corr(res,dq/dt) | corr(res,ddiv/dt) | corr(res,-dsrc/dt) |",
                "|:---|---:|:---|---:|---:|---:|",
            ]
        )
        for row in ranked_baseline:
            lines.append(
                f"| {row['channel']} | {fmt(row['max_abs_residual'])} | {row['dominant_term']} | "
                f"{fmt(row['corr_residual_dqdt'])} | {fmt(row['corr_residual_ddivdt'])} | "
                f"{fmt(row['corr_residual_negdsrcdt'])} |"
            )

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{summary_csv}")
    print(f"channels_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
