#!/usr/bin/env python3
"""Split EM pi mode-0 source integrals into RHS and mixing partitions."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


CHANNELS = {
    "pi0": {
        "q_col": "m4d_pi0_mode_0",
        "div_col": "m4d_int_divpi0_mode0",
        "net_col": "m4d_int_srcpi0_mode0",
        "rhs_col": "m4d_int_srcpi0_mode0_rhs",
        "mix_col": "m4d_int_srcpi0_mode0_mix",
    },
    "piy": {
        "q_col": "m4d_piy_mode_0",
        "div_col": "m4d_int_divpiy_mode0",
        "net_col": "m4d_int_srcpiy_mode0",
        "rhs_col": "m4d_int_srcpiy_mode0_rhs",
        "mix_col": "m4d_int_srcpiy_mode0_mix",
    },
    "piw": {
        "q_col": "m4d_piw_mode_0",
        "div_col": "m4d_int_divpiw_mode0",
        "net_col": "m4d_int_srcpiw_mode0",
        "rhs_col": "m4d_int_srcpiw_mode0_rhs",
        "mix_col": "m4d_int_srcpiw_mode0_mix",
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


def rms_abs(values):
    if not values:
        return math.nan
    return math.sqrt(sum(v * v for v in values) / len(values))


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


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def safe_ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def dominant_partition(max_rhs: float, max_mix: float) -> str:
    if not math.isfinite(max_rhs) and not math.isfinite(max_mix):
        return "nan"
    if not math.isfinite(max_mix):
        return "rhs"
    if not math.isfinite(max_rhs):
        return "mix"
    return "rhs" if max_rhs >= max_mix else "mix"


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Case directory glob")
    parser.add_argument("--baseline-case", default="baseline_pi_partition")
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
        case_rows = []
        for channel, config in CHANNELS.items():
            t_q, q_rate = derivative_series(times, cols.get(config["q_col"]))
            t_div, div_rate = derivative_series(times, cols.get(config["div_col"]))
            t_net, net_rate = derivative_series(times, cols.get(config["net_col"]))
            t_rhs, rhs_rate = derivative_series(times, cols.get(config["rhs_col"]))
            t_mix, mix_rate = derivative_series(times, cols.get(config["mix_col"]))
            if not q_rate or not (
                len(t_q)
                == len(t_div)
                == len(t_net)
                == len(t_rhs)
                == len(t_mix)
                == len(q_rate)
                == len(div_rate)
                == len(net_rate)
                == len(rhs_rate)
                == len(mix_rate)
            ):
                continue

            partition_gap = [
                net - rhs - mix for net, rhs, mix in zip(net_rate, rhs_rate, mix_rate)
            ]
            legacy_residual = [
                q + div - net for q, div, net in zip(q_rate, div_rate, net_rate)
            ]
            corrected_residual = [
                q + div - rhs - mix
                for q, div, rhs, mix in zip(q_rate, div_rate, rhs_rate, mix_rate)
            ]

            max_net = max_abs(net_rate)
            max_rhs = max_abs(rhs_rate)
            max_mix = max_abs(mix_rate)
            max_gap = max_abs(partition_gap)
            max_legacy = max_abs(legacy_residual)
            max_corrected = max_abs(corrected_residual)
            row = {
                "case_name": case_dir.name,
                "channel": channel,
                "transport_status": full.get(f"{channel}_transport_status", ""),
                "ledger_status": full.get("em_bulk_ledger_status", ""),
                "max_abs_q_rate": max_abs(q_rate),
                "max_abs_div_rate": max_abs(div_rate),
                "max_abs_net_rate": max_net,
                "max_abs_rhs_rate": max_rhs,
                "max_abs_mix_rate": max_mix,
                "rhs_over_net": safe_ratio(max_rhs, max_net),
                "mix_over_net": safe_ratio(max_mix, max_net),
                "max_abs_partition_gap_rate": max_gap,
                "gap_over_net": safe_ratio(max_gap, max_net),
                "max_abs_legacy_transport_residual": max_legacy,
                "max_abs_corrected_transport_residual": max_corrected,
                "corrected_over_legacy": safe_ratio(max_corrected, max_legacy),
                "corr_net_rhs_plus_mix": pearson(
                    net_rate, [rhs + mix for rhs, mix in zip(rhs_rate, mix_rate)]
                ),
                "corr_legacy_negnet": pearson(legacy_residual, [-v for v in net_rate]),
                "corr_corrected_neggap": pearson(
                    corrected_residual, [-v for v in partition_gap]
                ),
                "dominant_partition": dominant_partition(max_rhs, max_mix),
            }
            rows.append(row)
            case_rows.append(row)

        if case_rows:
            top = max(case_rows, key=lambda row: row["max_abs_net_rate"])
            summary_rows.append(
                {
                    "case_name": case_dir.name,
                    "top_channel": top["channel"],
                    "top_transport_status": top["transport_status"],
                    "top_ledger_status": top["ledger_status"],
                    "top_max_abs_net_rate": top["max_abs_net_rate"],
                    "top_rhs_over_net": top["rhs_over_net"],
                    "top_mix_over_net": top["mix_over_net"],
                    "top_gap_over_net": top["gap_over_net"],
                    "top_dominant_partition": top["dominant_partition"],
                    "top_corrected_over_legacy": top["corrected_over_legacy"],
                }
            )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_name",
        "channel",
        "transport_status",
        "ledger_status",
        "max_abs_q_rate",
        "max_abs_div_rate",
        "max_abs_net_rate",
        "max_abs_rhs_rate",
        "max_abs_mix_rate",
        "rhs_over_net",
        "mix_over_net",
        "max_abs_partition_gap_rate",
        "gap_over_net",
        "max_abs_legacy_transport_residual",
        "max_abs_corrected_transport_residual",
        "corrected_over_legacy",
        "corr_net_rhs_plus_mix",
        "corr_legacy_negnet",
        "corr_corrected_neggap",
        "dominant_partition",
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
        "top_transport_status",
        "top_ledger_status",
        "top_max_abs_net_rate",
        "top_rhs_over_net",
        "top_mix_over_net",
        "top_gap_over_net",
        "top_dominant_partition",
        "top_corrected_over_legacy",
    ]
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    out_md = Path(args.out_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Step-153 EM Pi Source Partition",
        "",
        "This report splits the EM mode-0 `pi` source integrals into:",
        "",
        "- `net`: legacy post-update diagnostic (`pi_new - pi_old`)",
        "- `rhs`: explicit RHS contribution accumulated during the step",
        "- `mix`: dynamic-basis mixing contribution accumulated during the step",
        "",
        "The key check is whether `net ≈ rhs + mix`, and whether the dominant",
        "transport-driving source sits in the explicit RHS path or the mixing path.",
        "",
        "## Case Summary",
        "",
        "| Case | Top Channel | Net Max | RHS/Net | Mix/Net | Gap/Net | Dominant | Corrected/Legacy |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| {case_name} | {top_channel} | {top_max_abs_net_rate} | {top_rhs_over_net} | "
            "{top_mix_over_net} | {top_gap_over_net} | {top_dominant_partition} | "
            "{top_corrected_over_legacy} |".format(
                case_name=row["case_name"],
                top_channel=row["top_channel"],
                top_max_abs_net_rate=fmt(row["top_max_abs_net_rate"]),
                top_rhs_over_net=fmt(row["top_rhs_over_net"]),
                top_mix_over_net=fmt(row["top_mix_over_net"]),
                top_gap_over_net=fmt(row["top_gap_over_net"]),
                top_dominant_partition=row["top_dominant_partition"],
                top_corrected_over_legacy=fmt(row["top_corrected_over_legacy"]),
            )
        )

    lines.extend(
        [
            "",
            "## Channel Detail",
            "",
            "| Case | Channel | Net Max | RHS Max | Mix Max | Gap Max | RHS/Net | Mix/Net | Dominant | Corr(net,rhs+mix) |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
        ]
    )
    for row in rows:
        lines.append(
            "| {case_name} | {channel} | {max_abs_net_rate} | {max_abs_rhs_rate} | "
            "{max_abs_mix_rate} | {max_abs_partition_gap_rate} | {rhs_over_net} | "
            "{mix_over_net} | {dominant_partition} | {corr_net_rhs_plus_mix} |".format(
                case_name=row["case_name"],
                channel=row["channel"],
                max_abs_net_rate=fmt(row["max_abs_net_rate"]),
                max_abs_rhs_rate=fmt(row["max_abs_rhs_rate"]),
                max_abs_mix_rate=fmt(row["max_abs_mix_rate"]),
                max_abs_partition_gap_rate=fmt(row["max_abs_partition_gap_rate"]),
                rhs_over_net=fmt(row["rhs_over_net"]),
                mix_over_net=fmt(row["mix_over_net"]),
                dominant_partition=row["dominant_partition"],
                corr_net_rhs_plus_mix=fmt(row["corr_net_rhs_plus_mix"]),
            )
        )

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"detail_csv,{summary_csv}")


if __name__ == "__main__":
    main()
