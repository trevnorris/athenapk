#!/usr/bin/env python3
"""Decompose EM mode-0 pi mixing into w0-dynamic, w0-local-gradient, and lambda families."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


CHANNELS = {
    "pi0": {
        "net_col": "m4d_int_srcpi0_mode0",
        "mix_col": "m4d_int_srcpi0_mode0_mix",
        "w0_dynamic_col": "m4d_int_srcpi0_mode0_mix_w0_dynamic",
        "w0_local_gradient_col": "m4d_int_srcpi0_mode0_mix_w0_local_gradient",
        "lambda_col": "m4d_int_srcpi0_mode0_mix_lambda",
    },
    "piy": {
        "net_col": "m4d_int_srcpiy_mode0",
        "mix_col": "m4d_int_srcpiy_mode0_mix",
        "w0_dynamic_col": "m4d_int_srcpiy_mode0_mix_w0_dynamic",
        "w0_local_gradient_col": "m4d_int_srcpiy_mode0_mix_w0_local_gradient",
        "lambda_col": "m4d_int_srcpiy_mode0_mix_lambda",
    },
    "piw": {
        "net_col": "m4d_int_srcpiw_mode0",
        "mix_col": "m4d_int_srcpiw_mode0_mix",
        "w0_dynamic_col": "m4d_int_srcpiw_mode0_mix_w0_dynamic",
        "w0_local_gradient_col": "m4d_int_srcpiw_mode0_mix_w0_local_gradient",
        "lambda_col": "m4d_int_srcpiw_mode0_mix_lambda",
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


def dominant_family(max_w0_dynamic: float, max_w0_local_gradient: float, max_lambda: float) -> str:
    items = [
        ("w0_dynamic", max_w0_dynamic),
        ("w0_local_gradient", max_w0_local_gradient),
        ("lambda", max_lambda),
    ]
    ranked = [(name, value) for name, value in items if math.isfinite(value)]
    if not ranked:
        return "nan"
    return max(ranked, key=lambda item: item[1])[0]


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
    parser.add_argument("--baseline-case", default="baseline_mix_family")
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
            t_net, net_rate = derivative_series(times, cols.get(config["net_col"]))
            t_mix, mix_rate = derivative_series(times, cols.get(config["mix_col"]))
            t_w0_dynamic, w0_dynamic_rate = derivative_series(
                times, cols.get(config["w0_dynamic_col"])
            )
            t_w0_local_gradient, w0_local_gradient_rate = derivative_series(
                times, cols.get(config["w0_local_gradient_col"])
            )
            t_lambda, lambda_rate = derivative_series(times, cols.get(config["lambda_col"]))
            if not net_rate or not (
                len(t_net)
                == len(t_mix)
                == len(t_w0_dynamic)
                == len(t_w0_local_gradient)
                == len(t_lambda)
                == len(net_rate)
                == len(mix_rate)
                == len(w0_dynamic_rate)
                == len(w0_local_gradient_rate)
                == len(lambda_rate)
            ):
                continue

            family_sum = [
                wd + wg + lam
                for wd, wg, lam in zip(w0_dynamic_rate, w0_local_gradient_rate, lambda_rate)
            ]
            family_gap = [mix - fam for mix, fam in zip(mix_rate, family_sum)]

            max_net = max_abs(net_rate)
            max_mix = max_abs(mix_rate)
            max_w0_dynamic = max_abs(w0_dynamic_rate)
            max_w0_local_gradient = max_abs(w0_local_gradient_rate)
            max_lambda = max_abs(lambda_rate)
            max_gap = max_abs(family_gap)

            row = {
                "case_name": case_dir.name,
                "channel": channel,
                "transport_status": full.get(f"{channel}_transport_status", ""),
                "ledger_status": full.get("em_bulk_ledger_status", ""),
                "max_abs_net_rate": max_net,
                "max_abs_mix_rate": max_mix,
                "max_abs_w0_dynamic_rate": max_w0_dynamic,
                "max_abs_w0_local_gradient_rate": max_w0_local_gradient,
                "max_abs_lambda_rate": max_lambda,
                "w0_dynamic_over_mix": safe_ratio(max_w0_dynamic, max_mix),
                "w0_local_gradient_over_mix": safe_ratio(max_w0_local_gradient, max_mix),
                "lambda_over_mix": safe_ratio(max_lambda, max_mix),
                "max_abs_family_gap_rate": max_gap,
                "family_gap_over_mix": safe_ratio(max_gap, max_mix),
                "corr_mix_family_sum": pearson(mix_rate, family_sum),
                "dominant_family": dominant_family(
                    max_w0_dynamic, max_w0_local_gradient, max_lambda
                ),
            }
            rows.append(row)
            case_rows.append(row)

        if case_rows:
            top = max(case_rows, key=lambda row: row["max_abs_mix_rate"])
            summary_rows.append(
                {
                    "case_name": case_dir.name,
                    "top_channel": top["channel"],
                    "top_transport_status": top["transport_status"],
                    "top_ledger_status": top["ledger_status"],
                    "top_max_abs_mix_rate": top["max_abs_mix_rate"],
                    "top_w0_dynamic_over_mix": top["w0_dynamic_over_mix"],
                    "top_w0_local_gradient_over_mix": top["w0_local_gradient_over_mix"],
                    "top_lambda_over_mix": top["lambda_over_mix"],
                    "top_family_gap_over_mix": top["family_gap_over_mix"],
                    "top_dominant_family": top["dominant_family"],
                }
            )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_name",
        "channel",
        "transport_status",
        "ledger_status",
        "max_abs_net_rate",
        "max_abs_mix_rate",
        "max_abs_w0_dynamic_rate",
        "max_abs_w0_local_gradient_rate",
        "max_abs_lambda_rate",
        "w0_dynamic_over_mix",
        "w0_local_gradient_over_mix",
        "lambda_over_mix",
        "max_abs_family_gap_rate",
        "family_gap_over_mix",
        "corr_mix_family_sum",
        "dominant_family",
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
        "top_max_abs_mix_rate",
        "top_w0_dynamic_over_mix",
        "top_w0_local_gradient_over_mix",
        "top_lambda_over_mix",
        "top_family_gap_over_mix",
        "top_dominant_family",
    ]
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    out_md = Path(args.out_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Step-154 EM Pi Mixing Family Decomposition",
        "",
        "This report splits the mode-0 `pi` mixing path into:",
        "",
        "- `w0_dynamic`: explicit time-dependent `w0` dynamic mixing",
        "- `w0_local_gradient`: local-gradient `w0` mixing",
        "- `lambda`: dynamic `lambda` mixing",
        "",
        "The check is whether the total mixing diagnostic closes as",
        "`mix ≈ w0_dynamic + w0_local_gradient + lambda`, and which family dominates.",
        "",
        "## Case Summary",
        "",
        "| Case | Top Channel | Mix Max | w0 dynamic/mix | w0 grad/mix | lambda/mix | Gap/mix | Dominant |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| {case_name} | {top_channel} | {top_max_abs_mix_rate} | "
            "{top_w0_dynamic_over_mix} | {top_w0_local_gradient_over_mix} | "
            "{top_lambda_over_mix} | {top_family_gap_over_mix} | {top_dominant_family} |".format(
                case_name=row["case_name"],
                top_channel=row["top_channel"],
                top_max_abs_mix_rate=fmt(row["top_max_abs_mix_rate"]),
                top_w0_dynamic_over_mix=fmt(row["top_w0_dynamic_over_mix"]),
                top_w0_local_gradient_over_mix=fmt(row["top_w0_local_gradient_over_mix"]),
                top_lambda_over_mix=fmt(row["top_lambda_over_mix"]),
                top_family_gap_over_mix=fmt(row["top_family_gap_over_mix"]),
                top_dominant_family=row["top_dominant_family"],
            )
        )

    lines.extend(
        [
            "",
            "## Channel Detail",
            "",
            "| Case | Channel | Mix Max | w0 dynamic | w0 grad | lambda | dynamic/mix | grad/mix | lambda/mix | Gap/mix | Dominant | Corr(mix,sum) |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
        ]
    )
    for row in rows:
        lines.append(
            "| {case_name} | {channel} | {max_abs_mix_rate} | {max_abs_w0_dynamic_rate} | "
            "{max_abs_w0_local_gradient_rate} | {max_abs_lambda_rate} | "
            "{w0_dynamic_over_mix} | {w0_local_gradient_over_mix} | {lambda_over_mix} | "
            "{family_gap_over_mix} | {dominant_family} | {corr_mix_family_sum} |".format(
                case_name=row["case_name"],
                channel=row["channel"],
                max_abs_mix_rate=fmt(row["max_abs_mix_rate"]),
                max_abs_w0_dynamic_rate=fmt(row["max_abs_w0_dynamic_rate"]),
                max_abs_w0_local_gradient_rate=fmt(row["max_abs_w0_local_gradient_rate"]),
                max_abs_lambda_rate=fmt(row["max_abs_lambda_rate"]),
                w0_dynamic_over_mix=fmt(row["w0_dynamic_over_mix"]),
                w0_local_gradient_over_mix=fmt(row["w0_local_gradient_over_mix"]),
                lambda_over_mix=fmt(row["lambda_over_mix"]),
                family_gap_over_mix=fmt(row["family_gap_over_mix"]),
                dominant_family=row["dominant_family"],
                corr_mix_family_sum=fmt(row["corr_mix_family_sum"]),
            )
        )

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"detail_csv,{summary_csv}")


if __name__ == "__main__":
    main()
