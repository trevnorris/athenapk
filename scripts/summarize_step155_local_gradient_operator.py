#!/usr/bin/env python3
"""Summarize Step-155 local-gradient operator audit."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_harris_scan_matrix(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def max_val(values):
    if not values:
        return math.nan
    return max(values)


def final_val(values):
    if not values:
        return math.nan
    return values[-1]


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def safe_ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan
    return num / den


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Case directory glob")
    parser.add_argument("--baseline-case", default="baseline_local_grad_operator")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix(script_dir)
    case_dirs = [Path(path) for path in sorted(glob.glob(args.glob)) if Path(path).is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case dirs matched: {args.glob}")

    rows: list[dict[str, object]] = []
    for case_dir in case_dirs:
        run_log = case_dir / "run.log"
        hst = case_dir / "harris_full.out1.hst"
        if not run_log.exists() or not hst.exists():
            continue

        _, full = parse_case_rows(run_log)
        cols = hsm.parse_hst(hst)
        times = cols.get("time", [])

        mix_rate_t, mix_rate = derivative_series(times, cols.get("m4d_int_response_w0_local_gradient_mix_rate_abs"))
        mix_pi_t, mix_pi = derivative_series(times, cols.get("m4d_int_response_w0_local_gradient_mix_pi_abs"))
        del mix_rate_t, mix_pi_t

        max_mode1_grad = max_val(cols.get("m4d_response_w0_local_mode1_grad_abs", []))
        max_mode2_grad = max_val(cols.get("m4d_response_w0_local_mode2_grad_abs", []))
        max_dx_abs = max_val(cols.get("m4d_response_w0_local_dx_abs", []))
        max_dz_abs = max_val(cols.get("m4d_response_w0_local_dz_abs", []))
        dominant_mode = "mode1" if max_mode1_grad >= max_mode2_grad else "mode2"
        dominant_direction = "dx" if max_dx_abs >= max_dz_abs else "dz"

        pi0_transport = safe_float(full.get("pi0_transport_max_abs_rate", "nan"))
        piy_transport = safe_float(full.get("piy_transport_max_abs_rate", "nan"))
        piw_transport = safe_float(full.get("piw_transport_max_abs_rate", "nan"))
        transport_values = {
            "pi0": pi0_transport,
            "piy": piy_transport,
            "piw": piw_transport,
        }
        top_channel = max(transport_values, key=lambda key: transport_values[key])

        row = {
            "case_name": case_dir.name,
            "transport_status": full.get("transport_closure_status", ""),
            "ledger_status": full.get("em_bulk_ledger_status", ""),
            "pi0_transport_max_abs_rate": pi0_transport,
            "piy_transport_max_abs_rate": piy_transport,
            "piw_transport_max_abs_rate": piw_transport,
            "top_pi_channel": top_channel,
            "top_pi_transport_max_abs_rate": transport_values[top_channel],
            "max_response_w0_local_abs": max_val(cols.get("m4d_response_w0_local_abs", [])),
            "max_response_w0_local_mode1_abs": max_val(cols.get("m4d_response_w0_local_mode1_abs", [])),
            "max_response_w0_local_mode2_abs": max_val(cols.get("m4d_response_w0_local_mode2_abs", [])),
            "max_response_w0_local_dx_abs": max_dx_abs,
            "max_response_w0_local_dz_abs": max_dz_abs,
            "max_response_w0_local_mode1_dx_abs": max_val(cols.get("m4d_response_w0_local_mode1_dx_abs", [])),
            "max_response_w0_local_mode1_dz_abs": max_val(cols.get("m4d_response_w0_local_mode1_dz_abs", [])),
            "max_response_w0_local_mode1_grad_abs": max_mode1_grad,
            "max_response_w0_local_mode2_dx_abs": max_val(cols.get("m4d_response_w0_local_mode2_dx_abs", [])),
            "max_response_w0_local_mode2_dz_abs": max_val(cols.get("m4d_response_w0_local_mode2_dz_abs", [])),
            "max_response_w0_local_mode2_grad_abs": max_mode2_grad,
            "max_response_w0_local_grad_abs": max_val(cols.get("m4d_response_w0_local_grad_abs", [])),
            "dominant_local_mode": dominant_mode,
            "dominant_gradient_direction": dominant_direction,
            "max_local_gradient_mix_rate_abs": max_abs(mix_rate),
            "max_local_gradient_mix_pi_abs_rate": max_abs(mix_pi),
            "final_local_gradient_mix_rate_abs": final_val(cols.get("m4d_int_response_w0_local_gradient_mix_rate_abs", [])),
            "final_local_gradient_mix_pi_abs": final_val(cols.get("m4d_int_response_w0_local_gradient_mix_pi_abs", [])),
        }
        row["final_mix_pi_over_mix_rate"] = safe_ratio(
            float(row["final_local_gradient_mix_pi_abs"]),
            float(row["final_local_gradient_mix_rate_abs"]),
        )
        rows.append(row)

    if not rows:
        raise SystemExit("No valid case rows found")

    baseline = next((row for row in rows if row["case_name"] == args.baseline_case), None)
    if baseline is not None:
        for row in rows:
            row["top_pi_transport_over_baseline"] = safe_ratio(
                float(row["top_pi_transport_max_abs_rate"]),
                float(baseline["top_pi_transport_max_abs_rate"]),
            )
            row["mix_rate_over_baseline"] = safe_ratio(
                float(row["max_local_gradient_mix_rate_abs"]),
                float(baseline["max_local_gradient_mix_rate_abs"]),
            )
            row["mix_pi_over_baseline"] = safe_ratio(
                float(row["max_local_gradient_mix_pi_abs_rate"]),
                float(baseline["max_local_gradient_mix_pi_abs_rate"]),
            )
    else:
        for row in rows:
            row["top_pi_transport_over_baseline"] = math.nan
            row["mix_rate_over_baseline"] = math.nan
            row["mix_pi_over_baseline"] = math.nan

    fieldnames = [
        "case_name",
        "transport_status",
        "ledger_status",
        "pi0_transport_max_abs_rate",
        "piy_transport_max_abs_rate",
        "piw_transport_max_abs_rate",
        "top_pi_channel",
        "top_pi_transport_max_abs_rate",
        "top_pi_transport_over_baseline",
        "max_response_w0_local_abs",
        "max_response_w0_local_mode1_abs",
        "max_response_w0_local_mode2_abs",
        "max_response_w0_local_dx_abs",
        "max_response_w0_local_dz_abs",
        "max_response_w0_local_mode1_dx_abs",
        "max_response_w0_local_mode1_dz_abs",
        "max_response_w0_local_mode1_grad_abs",
        "max_response_w0_local_mode2_dx_abs",
        "max_response_w0_local_mode2_dz_abs",
        "max_response_w0_local_mode2_grad_abs",
        "max_response_w0_local_grad_abs",
        "dominant_local_mode",
        "dominant_gradient_direction",
        "max_local_gradient_mix_rate_abs",
        "max_local_gradient_mix_pi_abs_rate",
        "mix_rate_over_baseline",
        "mix_pi_over_baseline",
        "final_local_gradient_mix_rate_abs",
        "final_local_gradient_mix_pi_abs",
        "final_mix_pi_over_mix_rate",
    ]

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    lines = [
        "# Step-155 Local-Gradient Operator Audit",
        "",
        "This audit decomposes the `w0_local_gradient` operator into its local-warp mode content and directional content on the focused Step-154 Harris baseline.",
        "",
        "## Case Summary",
        "",
        "| Case | Transport | Ledger | Top pi channel | Top pi max abs rate | Top pi / baseline | Dominant local mode | Dominant dir | Max grad abs | Max mix-rate abs | Max mix-pi abs rate | Mix-pi / mix-rate (final) |",
        "| --- | --- | --- | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {case_name} | {transport_status} | {ledger_status} | {top_pi_channel} | {top_pi_transport_max_abs_rate} | {top_pi_transport_over_baseline} | {dominant_local_mode} | {dominant_gradient_direction} | {max_response_w0_local_grad_abs} | {max_local_gradient_mix_rate_abs} | {max_local_gradient_mix_pi_abs_rate} | {final_mix_pi_over_mix_rate} |".format(
                case_name=row["case_name"],
                transport_status=row["transport_status"],
                ledger_status=row["ledger_status"],
                top_pi_channel=row["top_pi_channel"],
                top_pi_transport_max_abs_rate=fmt(float(row["top_pi_transport_max_abs_rate"])),
                top_pi_transport_over_baseline=fmt(float(row["top_pi_transport_over_baseline"])),
                dominant_local_mode=row["dominant_local_mode"],
                dominant_gradient_direction=row["dominant_gradient_direction"],
                max_response_w0_local_grad_abs=fmt(float(row["max_response_w0_local_grad_abs"])),
                max_local_gradient_mix_rate_abs=fmt(float(row["max_local_gradient_mix_rate_abs"])),
                max_local_gradient_mix_pi_abs_rate=fmt(float(row["max_local_gradient_mix_pi_abs_rate"])),
                final_mix_pi_over_mix_rate=fmt(float(row["final_mix_pi_over_mix_rate"])),
            )
        )

    lines.extend([
        "",
        "## Local-Warp Content",
        "",
        "| Case | mode1 abs | mode2 abs | mode1 grad abs | mode2 grad abs | dx abs | dz abs |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in rows:
        lines.append(
            "| {case_name} | {max_response_w0_local_mode1_abs} | {max_response_w0_local_mode2_abs} | {max_response_w0_local_mode1_grad_abs} | {max_response_w0_local_mode2_grad_abs} | {max_response_w0_local_dx_abs} | {max_response_w0_local_dz_abs} |".format(
                case_name=row["case_name"],
                max_response_w0_local_mode1_abs=fmt(float(row["max_response_w0_local_mode1_abs"])),
                max_response_w0_local_mode2_abs=fmt(float(row["max_response_w0_local_mode2_abs"])),
                max_response_w0_local_mode1_grad_abs=fmt(float(row["max_response_w0_local_mode1_grad_abs"])),
                max_response_w0_local_mode2_grad_abs=fmt(float(row["max_response_w0_local_mode2_grad_abs"])),
                max_response_w0_local_dx_abs=fmt(float(row["max_response_w0_local_dx_abs"])),
                max_response_w0_local_dz_abs=fmt(float(row["max_response_w0_local_dz_abs"])),
            )
        )

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
