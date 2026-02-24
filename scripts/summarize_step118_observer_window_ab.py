#!/usr/bin/env python3
"""Summarize Step-118 observer-window A/B runs."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
from pathlib import Path


def load_harris_scan_matrix_module(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def series_max_abs(series):
    if series is None or len(series) == 0:
        return math.nan
    finite = [abs(float(v)) for v in series if math.isfinite(float(v))]
    if not finite:
        return math.nan
    return max(finite)


def ratio(num, den):
    if not (math.isfinite(num) and math.isfinite(den)):
        return math.nan
    if den == 0.0:
        return math.nan
    return num / den


def parse_arg_from_log(log_path: Path, key: str):
    if not log_path.exists():
        return math.nan
    token = f"{key}="
    last_value = None
    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        for _ in range(96):
            line = f.readline()
            if not line:
                break
            if not line.startswith("command,"):
                continue
            start = 0
            while True:
                idx = line.find(token, start)
                if idx < 0:
                    break
                tail = line[idx + len(token) :].strip()
                last_value = tail.split()[0].split(",")[0]
                start = idx + len(token)
            break
    if last_value is None:
        return math.nan
    low = str(last_value).lower()
    if low in {"true", "false"}:
        return low
    try:
        return float(last_value)
    except ValueError:
        return last_value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", required=True, help="Glob for case directories")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--out-md", required=True, help="Output markdown path")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    rows = []
    for case_dir in sorted(Path(p) for p in glob.glob(args.glob)):
        controlled_hst = case_dir / "harris_controlled.out1.hst"
        full_hst = case_dir / "harris_full.out1.hst"
        if not controlled_hst.exists() or not full_hst.exists():
            continue

        controlled_cols = hsm.parse_hst(controlled_hst)
        full_cols = hsm.parse_hst(full_hst)
        controlled = hsm.analyze_case(
            "controlled",
            controlled_cols,
            5.0e-2,
            1.0e-8,
            1.0e-8,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )
        full = hsm.analyze_case(
            "full",
            full_cols,
            5.0e-2,
            1.0e-8,
            1.0e-8,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )
        log_path = case_dir / "run.log"

        psi0_ctrl = controlled.get("final_psi0_span", math.nan)
        psi0_full = full.get("final_psi0_span", math.nan)
        psiw0_minus_proj = full.get("final_psi_w0_minus_psi_proj_span", math.nan)
        jw_ew = full.get("final_jw_ew", math.nan)
        s_leak_abs = full.get("final_s_leak_abs", math.nan)
        em_leak_w = full.get("final_em_leak_w", math.nan)
        emf_vw_c_abs = full.get("final_emf_vw_c_abs", math.nan)

        rows.append(
            {
                "case_name": case_dir.name,
                "projection_enable": parse_arg_from_log(
                    log_path, "modes4d/response_w0_projection_enable"
                ),
                "geometry_shift_enable": parse_arg_from_log(
                    log_path, "modes4d/response_w0_geometry_shift_enable"
                ),
                "projection_gain": parse_arg_from_log(
                    log_path, "modes4d/response_w0_projection_gain"
                ),
                "projection_max_abs": parse_arg_from_log(
                    log_path, "modes4d/response_w0_projection_max_abs"
                ),
                "full_final_projection_center_w": full.get(
                    "final_projection_center_w", math.nan
                ),
                "ratio_full_over_controlled_psi0": ratio(psi0_full, psi0_ctrl),
                "full_final_psi_w0_minus_psi_proj_span": psiw0_minus_proj,
                "full_final_jw_ew": jw_ew,
                "full_max_abs_jw_ew": series_max_abs(full_cols.get("m4d_int_jw_ew")),
                "full_final_s_leak_abs": s_leak_abs,
                "full_max_abs_s_leak_abs": series_max_abs(full_cols.get("m4d_int_s_leak_abs")),
                "full_final_em_leak_w": em_leak_w,
                "full_max_abs_em_leak_w": series_max_abs(full_cols.get("m4d_em_leak_w")),
                "full_final_emf_vw_c_abs": emf_vw_c_abs,
                "full_max_abs_emf_vw_c_abs": series_max_abs(full_cols.get("m4d_emf_vw_c_abs")),
                "full_closure_projected_status": full.get("closure_projected_status", "N/A"),
                "full_closure_local_mode0_status": full.get(
                    "closure_local_mode0_status", "N/A"
                ),
                "full_transport_closure_status": full.get("transport_closure_status", "N/A"),
                "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
            }
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        else:
            f.write("")

    fixed = None
    comoving = None
    for r in rows:
        if str(r["projection_enable"]).lower() == "false":
            fixed = r
        elif str(r["projection_enable"]).lower() == "true":
            comoving = r

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-118 Observer Window A/B Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid case outputs found.\n")
        else:
            f.write(
                "| case | projection enable/gain/max | projection_center_w | psi0 ratio | psi_w0-psi_proj | JwEw(final/max) | S_leak_abs(final/max) | S_EMw(final/max) | emf_vw_c(final/max) | projected/local/transport/ledger |\n"
            )
            f.write(
                "|:---|:---|---:|---:|---:|---:|---:|---:|---:|:---|\n"
            )
            for r in rows:
                status = (
                    f"{r['full_closure_projected_status']}/"
                    f"{r['full_closure_local_mode0_status']}/"
                    f"{r['full_transport_closure_status']}/"
                    f"{r['full_em_bulk_ledger_status']}"
                )
                f.write(
                    "| "
                    + " | ".join(
                        [
                            r["case_name"],
                            f"{r['projection_enable']}/{fmt(r['projection_gain'])}/{fmt(r['projection_max_abs'])}",
                            fmt(r["full_final_projection_center_w"]),
                            fmt(r["ratio_full_over_controlled_psi0"]),
                            fmt(r["full_final_psi_w0_minus_psi_proj_span"]),
                            f"{fmt(r['full_final_jw_ew'])}/{fmt(r['full_max_abs_jw_ew'])}",
                            f"{fmt(r['full_final_s_leak_abs'])}/{fmt(r['full_max_abs_s_leak_abs'])}",
                            f"{fmt(r['full_final_em_leak_w'])}/{fmt(r['full_max_abs_em_leak_w'])}",
                            f"{fmt(r['full_final_emf_vw_c_abs'])}/{fmt(r['full_max_abs_emf_vw_c_abs'])}",
                            status,
                        ]
                    )
                    + " |\n"
                )

            if fixed is not None and comoving is not None:
                f.write("\n## A/B Delta (comoving - fixed)\n\n")
                f.write(
                    "| metric | delta |\n|:---|---:|\n"
                )
                for key in (
                    "full_final_projection_center_w",
                    "full_final_psi_w0_minus_psi_proj_span",
                    "full_final_jw_ew",
                    "full_max_abs_jw_ew",
                    "full_final_s_leak_abs",
                    "full_max_abs_s_leak_abs",
                    "full_final_em_leak_w",
                    "full_max_abs_em_leak_w",
                    "full_final_emf_vw_c_abs",
                    "full_max_abs_emf_vw_c_abs",
                ):
                    a = comoving.get(key, math.nan)
                    b = fixed.get(key, math.nan)
                    delta = math.nan
                    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                        if math.isfinite(a) and math.isfinite(b):
                            delta = a - b
                    f.write(f"| {key} | {fmt(delta)} |\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
