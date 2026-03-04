#!/usr/bin/env python3
"""Summarize Step-140 mode-0 closure / ledger term-isolation sweep."""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import math
import re
from pathlib import Path
from typing import Dict, List


def load_harris_scan_matrix_module(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def series_max_abs(series):
    if series is None or len(series) == 0:
        return math.nan
    vals = [abs(float(v)) for v in series if math.isfinite(float(v))]
    if not vals:
        return math.nan
    return max(vals)


def safe_ratio(numer: float, denom: float) -> float:
    if not (math.isfinite(numer) and math.isfinite(denom)):
        return math.nan
    if abs(denom) == 0.0:
        return math.nan
    return numer / denom


def parse_arg_from_log(log_path: Path, key: str) -> float:
    if not log_path.exists():
        return math.nan
    pattern = re.compile(rf"{re.escape(key)}=([^\s,]+)")
    try:
        text = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return math.nan
    matches = list(pattern.finditer(text))
    if not matches:
        return math.nan
    # Use the last match so explicit case-specific overrides win over base args.
    raw = matches[-1].group(1)
    low = raw.strip().lower()
    if low in {"true", "on", "yes"}:
        return 1.0
    if low in {"false", "off", "no"}:
        return 0.0
    try:
        return float(raw)
    except ValueError:
        return math.nan
    return math.nan


def fmt(v) -> str:
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "nan"
    return f"{v:.6e}"


def status_floor(value: float, floor: float) -> str:
    if not math.isfinite(value):
        return "FAIL"
    return "PASS" if value >= floor else "FAIL"


def best_row(rows: List[Dict[str, object]], key: str, require_active: bool = False):
    cand = []
    for row in rows:
        if require_active and row.get("active_transfer_status") != "PASS":
            continue
        value = row.get(key, math.nan)
        if isinstance(value, (int, float)) and math.isfinite(float(value)):
            cand.append(row)
    if not cand:
        return None
    return min(cand, key=lambda r: float(r[key]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Case directory glob (e.g., /path/out/*)")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--out-md", required=True, help="Output markdown path")
    parser.add_argument("--min-max-jw-ew", type=float, default=1.0e-12)
    parser.add_argument("--min-max-s-leak-abs", type=float, default=1.0e-12)
    parser.add_argument("--min-max-emf-vw-c-abs", type=float, default=1.0e-4)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    case_dirs = [Path(p) for p in sorted(glob.glob(args.glob)) if Path(p).is_dir()]
    rows: List[Dict[str, object]] = []

    for case_dir in case_dirs:
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

        max_abs_jw_ew = series_max_abs(full_cols.get("m4d_int_jw_ew"))
        max_abs_s_leak_abs = series_max_abs(full_cols.get("m4d_int_s_leak_abs"))
        max_abs_emf_vw_c_abs = series_max_abs(full_cols.get("m4d_emf_vw_c_abs"))

        gate_jw_ew = status_floor(max_abs_jw_ew, args.min_max_jw_ew)
        gate_s_leak = status_floor(max_abs_s_leak_abs, args.min_max_s_leak_abs)
        gate_emf = status_floor(max_abs_emf_vw_c_abs, args.min_max_emf_vw_c_abs)
        active = "PASS" if gate_jw_ew == "PASS" and gate_s_leak == "PASS" and gate_emf == "PASS" else "FAIL"

        log_path = case_dir / "run.log"
        row: Dict[str, object] = {
            "case_name": case_dir.name,
            "mass_matrix_projection_enable": parse_arg_from_log(
                log_path, "modes4d/mass_matrix_projection_enable"
            ),
            "em_source_mass_gain": parse_arg_from_log(log_path, "modes4d/em_source_mass_gain"),
            "em_source_laplacian_gain": parse_arg_from_log(log_path, "modes4d/em_source_laplacian_gain"),
            "em_source_current_gain": parse_arg_from_log(log_path, "modes4d/em_source_current_gain"),
            "em_source_timelike_gain": parse_arg_from_log(log_path, "modes4d/em_source_timelike_gain"),
            "em_source_damping_gain": parse_arg_from_log(log_path, "modes4d/em_source_damping_gain"),
            "plasma_momw_source_gain": parse_arg_from_log(log_path, "modes4d/plasma_momw_source_gain"),
            "plasma_momw_pressure_source_gain": parse_arg_from_log(
                log_path, "modes4d/plasma_momw_pressure_source_gain"
            ),
            "controlled_final_psi0_span": controlled.get("final_psi0_span", math.nan),
            "full_final_psi0_span": full.get("final_psi0_span", math.nan),
            "ratio_full_over_controlled_psi0": safe_ratio(
                float(full.get("final_psi0_span", math.nan)),
                float(controlled.get("final_psi0_span", math.nan)),
            ),
            "full_max_abs_jw_ew": max_abs_jw_ew,
            "full_max_abs_s_leak_abs": max_abs_s_leak_abs,
            "full_max_abs_emf_vw_c_abs": max_abs_emf_vw_c_abs,
            "gate_jw_ew": gate_jw_ew,
            "gate_s_leak_abs": gate_s_leak,
            "gate_emf_vw_c_abs": gate_emf,
            "active_transfer_status": active,
            "full_final_int_src_mode0_total_abs": full.get("final_int_src_mode0_total_abs", math.nan),
            "full_final_int_src_mode0_em_abs": full.get("final_int_src_mode0_em_abs", math.nan),
            "full_final_int_src_mode0_timelike_abs": full.get("final_int_src_mode0_timelike_abs", math.nan),
            "full_closure_projected_status": full.get("closure_projected_status", "N/A"),
            "full_closure_local_mode0_status": full.get("closure_local_mode0_status", "N/A"),
            "full_transport_closure_status": full.get("transport_closure_status", "N/A"),
            "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
            "full_closure_local_mode0_max_abs_rate": full.get(
                "closure_local_mode0_max_abs_rate", math.nan
            ),
            "full_closure_local_mode0_final_max_abs_rate": full.get(
                "closure_local_mode0_final_max_abs_rate", math.nan
            ),
            "full_closure_local_mode0_final_l1": full.get("closure_local_mode0_final_l1", math.nan),
            "full_em_bulk_ledger_max_abs_rate": full.get("em_bulk_ledger_max_abs_rate", math.nan),
            "full_em_bulk_ledger_rms_abs_rate": full.get("em_bulk_ledger_rms_abs_rate", math.nan),
            "full_em_bulk_ledger_final_rate": full.get("em_bulk_ledger_final_rate", math.nan),
        }
        rows.append(row)

    rows.sort(key=lambda r: str(r["case_name"]))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["case_name"])
        writer.writeheader()
        if rows:
            writer.writerows(rows)

    best_local = best_row(rows, "full_closure_local_mode0_max_abs_rate", require_active=True)
    best_ledger = best_row(rows, "full_em_bulk_ledger_max_abs_rate", require_active=True)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Step-140 Mode-0 Closure/Ledger Term-Isolation Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        f.write("Activation floors:\n")
        f.write(f"- `max|JwEw| >= {args.min_max_jw_ew:.3e}`\n")
        f.write(f"- `max|S_leak_abs| >= {args.min_max_s_leak_abs:.3e}`\n")
        f.write(f"- `max|emf_vw_c_abs| >= {args.min_max_emf_vw_c_abs:.3e}`\n\n")

        if rows:
            f.write(
                "| case | mm_proj | mass/lap/current/time/damp | momw/press | active | max|JwEw| | max|S_leak| | max|EMF| | local status | local max rate | local final l1 | transport | ledger | ledger max rate | src_mode0 (total/em/time) | psi0 ratio |\n"
            )
            f.write(
                "|:---|---:|---:|---:|:---|---:|---:|---:|:---|---:|---:|:---|:---|---:|---:|---:|\n"
            )
            for r in rows:
                gain_pack = (
                    f"{fmt(r['em_source_mass_gain'])}/"
                    f"{fmt(r['em_source_laplacian_gain'])}/"
                    f"{fmt(r['em_source_current_gain'])}/"
                    f"{fmt(r['em_source_timelike_gain'])}/"
                    f"{fmt(r['em_source_damping_gain'])}"
                )
                momw_pack = (
                    f"{fmt(r['plasma_momw_source_gain'])}/"
                    f"{fmt(r['plasma_momw_pressure_source_gain'])}"
                )
                src_pack = (
                    f"{fmt(r['full_final_int_src_mode0_total_abs'])}/"
                    f"{fmt(r['full_final_int_src_mode0_em_abs'])}/"
                    f"{fmt(r['full_final_int_src_mode0_timelike_abs'])}"
                )
                f.write(
                    "| "
                    + " | ".join(
                        [
                            str(r["case_name"]),
                            fmt(r["mass_matrix_projection_enable"]),
                            gain_pack,
                            momw_pack,
                            str(r["active_transfer_status"]),
                            fmt(r["full_max_abs_jw_ew"]),
                            fmt(r["full_max_abs_s_leak_abs"]),
                            fmt(r["full_max_abs_emf_vw_c_abs"]),
                            str(r["full_closure_local_mode0_status"]),
                            fmt(r["full_closure_local_mode0_max_abs_rate"]),
                            fmt(r["full_closure_local_mode0_final_l1"]),
                            str(r["full_transport_closure_status"]),
                            str(r["full_em_bulk_ledger_status"]),
                            fmt(r["full_em_bulk_ledger_max_abs_rate"]),
                            src_pack,
                            fmt(r["ratio_full_over_controlled_psi0"]),
                        ]
                    )
                    + " |\n"
                )

            f.write("\n")
            if best_local is not None:
                f.write(
                    "best_active_local_mode0_case,"
                    f"{best_local['case_name']},"
                    f"max_abs_rate={fmt(best_local['full_closure_local_mode0_max_abs_rate'])},"
                    f"status={best_local['full_closure_local_mode0_status']}\n"
                )
            else:
                f.write("best_active_local_mode0_case,none\n")

            if best_ledger is not None:
                f.write(
                    "best_active_ledger_case,"
                    f"{best_ledger['case_name']},"
                    f"max_abs_rate={fmt(best_ledger['full_em_bulk_ledger_max_abs_rate'])},"
                    f"status={best_ledger['full_em_bulk_ledger_status']}\n"
                )
            else:
                f.write("best_active_ledger_case,none\n")
        else:
            f.write("No valid case outputs found.\n")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
