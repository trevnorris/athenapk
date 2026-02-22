#!/usr/bin/env python3
"""Summarize Step-100 static-lambda sweep outputs."""

import argparse
import csv
import glob
import importlib.util
import math
import re
from pathlib import Path


CASE_RE = re.compile(r"lambda_(?P<lam>[0-9]+p[0-9]+)_t(?P<tlim>[0-9]+p[0-9]+)")
CASE_RESP_RE = re.compile(
    r"lambda_resp_g(?P<gain>[0-9]+(?:p[0-9]+)?)_l(?P<lam>[0-9]+p[0-9]+)_t(?P<tlim>[0-9]+p[0-9]+)"
)


def load_harris_scan_matrix_module(script_dir: Path):
    module_path = script_dir / "harris_scan_matrix.py"
    spec = importlib.util.spec_from_file_location("harris_scan_matrix", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_float_token(token: str):
    return float(token.replace("p", "."))


def parse_case_tokens(case_name: str):
    match = CASE_RE.match(case_name)
    if match:
        lam = parse_float_token(match.group("lam"))
        tlim = parse_float_token(match.group("tlim"))
        return lam, tlim, math.nan

    match = CASE_RESP_RE.match(case_name)
    if match:
        lam = parse_float_token(match.group("lam"))
        tlim = parse_float_token(match.group("tlim"))
        gain = parse_float_token(match.group("gain"))
        return lam, tlim, gain

    return math.nan, math.nan, math.nan


def ratio(numer, denom):
    if (
        numer is None
        or denom is None
        or not math.isfinite(numer)
        or not math.isfinite(denom)
        or abs(denom) == 0.0
    ):
        return math.nan
    return numer / denom


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        default="/projects/fluid-engine/out/step100_lambda_static_sweep_*/*",
        help="Glob for output case directories",
    )
    parser.add_argument(
        "--out-csv",
        default="/projects/fluid-engine/out/step100_lambda_static_sweep_summary.csv",
        help="Path to write summary CSV",
    )
    parser.add_argument(
        "--out-md",
        default="/projects/fluid-engine/out/step100_lambda_static_sweep_summary.md",
        help="Path to write summary markdown",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
    if not dirs:
        raise SystemExit(f"No directories matched: {args.glob}")

    rows = []
    for d in dirs:
        controlled_hst = d / "harris_controlled.out1.hst"
        full_hst = d / "harris_full.out1.hst"
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

        lambda_value, tlim, response_gain = parse_case_tokens(d.name)
        row = {
            "case_dir": str(d),
            "case_name": d.name,
            "lambda": lambda_value,
            "tlim": tlim,
            "response_gain": response_gain,
            "controlled_final_psi0_span": controlled.get("final_psi0_span", math.nan),
            "controlled_final_psi_proj_span": controlled.get("final_psi_proj_span", math.nan),
            "controlled_final_psi_w0_span": controlled.get("final_psi_w0_span", math.nan),
            "full_final_psi0_span": full.get("final_psi0_span", math.nan),
            "full_final_psi_proj_span": full.get("final_psi_proj_span", math.nan),
            "full_final_psi_proj_matched_span": full.get(
                "final_psi_proj_matched_span", math.nan
            ),
            "full_final_psi_proj_point_span": full.get("final_psi_proj_point_span", math.nan),
            "full_final_psi_proj_gaussian_span": full.get(
                "final_psi_proj_gaussian_span", math.nan
            ),
            "full_final_psi_w0_span": full.get("final_psi_w0_span", math.nan),
            "full_final_psi_w0_minus_psi_proj_span": full.get(
                "final_psi_w0_minus_psi_proj_span", math.nan
            ),
            "ratio_full_over_controlled_psi0": ratio(
                full.get("final_psi0_span"), controlled.get("final_psi0_span")
            ),
            "ratio_full_over_controlled_psi_proj": ratio(
                full.get("final_psi_proj_span"), controlled.get("final_psi_proj_span")
            ),
            "ratio_full_over_controlled_psi_w0": ratio(
                full.get("final_psi_w0_span"), controlled.get("final_psi_w0_span")
            ),
            "full_final_mixed_ew2": full.get("final_mixed_ew2", math.nan),
            "full_final_emf_vw_c_abs": full.get("final_emf_vw_c_abs", math.nan),
            "full_final_emf_cov_vxb_abs": full.get("final_emf_cov_vxb_abs", math.nan),
            "full_final_jw_ew": full.get("final_jw_ew", math.nan),
            "full_final_s_leak_abs": full.get("final_s_leak_abs", math.nan),
            "full_final_em_a2_even": full.get("final_em_a2_even", math.nan),
            "full_final_em_a2_odd": full.get("final_em_a2_odd", math.nan),
            "full_final_em_a2_odd_fraction": full.get(
                "final_em_a2_odd_fraction", math.nan
            ),
            "full_final_jw_l2_even": full.get("final_jw_l2_even", math.nan),
            "full_final_jw_l2_odd": full.get("final_jw_l2_odd", math.nan),
            "full_final_jw_l2_odd_fraction": full.get(
                "final_jw_l2_odd_fraction", math.nan
            ),
            "full_final_int_src_mode0_total_abs": full.get(
                "final_int_src_mode0_total_abs", math.nan
            ),
            "full_corr_psi0_emf_vw_c_abs": full.get("corr_psi0_emf_vw_c_abs", math.nan),
            "full_corr_psiw0_emf_vw_c_abs": full.get("corr_psiw0_emf_vw_c_abs", math.nan),
            "full_closure_status": full.get("closure_status", "N/A"),
            "full_closure_local_mode0_status": full.get(
                "closure_local_mode0_status", "N/A"
            ),
            "full_transport_closure_status": full.get("transport_closure_status", "N/A"),
            "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
            "full_em_bulk_ledger_best_label": full.get(
                "em_bulk_ledger_bridge_best_label", "N/A"
            ),
            "full_em_bulk_ledger_best_max_abs_rate": full.get(
                "em_bulk_ledger_bridge_best_max_abs_rate", math.nan
            ),
        }
        rows.append(row)

    rows.sort(
        key=lambda r: (
            (math.inf if not math.isfinite(r["lambda"]) else r["lambda"]),
            (math.inf if not math.isfinite(r["tlim"]) else r["tlim"]),
            (math.inf if not math.isfinite(r["response_gain"]) else r["response_gain"]),
            r["case_name"],
        )
    )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Lambda Sweep Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid cases found.\n")
        else:
            f.write(
                "| case | lambda | tlim | gain | psi0 ratio | psi_proj ratio | psi_w0 ratio | psi_w0-psi_proj | mixed_ew2 | emf_vw_c | emf_cov_vxb | src_mode0_total | closure_local | transport | ledger |\n"
            )
            f.write(
                "|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|\n"
            )
            for r in rows:
                f.write(
                    "| "
                    + " | ".join(
                        [
                            r["case_name"],
                            fmt(r["lambda"]),
                            fmt(r["tlim"]),
                            fmt(r["response_gain"]),
                            fmt(r["ratio_full_over_controlled_psi0"]),
                            fmt(r["ratio_full_over_controlled_psi_proj"]),
                            fmt(r["ratio_full_over_controlled_psi_w0"]),
                            fmt(r["full_final_psi_w0_minus_psi_proj_span"]),
                            fmt(r["full_final_mixed_ew2"]),
                            fmt(r["full_final_emf_vw_c_abs"]),
                            fmt(r["full_final_emf_cov_vxb_abs"]),
                            fmt(r["full_final_int_src_mode0_total_abs"]),
                            str(r["full_closure_local_mode0_status"]),
                            str(r["full_transport_closure_status"]),
                            str(r["full_em_bulk_ledger_status"]),
                        ]
                    )
                    + " |\n"
                )
            f.write("\n")
            f.write("## Ranking By |psi_w0 - psi_proj|\n\n")
            ranked = sorted(
                [r for r in rows if math.isfinite(r["full_final_psi_w0_minus_psi_proj_span"])],
                key=lambda r: abs(r["full_final_psi_w0_minus_psi_proj_span"]),
                reverse=True,
            )
            for idx, r in enumerate(ranked, start=1):
                f.write(
                    f"{idx}. {r['case_name']}: "
                    f"|psi_w0-psi_proj|={fmt(abs(r['full_final_psi_w0_minus_psi_proj_span']))}, "
                    f"psi0_ratio={fmt(r['ratio_full_over_controlled_psi0'])}, "
                    f"emf_vw_c={fmt(r['full_final_emf_vw_c_abs'])}, "
                    f"closure_local={r['full_closure_local_mode0_status']}, "
                    f"transport={r['full_transport_closure_status']}, "
                    f"ledger={r['full_em_bulk_ledger_status']}\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")


if __name__ == "__main__":
    main()
