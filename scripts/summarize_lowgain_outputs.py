#!/usr/bin/env python3
"""Summarize low-gain Harris scan outputs from precomputed .hst files."""

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


def parse_gain_tag(tag: str):
    sign = 1.0
    if tag.startswith("m"):
        sign = -1.0
        tag = tag[1:]
    elif tag.startswith("minus"):
        sign = -1.0
        tag = tag[5:]

    # Examples:
    #   "0" -> 0
    #   "1" -> 1
    #   "1em4" -> 1e-4
    #   "5em2" -> 5e-2
    if "em" in tag:
        base, exp = tag.split("em", 1)
        if not base or not exp:
            return math.nan
        try:
            return sign * float(base) * (10.0 ** (-int(exp)))
        except ValueError:
            return math.nan
    try:
        return sign * float(tag)
    except ValueError:
        return math.nan


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
        default="/projects/fluid-engine/out/step63_lowgain_g*",
        help="Glob for output case directories",
    )
    parser.add_argument(
        "--out-csv",
        default="/projects/fluid-engine/out/step63_lowgain_summary.csv",
        help="Path to write summary CSV",
    )
    parser.add_argument(
        "--out-md",
        default="/projects/fluid-engine/out/step63_lowgain_summary.md",
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
            1.0,
            1.0,
            1.0e-6,
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
            1.0,
            1.0,
            1.0e-6,
            1.0e-12,
            1.0e-3,
            1.0e-12,
        )

        tag = d.name.split("step63_lowgain_g", 1)[-1]
        gain = parse_gain_tag(tag)
        baseline = full.get("em_bulk_ledger_max_abs_rate", math.nan)
        plus = full.get("em_bulk_ledger_with_bridge_plus_max_abs_rate", math.nan)
        minus = full.get("em_bulk_ledger_with_bridge_minus_max_abs_rate", math.nan)
        best = full.get("em_bulk_ledger_bridge_best_max_abs_rate", math.nan)

        row = {
            "case_dir": str(d),
            "gain_tag": tag,
            "gain": gain,
            "ratio_full_over_controlled_psi0": ratio(
                full.get("final_psi0_span"), controlled.get("final_psi0_span")
            ),
            "ratio_full_over_controlled_psi_proj": ratio(
                full.get("final_psi_proj_span"), controlled.get("final_psi_proj_span")
            ),
            "ratio_full_over_controlled_psi_w0": ratio(
                full.get("final_psi_w0_span"), controlled.get("final_psi_w0_span")
            ),
            "full_final_psi0_span": full.get("final_psi0_span", math.nan),
            "full_final_psi_proj_span": full.get("final_psi_proj_span", math.nan),
            "full_final_psi_w0_span": full.get("final_psi_w0_span", math.nan),
            "full_final_int_rhs_ay_mode0_even_bridge": full.get(
                "final_int_rhs_ay_mode0_even_bridge", math.nan
            ),
            "full_final_int_rhs_ay_mode0_even_bridge_abs": full.get(
                "final_int_rhs_ay_mode0_even_bridge_abs", math.nan
            ),
            "full_final_int_bridge_power_mode0": full.get(
                "final_int_bridge_power_mode0", math.nan
            ),
            "full_final_int_bridge_power_mode0_abs": full.get(
                "final_int_bridge_power_mode0_abs", math.nan
            ),
            "full_closure_status": full.get("closure_status", "N/A"),
            "full_closure_local_mode0_status": full.get(
                "closure_local_mode0_status", "N/A"
            ),
            "full_transport_closure_status": full.get(
                "transport_closure_status", "N/A"
            ),
            "full_em_bulk_ledger_status": full.get("em_bulk_ledger_status", "N/A"),
            "em_bulk_ledger_baseline_max_abs_rate": baseline,
            "em_bulk_ledger_with_bridge_plus_max_abs_rate": plus,
            "em_bulk_ledger_with_bridge_minus_max_abs_rate": minus,
            "em_bulk_ledger_bridge_best_label": full.get(
                "em_bulk_ledger_bridge_best_label", "N/A"
            ),
            "em_bulk_ledger_bridge_best_max_abs_rate": best,
            "em_bulk_ledger_improvement_vs_baseline": (
                baseline / best
                if math.isfinite(baseline)
                and math.isfinite(best)
                and best > 0.0
                else math.nan
            ),
        }
        rows.append(row)

    rows.sort(
        key=lambda r: (
            (math.inf if not math.isfinite(r["gain"]) else r["gain"]),
            r["gain_tag"],
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
        f.write("# Step63 Low-Gain Summary\n\n")
        f.write(f"Parsed `{len(rows)}` cases from glob `{args.glob}`.\n\n")
        if not rows:
            f.write("No valid cases found.\n")
        else:
            f.write(
                "| gain | psi0 ratio | psi_proj ratio | psi_w0 ratio | ledger baseline | ledger best | best label | improvement |\n"
            )
            f.write(
                "|---:|---:|---:|---:|---:|---:|:---|---:|\n"
            )
            for r in rows:
                f.write(
                    "| "
                    + " | ".join(
                        [
                            fmt(r["gain"]),
                            fmt(r["ratio_full_over_controlled_psi0"]),
                            fmt(r["ratio_full_over_controlled_psi_proj"]),
                            fmt(r["ratio_full_over_controlled_psi_w0"]),
                            fmt(r["em_bulk_ledger_baseline_max_abs_rate"]),
                            fmt(r["em_bulk_ledger_bridge_best_max_abs_rate"]),
                            str(r["em_bulk_ledger_bridge_best_label"]),
                            fmt(r["em_bulk_ledger_improvement_vs_baseline"]),
                        ]
                    )
                    + " |\n"
                )
            f.write("\n")
            f.write(
                "## Candidate Ranking (lowest best-ledger residual)\n\n"
            )
            ranked = sorted(
                [r for r in rows if math.isfinite(r["em_bulk_ledger_bridge_best_max_abs_rate"])],
                key=lambda r: r["em_bulk_ledger_bridge_best_max_abs_rate"],
            )
            if not ranked:
                f.write("No finite ledger residuals found.\n")
            else:
                for idx, r in enumerate(ranked, start=1):
                    f.write(
                        f"{idx}. gain={fmt(r['gain'])}, "
                        f"best={fmt(r['em_bulk_ledger_bridge_best_max_abs_rate'])}, "
                        f"label={r['em_bulk_ledger_bridge_best_label']}, "
                        f"psi0_ratio={fmt(r['ratio_full_over_controlled_psi0'])}, "
                        f"closure_local={r['full_closure_local_mode0_status']}, "
                        f"transport={r['full_transport_closure_status']}\n"
                    )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")
    if rows:
        best_row = min(
            rows,
            key=lambda r: (
                math.inf
                if not math.isfinite(r["em_bulk_ledger_bridge_best_max_abs_rate"])
                else r["em_bulk_ledger_bridge_best_max_abs_rate"]
            ),
        )
        print(f"best_gain,{fmt(best_row['gain'])}")
        print(
            "best_ledger_max_abs_rate,"
            f"{fmt(best_row['em_bulk_ledger_bridge_best_max_abs_rate'])}"
        )
        print(
            "best_psi0_ratio,"
            f"{fmt(best_row['ratio_full_over_controlled_psi0'])}"
        )


if __name__ == "__main__":
    main()
