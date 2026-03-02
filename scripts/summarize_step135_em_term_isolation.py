#!/usr/bin/env python3
"""Summarize Step-135 EM-term isolation vacuum runs."""

from __future__ import annotations

import argparse
import csv
import glob
import math
from pathlib import Path


def parse_hst(path: Path):
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
        raise RuntimeError(f"No history rows in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


def latest_hst(case_dir: Path) -> Path | None:
    for name in ("parthenon.out0.hst", "parthenon.out1.hst"):
        p = case_dir / name
        if p.exists():
            return p
    return None


def last(cols, name, default=math.nan):
    seq = cols.get(name)
    if not seq:
        return default
    return seq[-1]


def safe_ratio(a: float, b: float) -> float:
    if not math.isfinite(a) or not math.isfinite(b) or b == 0.0:
        return math.nan
    return a / b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", required=True, help="case directory glob")
    ap.add_argument("--baseline-case", default="term_baseline")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()

    case_dirs = sorted(Path(p) for p in glob.glob(args.glob))
    rows = []

    for case_dir in case_dirs:
        if not case_dir.is_dir():
            continue
        rc = math.nan
        rc_path = case_dir / "exit_code.txt"
        if rc_path.exists():
            try:
                rc = float(rc_path.read_text(encoding="utf-8").strip())
            except Exception:
                rc = math.nan

        hst = latest_hst(case_dir)
        if hst is None:
            rows.append(
                {
                    "case": case_dir.name,
                    "status": "FAIL_NO_HST",
                    "exit_code": rc,
                    "final_psi0_span": math.nan,
                    "final_psi_proj_span": math.nan,
                    "final_psi_w0_span": math.nan,
                    "em_u_bulk_initial": math.nan,
                    "em_u_bulk_final": math.nan,
                    "em_u_bulk_growth": math.nan,
                    "em_u_bulk_max": math.nan,
                    "final_rhs_a0_mass": math.nan,
                    "final_rhs_a0_damping": math.nan,
                    "final_rhs_a0_gauge": math.nan,
                    "final_rhs_a0_total": math.nan,
                    "final_src_em_mass_abs": math.nan,
                    "final_src_em_damping_abs": math.nan,
                    "final_src_em_gauge_abs": math.nan,
                    "notes": "missing history file",
                }
            )
            continue

        try:
            cols = parse_hst(hst)
            em = cols.get("m4d_em_u_bulk", [math.nan])
            row = {
                "case": case_dir.name,
                "status": "PASS_RUN" if rc == 0.0 else "FAIL_RUN",
                "exit_code": rc,
                "final_psi0_span": last(cols, "m4d_psi0_span"),
                "final_psi_proj_span": last(cols, "m4d_psi_proj_span"),
                "final_psi_w0_span": last(cols, "m4d_psi_w0_span"),
                "em_u_bulk_initial": em[0] if em else math.nan,
                "em_u_bulk_final": em[-1] if em else math.nan,
                "em_u_bulk_growth": safe_ratio(em[-1], em[0]) if em else math.nan,
                "em_u_bulk_max": max((abs(v) for v in em), default=math.nan),
                "final_rhs_a0_mass": last(cols, "m4d_int_rhs_a0_mode0_mass"),
                "final_rhs_a0_damping": last(cols, "m4d_int_rhs_a0_mode0_damping"),
                "final_rhs_a0_gauge": last(cols, "m4d_int_rhs_a0_mode0_gauge"),
                "final_rhs_a0_total": last(cols, "m4d_int_rhs_a0_mode0_total"),
                "final_src_em_mass_abs": last(cols, "m4d_int_src_em_mass_abs"),
                "final_src_em_damping_abs": last(cols, "m4d_int_src_em_damping_abs"),
                "final_src_em_gauge_abs": last(cols, "m4d_int_src_em_gauge_abs"),
                "notes": "",
            }
            rows.append(row)
        except Exception as exc:
            rows.append(
                {
                    "case": case_dir.name,
                    "status": "FAIL_PARSE",
                    "exit_code": rc,
                    "final_psi0_span": math.nan,
                    "final_psi_proj_span": math.nan,
                    "final_psi_w0_span": math.nan,
                    "em_u_bulk_initial": math.nan,
                    "em_u_bulk_final": math.nan,
                    "em_u_bulk_growth": math.nan,
                    "em_u_bulk_max": math.nan,
                    "final_rhs_a0_mass": math.nan,
                    "final_rhs_a0_damping": math.nan,
                    "final_rhs_a0_gauge": math.nan,
                    "final_rhs_a0_total": math.nan,
                    "final_src_em_mass_abs": math.nan,
                    "final_src_em_damping_abs": math.nan,
                    "final_src_em_gauge_abs": math.nan,
                    "notes": str(exc),
                }
            )

    rows.sort(key=lambda r: r["case"])
    baseline = next((r for r in rows if r["case"] == args.baseline_case), None)
    base_emmax = baseline["em_u_bulk_max"] if baseline else math.nan
    base_psi0 = baseline["final_psi0_span"] if baseline else math.nan
    for r in rows:
        r["emmax_vs_baseline"] = safe_ratio(r["em_u_bulk_max"], base_emmax)
        r["psi0_vs_baseline"] = safe_ratio(r["final_psi0_span"], base_psi0)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case",
        "status",
        "exit_code",
        "final_psi0_span",
        "final_psi_proj_span",
        "final_psi_w0_span",
        "em_u_bulk_initial",
        "em_u_bulk_final",
        "em_u_bulk_growth",
        "em_u_bulk_max",
        "final_rhs_a0_mass",
        "final_rhs_a0_damping",
        "final_rhs_a0_gauge",
        "final_rhs_a0_total",
        "final_src_em_mass_abs",
        "final_src_em_damping_abs",
        "final_src_em_gauge_abs",
        "emmax_vs_baseline",
        "psi0_vs_baseline",
        "notes",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out_md = Path(args.out_md)
    lines = [
        "# Step-135 EM-Term Isolation Vacuum Summary",
        "",
        f"- cases: `{len(rows)}`",
        f"- baseline_case: `{args.baseline_case}`",
        "",
        "| case | status | final psi0 | em_u_bulk growth | em_u_bulk max | rhs_a0(mass/damp/gauge/total) | emmax vs baseline |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        rhs = (
            f"{r['final_rhs_a0_mass']:.3e}/"
            f"{r['final_rhs_a0_damping']:.3e}/"
            f"{r['final_rhs_a0_gauge']:.3e}/"
            f"{r['final_rhs_a0_total']:.3e}"
        )
        lines.append(
            f"| {r['case']} | {r['status']} | {r['final_psi0_span']:.6e} | "
            f"{r['em_u_bulk_growth']:.6e} | {r['em_u_bulk_max']:.6e} | {rhs} | "
            f"{r['emmax_vs_baseline']:.6e} |"
        )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")


if __name__ == "__main__":
    main()
