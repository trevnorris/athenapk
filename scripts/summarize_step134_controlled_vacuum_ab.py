#!/usr/bin/env python3
"""Summarize Step-134 controlled vacuum A/B runs."""

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
        parts = line.split()
        rows.append([float(x) for x in parts])

    if not rows:
        raise RuntimeError(f"No history rows in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


def first_existing(cols, names, default=math.nan):
    for name in names:
        if name in cols:
            return cols[name]
    return [default]


def latest_hst(case_dir: Path) -> Path | None:
    for name in ("parthenon.out1.hst", "parthenon.out0.hst"):
        p = case_dir / name
        if p.exists():
            return p
    return None


def safe_ratio(a: float, b: float) -> float:
    if not math.isfinite(a) or not math.isfinite(b) or b == 0.0:
        return math.nan
    return a / b


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Case directory glob")
    p.add_argument("--baseline-case", default="vacuum_nohard")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-md", required=True)
    args = p.parse_args()

    case_dirs = sorted(Path(p) for p in glob.glob(args.glob))
    rows = []

    for case_dir in case_dirs:
        if not case_dir.is_dir():
            continue
        case_name = case_dir.name
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
                    "case": case_name,
                    "status": "FAIL_NO_HST",
                    "exit_code": rc,
                    "final_psi0_span": math.nan,
                    "final_psi_proj_span": math.nan,
                    "final_psi_w0_span": math.nan,
                    "em_u_bulk_initial": math.nan,
                    "em_u_bulk_final": math.nan,
                    "em_u_bulk_growth": math.nan,
                    "em_u_bulk_max": math.nan,
                    "psi0_max_abs": math.nan,
                    "notes": "missing history output",
                }
            )
            continue

        try:
            cols = parse_hst(hst)
            psi0 = first_existing(cols, ["m4d_psi0_span"])
            psiproj = first_existing(cols, ["m4d_psi_proj_span"])
            psiw0 = first_existing(cols, ["m4d_psi_w0_span"])
            em_bulk = first_existing(cols, ["m4d_em_u_bulk", "m4d_em_u_resolved"])

            final_psi0 = psi0[-1]
            final_psiproj = psiproj[-1]
            final_psiw0 = psiw0[-1]
            em_init = em_bulk[0]
            em_final = em_bulk[-1]
            em_growth = safe_ratio(em_final, em_init)
            em_max = max((abs(v) for v in em_bulk), default=math.nan)
            psi0_max = max((abs(v) for v in psi0), default=math.nan)

            finite = all(
                math.isfinite(v)
                for v in (
                    final_psi0,
                    final_psiproj,
                    final_psiw0,
                    em_init,
                    em_final,
                    em_max,
                    psi0_max,
                )
            )
            status = "PASS_RUN" if (rc == 0.0 and finite) else "FAIL_RUN"
            rows.append(
                {
                    "case": case_name,
                    "status": status,
                    "exit_code": rc,
                    "final_psi0_span": final_psi0,
                    "final_psi_proj_span": final_psiproj,
                    "final_psi_w0_span": final_psiw0,
                    "em_u_bulk_initial": em_init,
                    "em_u_bulk_final": em_final,
                    "em_u_bulk_growth": em_growth,
                    "em_u_bulk_max": em_max,
                    "psi0_max_abs": psi0_max,
                    "notes": "",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "case": case_name,
                    "status": "FAIL_PARSE",
                    "exit_code": rc,
                    "final_psi0_span": math.nan,
                    "final_psi_proj_span": math.nan,
                    "final_psi_w0_span": math.nan,
                    "em_u_bulk_initial": math.nan,
                    "em_u_bulk_final": math.nan,
                    "em_u_bulk_growth": math.nan,
                    "em_u_bulk_max": math.nan,
                    "psi0_max_abs": math.nan,
                    "notes": str(exc),
                }
            )

    rows = sorted(rows, key=lambda r: r["case"])

    baseline = next((r for r in rows if r["case"] == args.baseline_case), None)
    baseline_psi0 = baseline["final_psi0_span"] if baseline else math.nan
    baseline_emmax = baseline["em_u_bulk_max"] if baseline else math.nan

    for r in rows:
        r["psi0_vs_baseline"] = safe_ratio(r["final_psi0_span"], baseline_psi0)
        r["emmax_vs_baseline"] = safe_ratio(r["em_u_bulk_max"], baseline_emmax)

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
        "psi0_max_abs",
        "psi0_vs_baseline",
        "emmax_vs_baseline",
        "notes",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out_md = Path(args.out_md)
    lines = [
        "# Step-134 Controlled Vacuum A/B Summary",
        "",
        f"- cases: `{len(rows)}`",
        f"- baseline_case: `{args.baseline_case}`",
        "",
        "| case | status | final psi0 | final psi_proj | final psi_w0 | em_u_bulk growth | em_u_bulk max | psi0 vs baseline | emmax vs baseline |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['case']} | {r['status']} | "
            f"{r['final_psi0_span']:.6e} | {r['final_psi_proj_span']:.6e} | "
            f"{r['final_psi_w0_span']:.6e} | {r['em_u_bulk_growth']:.6e} | "
            f"{r['em_u_bulk_max']:.6e} | {r['psi0_vs_baseline']:.6e} | "
            f"{r['emmax_vs_baseline']:.6e} |"
        )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"cases,{len(rows)}")


if __name__ == "__main__":
    main()
