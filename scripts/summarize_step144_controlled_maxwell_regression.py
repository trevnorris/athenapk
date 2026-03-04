#!/usr/bin/env python3
"""Summarize Step-144 controlled-limit Maxwell regression runs."""

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
    for name in ("parthenon.out1.hst", "parthenon.out0.hst"):
        p = case_dir / name
        if p.exists():
            return p
    return None


def safe_ratio(a: float, b: float) -> float:
    if not math.isfinite(a) or not math.isfinite(b) or b == 0.0:
        return math.nan
    return a / b


def series_max_abs(cols: dict[str, list[float]], key: str) -> float:
    values = cols.get(key, [])
    if not values:
        return math.nan
    return max(abs(v) for v in values)


def final_value(cols: dict[str, list[float]], key: str) -> float:
    values = cols.get(key)
    if not values:
        return math.nan
    return values[-1]


def em_bulk_ledger_max_abs_rate(cols: dict[str, list[float]]) -> float:
    times = cols.get("time")
    em_u = cols.get("m4d_em_u_bulk", cols.get("m4d_em_u_resolved"))
    ja = cols.get("m4d_int_ja_ea")
    jw = cols.get("m4d_int_jw_ew")
    if (
        times is None
        or em_u is None
        or ja is None
        or jw is None
        or len(times) < 2
        or len(em_u) != len(times)
        or len(ja) != len(times)
        or len(jw) != len(times)
    ):
        return math.nan

    max_abs = 0.0
    found = False
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        du_dt = (em_u[i] - em_u[i - 1]) / dt
        dja_dt = (ja[i] - ja[i - 1]) / dt
        djw_dt = (jw[i] - jw[i - 1]) / dt
        rate = du_dt + dja_dt + djw_dt
        max_abs = max(max_abs, abs(rate))
        found = True
    return max_abs if found else math.nan


def load_check_statuses(path: Path | None):
    if path is None or not path.exists():
        return []
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--glob", required=True, help="Case directory glob")
    p.add_argument("--baseline-case", required=True)
    p.add_argument("--checks-csv", default=None)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--out-md", required=True)
    p.add_argument("--max-em-growth", type=float, default=2.0)
    p.add_argument("--max-ledger-rate", type=float, default=1.0e-1)
    p.add_argument("--max-cont-mode0", type=float, default=1.0e-4)
    p.add_argument("--quad-invariance-rel-tol", type=float, default=5.0e-2)
    p.add_argument("--quad-invariance-abs-tol", type=float, default=1.0e-12)
    args = p.parse_args()

    case_dirs = sorted(Path(x) for x in glob.glob(args.glob))
    rows = []
    for case_dir in case_dirs:
        if not case_dir.is_dir():
            continue
        case = case_dir.name
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
                    "case": case,
                    "status": "FAIL_NO_HST",
                    "exit_code": rc,
                    "final_psi0_span": math.nan,
                    "final_psi_proj_span": math.nan,
                    "final_psi_w0_span": math.nan,
                    "em_u_bulk_growth": math.nan,
                    "em_u_bulk_max": math.nan,
                    "cont_mode0_max_abs": math.nan,
                    "div_mode0_total_abs_max": math.nan,
                    "src_mode0_total_abs_max": math.nan,
                    "gauge_mode0_l2_max": math.nan,
                    "ledger_max_abs_rate": math.nan,
                    "notes": "missing history output",
                }
            )
            continue

        try:
            cols = parse_hst(hst)
            psi0 = final_value(cols, "m4d_psi0_span")
            psip = final_value(cols, "m4d_psi_proj_span")
            psiw = final_value(cols, "m4d_psi_w0_span")

            em_u_series = cols.get("m4d_em_u_bulk", cols.get("m4d_em_u_resolved"))
            if not em_u_series:
                raise RuntimeError("Missing m4d_em_u_bulk/m4d_em_u_resolved")
            em_growth = safe_ratio(em_u_series[-1], em_u_series[0])
            em_max = max(abs(v) for v in em_u_series)

            cont_mode0 = series_max_abs(cols, "m4d_cont_mode0_max_abs")
            div_mode0 = series_max_abs(cols, "m4d_int_div_mode0_total_abs")
            src_mode0 = series_max_abs(cols, "m4d_int_src_mode0_total_abs")
            gauge_mode0 = series_max_abs(cols, "m4d_gauge_mode0_l2")
            ledger = em_bulk_ledger_max_abs_rate(cols)

            finite = all(
                math.isfinite(v)
                for v in (
                    psi0,
                    psip,
                    psiw,
                    em_growth,
                    em_max,
                    cont_mode0,
                    div_mode0,
                    src_mode0,
                    gauge_mode0,
                    ledger,
                )
            )
            checks_ok = (
                em_growth <= args.max_em_growth
                and ledger <= args.max_ledger_rate
                and cont_mode0 <= args.max_cont_mode0
            )
            status = "PASS_RUN" if (rc == 0.0 and finite and checks_ok) else "FAIL_RUN"
            rows.append(
                {
                    "case": case,
                    "status": status,
                    "exit_code": rc,
                    "final_psi0_span": psi0,
                    "final_psi_proj_span": psip,
                    "final_psi_w0_span": psiw,
                    "em_u_bulk_growth": em_growth,
                    "em_u_bulk_max": em_max,
                    "cont_mode0_max_abs": cont_mode0,
                    "div_mode0_total_abs_max": div_mode0,
                    "src_mode0_total_abs_max": src_mode0,
                    "gauge_mode0_l2_max": gauge_mode0,
                    "ledger_max_abs_rate": ledger,
                    "notes": "",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "case": case,
                    "status": "FAIL_PARSE",
                    "exit_code": rc,
                    "final_psi0_span": math.nan,
                    "final_psi_proj_span": math.nan,
                    "final_psi_w0_span": math.nan,
                    "em_u_bulk_growth": math.nan,
                    "em_u_bulk_max": math.nan,
                    "cont_mode0_max_abs": math.nan,
                    "div_mode0_total_abs_max": math.nan,
                    "src_mode0_total_abs_max": math.nan,
                    "gauge_mode0_l2_max": math.nan,
                    "ledger_max_abs_rate": math.nan,
                    "notes": str(exc),
                }
            )

    rows = sorted(rows, key=lambda r: r["case"])
    by_case = {r["case"]: r for r in rows}
    baseline = by_case.get(args.baseline_case)

    invariance_rows = []
    invariance_ok = True
    compare_metrics = [
        "final_psi0_span",
        "final_psi_proj_span",
        "final_psi_w0_span",
        "em_u_bulk_max",
        "cont_mode0_max_abs",
        "div_mode0_total_abs_max",
        "src_mode0_total_abs_max",
        "ledger_max_abs_rate",
    ]
    if baseline:
        for row in rows:
            if row["case"] == args.baseline_case:
                continue
            for metric in compare_metrics:
                base_val = baseline[metric]
                val = row[metric]
                rel = math.nan
                abs_diff = math.nan
                metric_ok = False
                if math.isfinite(base_val) and math.isfinite(val):
                    abs_diff = abs(val - base_val)
                    if abs(base_val) > args.quad_invariance_abs_tol:
                        rel = abs_diff / abs(base_val)
                        metric_ok = rel <= args.quad_invariance_rel_tol
                    else:
                        metric_ok = abs_diff <= args.quad_invariance_abs_tol
                invariance_rows.append(
                    {
                        "case": row["case"],
                        "metric": metric,
                        "baseline": base_val,
                        "value": val,
                        "abs_diff": abs_diff,
                        "rel_diff": rel,
                        "status": "PASS" if metric_ok else "FAIL",
                    }
                )
                if not metric_ok:
                    invariance_ok = False

    checks_rows = load_check_statuses(Path(args.checks_csv) if args.checks_csv else None)
    checks_ok = all(r.get("status") == "PASS" for r in checks_rows)
    runs_ok = all(r["status"] == "PASS_RUN" for r in rows)
    overall = "PASS" if (checks_ok and runs_ok and invariance_ok) else "FAIL"

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case",
        "status",
        "exit_code",
        "final_psi0_span",
        "final_psi_proj_span",
        "final_psi_w0_span",
        "em_u_bulk_growth",
        "em_u_bulk_max",
        "cont_mode0_max_abs",
        "div_mode0_total_abs_max",
        "src_mode0_total_abs_max",
        "gauge_mode0_l2_max",
        "ledger_max_abs_rate",
        "notes",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    inv_csv = out_csv.with_name(out_csv.stem + "_invariance.csv")
    with inv_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "metric",
                "baseline",
                "value",
                "abs_diff",
                "rel_diff",
                "status",
            ],
        )
        w.writeheader()
        for r in invariance_rows:
            w.writerow(r)

    out_md = Path(args.out_md)
    lines = [
        "# Step-144 Controlled Maxwell Regression Summary",
        "",
        f"- cases: `{len(rows)}`",
        f"- baseline_case: `{args.baseline_case}`",
        f"- thresholds: em_growth<={args.max_em_growth:.3e}, ledger_rate<={args.max_ledger_rate:.3e}, cont_mode0<={args.max_cont_mode0:.3e}",
        f"- invariance tolerance: rel<={args.quad_invariance_rel_tol:.3e}, abs<={args.quad_invariance_abs_tol:.3e}",
        f"- overall_status: `{overall}`",
        "",
        "## Runtime Cases",
        "",
        "| case | status | final psi0 | final psi_proj | final psi_w0 | em growth | em max | cont_mode0 max | div_mode0 max | src_mode0 max | gauge_mode0_l2 max | ledger max rate |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['case']} | {r['status']} | "
            f"{r['final_psi0_span']:.6e} | {r['final_psi_proj_span']:.6e} | {r['final_psi_w0_span']:.6e} | "
            f"{r['em_u_bulk_growth']:.6e} | {r['em_u_bulk_max']:.6e} | {r['cont_mode0_max_abs']:.6e} | "
            f"{r['div_mode0_total_abs_max']:.6e} | {r['src_mode0_total_abs_max']:.6e} | "
            f"{r['gauge_mode0_l2_max']:.6e} | {r['ledger_max_abs_rate']:.6e} |"
        )

    lines.extend(
        [
            "",
            "## Transform/Quadrature Invariance vs Baseline",
            "",
            f"- status: `{'PASS' if invariance_ok else 'FAIL'}`",
            "",
            "| case | metric | baseline | value | abs diff | rel diff | status |",
            "|---|---|---:|---:|---:|---:|---|",
        ]
    )
    if invariance_rows:
        for r in invariance_rows:
            lines.append(
                f"| {r['case']} | {r['metric']} | {r['baseline']:.6e} | {r['value']:.6e} | "
                f"{r['abs_diff']:.6e} | {r['rel_diff']:.6e} | {r['status']} |"
            )
    else:
        lines.append("| n/a | n/a | n/a | n/a | n/a | n/a | n/a |")

    lines.extend(
        [
            "",
            "## Auxiliary Checks",
            "",
            f"- status: `{'PASS' if checks_ok else 'FAIL'}`",
        ]
    )
    if checks_rows:
        lines.extend(
            [
                "",
                "| check | status | log |",
                "|---|---|---|",
            ]
        )
        for r in checks_rows:
            lines.append(f"| `{r.get('check', '')}` | `{r.get('status', '')}` | `{r.get('log', '')}` |")
    else:
        lines.append("")
        lines.append("- no auxiliary checks provided")

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"summary_csv,{out_csv}")
    print(f"invariance_csv,{inv_csv}")
    print(f"summary_md,{out_md}")
    print(f"overall_status,{overall}")
    print(f"cases,{len(rows)}")


if __name__ == "__main__":
    main()
