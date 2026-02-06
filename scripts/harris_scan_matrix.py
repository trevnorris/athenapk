#!/usr/bin/env python3
# AthenaPK 4D modes Harris scan helper:
# runs controlled/full input decks and reports channel correlations.

import argparse
import math
import os
import shutil
import subprocess
from pathlib import Path


def parse_hst(path):
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
        raise RuntimeError(f"No history data rows in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


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


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def maybe_col(cols, key):
    return cols.get(key)


def continuity_closure_metrics(times, charge_mode0, int_s_leak, int_divj_mode0, abs_rate_tol):
    if len(times) < 2 or len(charge_mode0) != len(times) or len(int_s_leak) != len(times):
        return (math.nan, math.nan, math.nan, math.nan)
    if int_divj_mode0 is not None and len(int_divj_mode0) != len(times):
        return (math.nan, math.nan, math.nan, math.nan)

    norm_residuals = []
    abs_residuals = []
    final_residual = math.nan
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        dcharge_dt = (charge_mode0[i] - charge_mode0[i - 1]) / dt
        dsleak_dt = (int_s_leak[i] - int_s_leak[i - 1]) / dt
        ddivj_dt = 0.0
        if int_divj_mode0 is not None:
            ddivj_dt = (int_divj_mode0[i] - int_divj_mode0[i - 1]) / dt
        residual = dcharge_dt - dsleak_dt + ddivj_dt
        scale = max(abs(dcharge_dt) + abs(dsleak_dt) + abs(ddivj_dt), abs_rate_tol)
        norm_residuals.append(abs(residual) / scale)
        abs_residuals.append(abs(residual))
        final_residual = residual

    if not norm_residuals:
        return (math.nan, math.nan, math.nan, final_residual)

    max_norm = max(norm_residuals)
    rms_norm = math.sqrt(sum(x * x for x in norm_residuals) / len(norm_residuals))
    max_abs = max(abs_residuals)
    return (max_norm, rms_norm, max_abs, final_residual)


def analyze_case(
    case_name, cols, closure_norm_tol, closure_abs_rate_tol, local_mode0_abs_rate_tol
):
    times = cols["time"]
    psi = cols["m4d_psi0_span"]
    leak = cols["m4d_int_s_leak"]
    leak_abs = cols.get("m4d_int_s_leak_abs", leak)
    jw_ew = cols["m4d_int_jw_ew"]
    epar2 = cols["m4d_brane_epar2"]
    mixed_ew2 = cols["m4d_mixed_ew2"]
    mixed_c2 = cols["m4d_mixed_c2"]
    em_a2_mode1 = maybe_col(cols, "m4d_em_a2_mode_1")
    jw_mode1_l2 = maybe_col(cols, "m4d_jw_mode_l2_1")
    jw_mode1 = maybe_col(cols, "m4d_jw_mode_1")
    charge_mode0 = maybe_col(cols, "m4d_charge_mode_0")
    momx_mode0 = maybe_col(cols, "m4d_momx_mode_0")
    momy_mode0 = maybe_col(cols, "m4d_momy_mode_0")
    momz_mode0 = maybe_col(cols, "m4d_momz_mode_0")
    momw_mode0 = maybe_col(cols, "m4d_momw_mode_0")
    energy_mode0 = maybe_col(cols, "m4d_energy_mode_0")
    int_divj_mode0 = maybe_col(cols, "m4d_int_divj_mode0")
    int_divmomx_mode0 = maybe_col(cols, "m4d_int_divmomx_mode0")
    int_divmomy_mode0 = maybe_col(cols, "m4d_int_divmomy_mode0")
    int_divmomz_mode0 = maybe_col(cols, "m4d_int_divmomz_mode0")
    int_divmomw_mode0 = maybe_col(cols, "m4d_int_divmomw_mode0")
    int_divenergy_mode0 = maybe_col(cols, "m4d_int_divenergy_mode0")
    cont_local_mode0_max_abs = maybe_col(cols, "m4d_cont_mode0_max_abs")
    cont_local_mode0_l1 = maybe_col(cols, "m4d_cont_mode0_l1")
    cont_local_mode0_l2 = maybe_col(cols, "m4d_cont_mode0_l2")

    closure_max_norm = math.nan
    closure_rms_norm = math.nan
    closure_max_abs_rate = math.nan
    closure_final_residual = math.nan
    closure_status = "N/A"
    closure_local_mode0_status = "N/A"
    if charge_mode0 is not None:
        (
            closure_max_norm,
            closure_rms_norm,
            closure_max_abs_rate,
            closure_final_residual,
        ) = continuity_closure_metrics(
            times, charge_mode0, leak, int_divj_mode0, closure_abs_rate_tol
        )
        if not math.isnan(closure_max_norm):
            closure_status = (
                "PASS"
                if (
                    closure_max_norm <= closure_norm_tol
                    or closure_max_abs_rate <= closure_abs_rate_tol
                )
                else "FAIL"
            )
    if cont_local_mode0_max_abs is not None:
        closure_local_mode0_status = (
            "PASS"
            if max(cont_local_mode0_max_abs) <= local_mode0_abs_rate_tol
            else "FAIL"
        )

    result = {
        "case": case_name,
        "final_psi0_span": psi[-1],
        "final_s_leak": leak[-1],
        "final_s_leak_abs": leak_abs[-1],
        "final_jw_ew": jw_ew[-1],
        "final_brane_epar2": epar2[-1],
        "final_mixed_ew2": mixed_ew2[-1],
        "final_mixed_c2": mixed_c2[-1],
        "corr_psi0_s_leak": pearson(psi, leak),
        "corr_psi0_s_leak_abs": pearson(psi, leak_abs),
        "corr_psi0_jw_ew": pearson(psi, jw_ew),
        "corr_psi0_brane_epar2": pearson(psi, epar2),
        "corr_psi0_mixed_ew2": pearson(psi, mixed_ew2),
        "corr_psi0_mixed_c2": pearson(psi, mixed_c2),
        "final_em_a2_mode_1": em_a2_mode1[-1] if em_a2_mode1 is not None else math.nan,
        "final_jw_mode_l2_1": jw_mode1_l2[-1] if jw_mode1_l2 is not None else math.nan,
        "final_jw_mode_1": jw_mode1[-1] if jw_mode1 is not None else math.nan,
        "final_charge_mode_0": charge_mode0[-1] if charge_mode0 is not None else math.nan,
        "final_momx_mode_0": momx_mode0[-1] if momx_mode0 is not None else math.nan,
        "final_momy_mode_0": momy_mode0[-1] if momy_mode0 is not None else math.nan,
        "final_momz_mode_0": momz_mode0[-1] if momz_mode0 is not None else math.nan,
        "final_momw_mode_0": momw_mode0[-1] if momw_mode0 is not None else math.nan,
        "final_energy_mode_0": energy_mode0[-1] if energy_mode0 is not None else math.nan,
        "final_int_divj_mode0": int_divj_mode0[-1] if int_divj_mode0 is not None else math.nan,
        "final_int_divmomx_mode0": int_divmomx_mode0[-1]
        if int_divmomx_mode0 is not None
        else math.nan,
        "final_int_divmomy_mode0": int_divmomy_mode0[-1]
        if int_divmomy_mode0 is not None
        else math.nan,
        "final_int_divmomz_mode0": int_divmomz_mode0[-1]
        if int_divmomz_mode0 is not None
        else math.nan,
        "final_int_divmomw_mode0": int_divmomw_mode0[-1]
        if int_divmomw_mode0 is not None
        else math.nan,
        "final_int_divenergy_mode0": int_divenergy_mode0[-1]
        if int_divenergy_mode0 is not None
        else math.nan,
        "corr_psi0_jw_mode_l2_1": pearson(psi, jw_mode1_l2)
        if jw_mode1_l2 is not None
        else math.nan,
        "closure_max_norm": closure_max_norm,
        "closure_rms_norm": closure_rms_norm,
        "closure_max_abs_rate": closure_max_abs_rate,
        "closure_final_residual": closure_final_residual,
        "closure_status": closure_status,
        "closure_local_mode0_max_abs_rate": max(cont_local_mode0_max_abs)
        if cont_local_mode0_max_abs is not None
        else math.nan,
        "closure_local_mode0_final_max_abs_rate": cont_local_mode0_max_abs[-1]
        if cont_local_mode0_max_abs is not None
        else math.nan,
        "closure_local_mode0_final_l1": cont_local_mode0_l1[-1]
        if cont_local_mode0_l1 is not None
        else math.nan,
        "closure_local_mode0_final_l2": cont_local_mode0_l2[-1]
        if cont_local_mode0_l2 is not None
        else math.nan,
        "closure_local_mode0_status": closure_local_mode0_status,
    }
    return result


def run_case(binary, input_path, workdir, output_hst):
    for filename in ("parthenon.out0.hst", "parthenon.out1.hst"):
        try:
            (workdir / filename).unlink()
        except FileNotFoundError:
            pass

    cmd = [str(binary), "-i", str(input_path)]
    subprocess.run(cmd, cwd=workdir, check=True)

    produced = workdir / "parthenon.out1.hst"
    if not produced.exists():
        raise RuntimeError(f"Expected history output not found: {produced}")
    shutil.copyfile(produced, output_hst)
    return parse_hst(output_hst)


def print_table(results):
    keys = [
        "case",
        "final_psi0_span",
        "final_s_leak",
        "final_s_leak_abs",
        "final_jw_ew",
        "final_brane_epar2",
        "final_mixed_ew2",
        "final_mixed_c2",
        "final_em_a2_mode_1",
        "final_jw_mode_l2_1",
        "final_jw_mode_1",
        "final_charge_mode_0",
        "final_momx_mode_0",
        "final_momy_mode_0",
        "final_momz_mode_0",
        "final_momw_mode_0",
        "final_energy_mode_0",
        "final_int_divj_mode0",
        "final_int_divmomx_mode0",
        "final_int_divmomy_mode0",
        "final_int_divmomz_mode0",
        "final_int_divmomw_mode0",
        "final_int_divenergy_mode0",
        "corr_psi0_s_leak",
        "corr_psi0_s_leak_abs",
        "corr_psi0_jw_ew",
        "corr_psi0_brane_epar2",
        "corr_psi0_mixed_ew2",
        "corr_psi0_mixed_c2",
        "corr_psi0_jw_mode_l2_1",
        "closure_max_norm",
        "closure_rms_norm",
        "closure_max_abs_rate",
        "closure_final_residual",
        "closure_status",
        "closure_local_mode0_max_abs_rate",
        "closure_local_mode0_final_max_abs_rate",
        "closure_local_mode0_final_l1",
        "closure_local_mode0_final_l2",
        "closure_local_mode0_status",
        "activity_status",
        "activity_failures",
    ]
    print(",".join(keys))
    for row in results:
        print(",".join(fmt(row[k]) for k in keys))


def evaluate_full_activity(
    row,
    min_jw_ew_abs,
    min_s_leak_abs,
    min_mixed_ew2,
    min_jw_mode_l2_1,
):
    failures = []
    if min_jw_ew_abs > 0.0 and abs(row["final_jw_ew"]) < min_jw_ew_abs:
        failures.append("jw_ew")
    if min_s_leak_abs > 0.0 and row["final_s_leak_abs"] < min_s_leak_abs:
        failures.append("s_leak_abs")
    if min_mixed_ew2 > 0.0 and row["final_mixed_ew2"] < min_mixed_ew2:
        failures.append("mixed_ew2")
    if min_jw_mode_l2_1 > 0.0 and row["final_jw_mode_l2_1"] < min_jw_mode_l2_1:
        failures.append("jw_mode_l2_1")

    if failures:
        return ("FAIL", "|".join(failures))
    return ("PASS", "none")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Directory where athenaPK runs")
    parser.add_argument(
        "--output-dir",
        default="harris_scan_outputs",
        help="Directory to store copied .hst files",
    )
    parser.add_argument(
        "--controlled-input",
        default="inputs/harris_4d_controlled.in",
        help="Controlled-limit input deck",
    )
    parser.add_argument(
        "--full-input",
        default="inputs/harris_4d_full.in",
        help="Full-channel input deck",
    )
    parser.add_argument(
        "--closure-norm-tol",
        type=float,
        default=5.0e-2,
        help="Maximum normalized continuity residual for PASS",
    )
    parser.add_argument(
        "--closure-abs-rate-tol",
        type=float,
        default=1.0e-8,
        help="Maximum absolute continuity residual rate for low-signal PASS",
    )
    parser.add_argument(
        "--closure-local-mode0-abs-rate-tol",
        type=float,
        default=1.0e-8,
        help="Maximum mode-0 local continuity residual rate for PASS",
    )
    parser.add_argument(
        "--fail-on-check",
        action="store_true",
        help="Exit nonzero if any case fails closure or activity checks",
    )
    parser.add_argument(
        "--full-min-jw-ew-abs",
        type=float,
        default=0.0,
        help="Minimum |final_jw_ew| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-s-leak-abs",
        type=float,
        default=0.0,
        help="Minimum final_s_leak_abs required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-mixed-ew2",
        type=float,
        default=0.0,
        help="Minimum final_mixed_ew2 required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-jw-mode-l2-1",
        type=float,
        default=0.0,
        help="Minimum final_jw_mode_l2_1 required for full-case activity PASS (disabled when 0)",
    )
    args = parser.parse_args()

    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    controlled_input = Path(args.controlled_input).resolve()
    full_input = Path(args.full_input).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    controlled_cols = run_case(
        binary, controlled_input, workdir, output_dir / "harris_controlled.out1.hst"
    )
    full_cols = run_case(binary, full_input, workdir, output_dir / "harris_full.out1.hst")

    results = [
        analyze_case(
            "controlled",
            controlled_cols,
            args.closure_norm_tol,
            args.closure_abs_rate_tol,
            args.closure_local_mode0_abs_rate_tol,
        ),
        analyze_case(
            "full",
            full_cols,
            args.closure_norm_tol,
            args.closure_abs_rate_tol,
            args.closure_local_mode0_abs_rate_tol,
        ),
    ]
    for row in results:
        row["activity_status"] = "N/A"
        row["activity_failures"] = "n/a"
        if row["case"] == "full":
            status, failures = evaluate_full_activity(
                row,
                args.full_min_jw_ew_abs,
                args.full_min_s_leak_abs,
                args.full_min_mixed_ew2,
                args.full_min_jw_mode_l2_1,
            )
            row["activity_status"] = status
            row["activity_failures"] = failures

    print_table(results)

    if args.fail_on_check:
        failing = [
            r["case"]
            for r in results
            if r["closure_status"] == "FAIL"
            or r["closure_local_mode0_status"] == "FAIL"
            or r["activity_status"] == "FAIL"
        ]
        if failing:
            raise SystemExit(f"scan checks failed for: {', '.join(failing)}")


if __name__ == "__main__":
    main()
