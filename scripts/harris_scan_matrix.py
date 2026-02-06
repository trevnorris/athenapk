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


def analyze_case(case_name, cols):
    psi = cols["m4d_psi0_span"]
    leak = cols["m4d_int_s_leak"]
    leak_abs = cols.get("m4d_int_s_leak_abs", leak)
    jw_ew = cols["m4d_int_jw_ew"]
    epar2 = cols["m4d_brane_epar2"]
    mixed_ew2 = cols["m4d_mixed_ew2"]
    mixed_c2 = cols["m4d_mixed_c2"]
    em_a2_mode1 = maybe_col(cols, "m4d_em_a2_mode_1")
    jw_mode1_l2 = maybe_col(cols, "m4d_jw_mode_l2_1")

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
        "corr_psi0_jw_mode_l2_1": pearson(psi, jw_mode1_l2)
        if jw_mode1_l2 is not None
        else math.nan,
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
        "corr_psi0_s_leak",
        "corr_psi0_s_leak_abs",
        "corr_psi0_jw_ew",
        "corr_psi0_brane_epar2",
        "corr_psi0_mixed_ew2",
        "corr_psi0_mixed_c2",
        "corr_psi0_jw_mode_l2_1",
    ]
    print(",".join(keys))
    for row in results:
        print(",".join(fmt(row[k]) for k in keys))


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
        analyze_case("controlled", controlled_cols),
        analyze_case("full", full_cols),
    ]
    print_table(results)


if __name__ == "__main__":
    main()
