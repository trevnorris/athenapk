#!/usr/bin/env python3
"""Run controlled-limit regression smoke checks for AthenaPK + modes4d."""

import argparse
import math
import shutil
import subprocess
import sys
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


def latest_hst(workdir):
    candidates = [workdir / "parthenon.out1.hst", workdir / "parthenon.out0.hst"]
    for path in candidates:
        if path.exists():
            return path
    return None


def sanitize_input_for_disabled_hdf5(
    input_path, sanitized_path, smoke_tlim_override=None, smoke_nlim_override=None
):
    lines = input_path.read_text(encoding="utf-8").splitlines()
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith("<") and stripped.endswith(">"):
            block_name = stripped[1:-1].strip().lower()
            out.append(line)
            i += 1
            block_lines = []
            while i < len(lines):
                nxt = lines[i]
                nxt_stripped = nxt.strip()
                if nxt_stripped.startswith("<") and nxt_stripped.endswith(">"):
                    break
                block_lines.append(nxt)
                i += 1

            if block_name.startswith("parthenon/output"):
                has_hdf5 = any(
                    ("file_type" in bl and ("hdf5" in bl.lower() or "phdf" in bl.lower()))
                    for bl in block_lines
                )
                if has_hdf5:
                    dt_set = False
                    for bl in block_lines:
                        if bl.strip().startswith("dt"):
                            out.append("dt = -1")
                            dt_set = True
                        else:
                            out.append(bl)
                    if not dt_set:
                        out.append("dt = -1")
                    continue

            if block_name == "parthenon/time" and (
                smoke_tlim_override is not None or smoke_nlim_override is not None
            ):
                tlim_set = False
                nlim_set = False
                for bl in block_lines:
                    stripped_bl = bl.strip()
                    if stripped_bl.startswith("tlim") and smoke_tlim_override is not None:
                        out.append(f"tlim = {smoke_tlim_override}")
                        tlim_set = True
                    elif stripped_bl.startswith("nlim") and smoke_nlim_override is not None:
                        out.append(f"nlim = {smoke_nlim_override}")
                        nlim_set = True
                    else:
                        out.append(bl)
                if smoke_tlim_override is not None and not tlim_set:
                    out.append(f"tlim = {smoke_tlim_override}")
                if smoke_nlim_override is not None and not nlim_set:
                    out.append(f"nlim = {smoke_nlim_override}")
                continue

            out.extend(block_lines)
            continue

        out.append(line)
        i += 1

    sanitized_path.write_text("\n".join(out) + "\n", encoding="utf-8")


def run_input(binary, input_path, workdir, copy_to):
    for filename in ("parthenon.out0.hst", "parthenon.out1.hst"):
        try:
            (workdir / filename).unlink()
        except FileNotFoundError:
            pass

    cmd = [str(binary), "-i", str(input_path)]
    subprocess.run(cmd, cwd=workdir, check=True)

    produced = latest_hst(workdir)
    if produced is None:
        return None
    copy_to.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(produced, copy_to)
    return copy_to


def max_abs(values):
    return max(abs(v) for v in values) if values else math.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_controlled_regression_outputs",
        help="Directory to store copied .hst files",
    )
    parser.add_argument(
        "--max-controlled-jw-ew-abs",
        type=float,
        default=1.0e-12,
        help="Maximum allowed |m4d_int_jw_ew| in controlled run",
    )
    parser.add_argument(
        "--max-controlled-s-leak-abs",
        type=float,
        default=1.0e-12,
        help="Maximum allowed m4d_int_s_leak_abs in controlled run",
    )
    parser.add_argument(
        "--max-controlled-mixed-ew2",
        type=float,
        default=1.0e-12,
        help="Maximum allowed m4d_mixed_ew2 in controlled run",
    )
    parser.add_argument(
        "--max-controlled-mixed-c2",
        type=float,
        default=1.0e-12,
        help="Maximum allowed m4d_mixed_c2 in controlled run",
    )
    parser.add_argument(
        "--smoke-tlim",
        type=float,
        default=5.0e-2,
        help="Override tlim for baseline smoke cases",
    )
    parser.add_argument(
        "--smoke-nlim",
        type=int,
        default=200,
        help="Override nlim for baseline smoke cases",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]

    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    smoke_inputs = [
        ("linear_wave3d", repo_root / "inputs" / "linear_wave3d.in"),
        ("sod", repo_root / "inputs" / "sod.in"),
        ("orszag_tang", repo_root / "inputs" / "orszag_tang.in"),
    ]

    for name, input_path in smoke_inputs:
        sanitized_input = output_dir / f"{name}.sanitized.in"
        sanitize_input_for_disabled_hdf5(
            input_path, sanitized_input, args.smoke_tlim, args.smoke_nlim
        )
        run_input(binary, sanitized_input, workdir, output_dir / f"{name}.hst")

    controlled_input = repo_root / "inputs" / "harris_4d_controlled.in"
    controlled_sanitized = output_dir / "harris_4d_controlled.sanitized.in"
    sanitize_input_for_disabled_hdf5(controlled_input, controlled_sanitized)
    controlled_hst = run_input(
        binary,
        controlled_sanitized,
        workdir,
        output_dir / "harris_controlled.hst",
    )
    if controlled_hst is None:
        raise RuntimeError("Controlled run did not produce a history output file")

    cols = parse_hst(controlled_hst)
    required = ("m4d_int_jw_ew", "m4d_int_s_leak_abs", "m4d_mixed_ew2", "m4d_mixed_c2")
    missing = [key for key in required if key not in cols]
    if missing:
        raise RuntimeError(f"Missing required controlled-run columns: {', '.join(missing)}")

    controlled_max_jw_ew_abs = max_abs(cols["m4d_int_jw_ew"])
    controlled_max_s_leak_abs = max(cols["m4d_int_s_leak_abs"])
    controlled_max_mixed_ew2 = max(cols["m4d_mixed_ew2"])
    controlled_max_mixed_c2 = max(cols["m4d_mixed_c2"])

    failures = []
    if controlled_max_jw_ew_abs > args.max_controlled_jw_ew_abs:
        failures.append("jw_ew")
    if controlled_max_s_leak_abs > args.max_controlled_s_leak_abs:
        failures.append("s_leak_abs")
    if controlled_max_mixed_ew2 > args.max_controlled_mixed_ew2:
        failures.append("mixed_ew2")
    if controlled_max_mixed_c2 > args.max_controlled_mixed_c2:
        failures.append("mixed_c2")

    status = "FAIL" if failures else "PASS"
    fail_text = "|".join(failures) if failures else "none"
    print(
        "case,max_abs_int_jw_ew,max_int_s_leak_abs,max_mixed_ew2,max_mixed_c2,status,failures"
    )
    print(
        f"controlled,{controlled_max_jw_ew_abs:.6e},{controlled_max_s_leak_abs:.6e},"
        f"{controlled_max_mixed_ew2:.6e},{controlled_max_mixed_c2:.6e},{status},{fail_text}"
    )

    if failures:
        raise SystemExit(f"controlled regression checks failed: {', '.join(failures)}")


if __name__ == "__main__":
    main()
