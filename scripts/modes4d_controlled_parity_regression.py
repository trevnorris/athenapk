#!/usr/bin/env python3
"""Quantitative controlled-limit parity checks for standard benchmarks."""

import argparse
import math
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
        raise RuntimeError(f"No history rows in {path}")

    cols = {name: [] for name in labels}
    for row in rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


def latest_hst(workdir):
    candidates = sorted(
        workdir.glob("parthenon.out*.hst"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        return None
    return candidates[-1]


def sanitize_input(
    input_path, output_path, smoke_tlim, smoke_nlim, force_modes4d_nw1=False
):
    lines = input_path.read_text(encoding="utf-8").splitlines()
    out = []
    i = 0
    has_hst = False
    has_modes4d = False
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
                has_hdf5_like = any(
                    ("file_type" in bl and ("hdf5" in bl.lower() or "phdf" in bl.lower()))
                    for bl in block_lines
                )
                has_hst_like = any(
                    ("file_type" in bl and "hst" in bl.lower()) for bl in block_lines
                )
                if has_hst_like:
                    has_hst = True
                    dt_set = False
                    forced_hst_dt = min(0.01, smoke_tlim)
                    for bl in block_lines:
                        if bl.strip().startswith("dt"):
                            out.append(f"dt = {forced_hst_dt}")
                            dt_set = True
                        else:
                            out.append(bl)
                    if not dt_set:
                        out.append(f"dt = {forced_hst_dt}")
                    continue
                if has_hdf5_like:
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

            if block_name == "parthenon/time":
                tlim_set = False
                nlim_set = False
                for bl in block_lines:
                    s = bl.strip()
                    if s.startswith("tlim"):
                        out.append(f"tlim = {smoke_tlim}")
                        tlim_set = True
                    elif s.startswith("nlim"):
                        out.append(f"nlim = {smoke_nlim}")
                        nlim_set = True
                    else:
                        out.append(bl)
                if not tlim_set:
                    out.append(f"tlim = {smoke_tlim}")
                if not nlim_set:
                    out.append(f"nlim = {smoke_nlim}")
                continue

            if block_name == "modes4d":
                has_modes4d = True
                if force_modes4d_nw1:
                    enabled_set = False
                    n_modes_set = False
                    n_quad_set = False
                    for bl in block_lines:
                        s = bl.strip()
                        if s.startswith("enabled"):
                            out.append("enabled = true")
                            enabled_set = True
                        elif s.startswith("n_modes"):
                            out.append("n_modes = 1")
                            n_modes_set = True
                        elif s.startswith("n_quadrature"):
                            out.append("n_quadrature = 3")
                            n_quad_set = True
                        else:
                            out.append(bl)
                    if not enabled_set:
                        out.append("enabled = true")
                    if not n_modes_set:
                        out.append("n_modes = 1")
                    if not n_quad_set:
                        out.append("n_quadrature = 3")
                    continue

            out.extend(block_lines)
            continue

        out.append(line)
        i += 1

    if not has_hst:
        out.extend(
            [
                "",
                "<parthenon/output2>",
                "file_type = hst",
                f"dt = {min(0.01, smoke_tlim)}",
            ]
        )

    if force_modes4d_nw1 and not has_modes4d:
        out.extend(
            [
                "",
                "<modes4d>",
                "enabled = true",
                "n_modes = 1",
                "n_quadrature = 3",
                "lambda = 1.0",
            ]
        )

    output_path.write_text("\n".join(out) + "\n", encoding="utf-8")


def run_case(binary, input_path, workdir, copy_to):
    for path in workdir.glob("parthenon.out*.hst"):
        path.unlink()
    subprocess.run([str(binary), "-i", str(input_path)], cwd=workdir, check=True)
    produced = latest_hst(workdir)
    if produced is None:
        raise RuntimeError(f"No history output produced for {input_path}")
    shutil.copyfile(produced, copy_to)
    return parse_hst(copy_to)


def rel_diff(a, b, floor):
    return abs(a - b) / max(abs(a), abs(b), floor)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_controlled_parity_outputs",
        help="Directory for generated inputs and .hst files",
    )
    parser.add_argument(
        "--smoke-tlim", type=float, default=5.0e-2, help="Override tlim for quick runs"
    )
    parser.add_argument(
        "--smoke-nlim", type=int, default=200, help="Override nlim for quick runs"
    )
    parser.add_argument(
        "--mass-rel-tol",
        type=float,
        default=1.0e-10,
        help="Relative tolerance for final mass parity",
    )
    parser.add_argument(
        "--mom-rel-tol",
        type=float,
        default=1.0e-10,
        help="Relative tolerance for final momentum parity (each component)",
    )
    parser.add_argument(
        "--energy-rel-tol",
        type=float,
        default=1.0e-10,
        help="Relative tolerance for final total-energy parity",
    )
    parser.add_argument(
        "--time-rel-tol",
        type=float,
        default=1.0e-12,
        help="Relative tolerance for final time parity",
    )
    parser.add_argument(
        "--rel-floor",
        type=float,
        default=1.0e-30,
        help="Floor used in relative-difference denominator",
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

    cases = [
        ("linear_wave3d", repo_root / "inputs" / "linear_wave3d.in"),
        ("sod", repo_root / "inputs" / "sod.in"),
        ("orszag_tang", repo_root / "inputs" / "orszag_tang.in"),
    ]

    keys = ["time", "mass", "1-mom", "2-mom", "3-mom", "tot-E"]
    rows = []
    failing_cases = []
    for name, input_path in cases:
        base_input = output_dir / f"{name}.base.in"
        mode_input = output_dir / f"{name}.modes_nw1.in"
        sanitize_input(
            input_path, base_input, args.smoke_tlim, args.smoke_nlim, force_modes4d_nw1=False
        )
        sanitize_input(
            input_path, mode_input, args.smoke_tlim, args.smoke_nlim, force_modes4d_nw1=True
        )

        base_cols = run_case(binary, base_input, workdir, output_dir / f"{name}.base.hst")
        mode_cols = run_case(binary, mode_input, workdir, output_dir / f"{name}.modes_nw1.hst")

        for k in keys:
            if k not in base_cols or k not in mode_cols:
                raise RuntimeError(f"Missing key '{k}' in case '{name}' history output")

        d_time = rel_diff(base_cols["time"][-1], mode_cols["time"][-1], args.rel_floor)
        d_mass = rel_diff(base_cols["mass"][-1], mode_cols["mass"][-1], args.rel_floor)
        d_mom1 = rel_diff(base_cols["1-mom"][-1], mode_cols["1-mom"][-1], args.rel_floor)
        d_mom2 = rel_diff(base_cols["2-mom"][-1], mode_cols["2-mom"][-1], args.rel_floor)
        d_mom3 = rel_diff(base_cols["3-mom"][-1], mode_cols["3-mom"][-1], args.rel_floor)
        d_energy = rel_diff(base_cols["tot-E"][-1], mode_cols["tot-E"][-1], args.rel_floor)

        failures = []
        if d_time > args.time_rel_tol:
            failures.append("time")
        if d_mass > args.mass_rel_tol:
            failures.append("mass")
        if d_mom1 > args.mom_rel_tol:
            failures.append("1-mom")
        if d_mom2 > args.mom_rel_tol:
            failures.append("2-mom")
        if d_mom3 > args.mom_rel_tol:
            failures.append("3-mom")
        if d_energy > args.energy_rel_tol:
            failures.append("tot-E")

        row = {
            "case": name,
            "d_time": d_time,
            "d_mass": d_mass,
            "d_1mom": d_mom1,
            "d_2mom": d_mom2,
            "d_3mom": d_mom3,
            "d_totE": d_energy,
            "status": "FAIL" if failures else "PASS",
            "failures": "|".join(failures) if failures else "none",
        }
        rows.append(row)
        if failures:
            failing_cases.append(name)

    print("case,d_time,d_mass,d_1mom,d_2mom,d_3mom,d_totE,status,failures")
    for row in rows:
        print(
            f"{row['case']},{row['d_time']:.6e},{row['d_mass']:.6e},{row['d_1mom']:.6e},"
            f"{row['d_2mom']:.6e},{row['d_3mom']:.6e},{row['d_totE']:.6e},"
            f"{row['status']},{row['failures']}"
        )

    if failing_cases:
        raise SystemExit(
            "controlled parity regression failed for: " + ", ".join(failing_cases)
        )


if __name__ == "__main__":
    main()
