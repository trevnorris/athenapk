#!/usr/bin/env python3
"""Regression checks for the optional conservative EM transport path."""

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


def require_col(cols, key):
    values = cols.get(key)
    if values is None:
        raise RuntimeError(f"Missing required history column: {key}")
    return values


def ensure_finite(values, label):
    for value in values:
        if not math.isfinite(value):
            raise RuntimeError(f"Non-finite value in {label}")


def run_case(binary, input_path, workdir, output_hst, extra_args):
    for filename in ("parthenon.out0.hst", "parthenon.out1.hst"):
        try:
            (workdir / filename).unlink()
        except FileNotFoundError:
            pass

    cmd = [str(binary), "-i", str(input_path)]
    cmd.extend(extra_args)
    subprocess.run(cmd, cwd=workdir, check=True)

    produced = None
    for candidate in ("parthenon.out1.hst", "parthenon.out0.hst"):
        path = workdir / candidate
        if path.exists():
            produced = path
            break
    if produced is None:
        raise RuntimeError("Expected history output not found (parthenon.out1.hst/out0.hst)")

    shutil.copyfile(produced, output_hst)
    return parse_hst(output_hst)


def fmt(value):
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def analyze_scalar(cols, min_mixed_ew2, min_mixed_c2, min_em_a2_mode1, min_pulse_xc_shift):
    mixed_ew2 = require_col(cols, "m4d_mixed_ew2")
    mixed_c2 = require_col(cols, "m4d_mixed_c2")
    em_a2_mode1 = require_col(cols, "m4d_em_a2_mode_1")
    pulse_xc = require_col(cols, "m4d_pulse_xc")

    ensure_finite(mixed_ew2, "scalar:m4d_mixed_ew2")
    ensure_finite(mixed_c2, "scalar:m4d_mixed_c2")
    ensure_finite(em_a2_mode1, "scalar:m4d_em_a2_mode_1")
    ensure_finite(pulse_xc, "scalar:m4d_pulse_xc")

    max_mixed_ew2 = max(mixed_ew2)
    max_mixed_c2 = max(mixed_c2)
    max_em_a2_mode1 = max(em_a2_mode1)
    pulse_xc_shift = abs(pulse_xc[-1] - pulse_xc[0])

    failures = []
    if min_mixed_ew2 > 0.0 and max_mixed_ew2 < min_mixed_ew2:
        failures.append("mixed_ew2")
    if min_mixed_c2 > 0.0 and max_mixed_c2 < min_mixed_c2:
        failures.append("mixed_c2")
    if min_em_a2_mode1 > 0.0 and max_em_a2_mode1 < min_em_a2_mode1:
        failures.append("em_a2_mode_1")
    if min_pulse_xc_shift > 0.0 and pulse_xc_shift < min_pulse_xc_shift:
        failures.append("pulse_xc_shift")

    return {
        "case": "scalar_pulse_conservative",
        "max_mixed_ew2": max_mixed_ew2,
        "max_mixed_c2": max_mixed_c2,
        "max_em_a2_mode_1": max_em_a2_mode1,
        "pulse_xc_shift": pulse_xc_shift,
        "status": "FAIL" if failures else "PASS",
        "failures": "|".join(failures) if failures else "none",
    }


def transport_balance_metrics(
    times,
    quantity_mode0,
    int_divquantity_mode0,
    int_srcquantity_mode0,
    abs_rate_tol,
    skip_initial_interval=True,
):
    if (
        len(times) < 2
        or quantity_mode0 is None
        or int_divquantity_mode0 is None
        or int_srcquantity_mode0 is None
        or len(quantity_mode0) != len(times)
        or len(int_divquantity_mode0) != len(times)
        or len(int_srcquantity_mode0) != len(times)
    ):
        return (math.nan, math.nan, math.nan, math.nan)

    norm_rates = []
    abs_rates = []
    final_rate = math.nan
    for i in range(1, len(times)):
        if skip_initial_interval and i == 1:
            continue
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        dquantity_dt = (quantity_mode0[i] - quantity_mode0[i - 1]) / dt
        ddiv_dt = (int_divquantity_mode0[i] - int_divquantity_mode0[i - 1]) / dt
        dsrc_dt = (int_srcquantity_mode0[i] - int_srcquantity_mode0[i - 1]) / dt
        rate = dquantity_dt + ddiv_dt - dsrc_dt
        scale = max(abs(dquantity_dt) + abs(ddiv_dt) + abs(dsrc_dt), abs_rate_tol)
        norm_rates.append(abs(rate) / scale)
        abs_rates.append(abs(rate))
        final_rate = rate

    if not norm_rates:
        return (math.nan, math.nan, math.nan, final_rate)

    max_norm = max(norm_rates)
    rms_norm = math.sqrt(sum(x * x for x in norm_rates) / len(norm_rates))
    max_abs = max(abs_rates)
    return (max_norm, rms_norm, max_abs, final_rate)


def analyze_harris(
    cols, min_mixed_ew2, min_mixed_c2, transport_norm_tol, transport_abs_rate_tol
):
    times = require_col(cols, "time")
    mixed_ew2 = require_col(cols, "m4d_mixed_ew2")
    mixed_c2 = require_col(cols, "m4d_mixed_c2")
    cont_mode0_max_abs = require_col(cols, "m4d_cont_mode0_max_abs")
    int_jw_ew = require_col(cols, "m4d_int_jw_ew")
    int_s_leak_abs = require_col(cols, "m4d_int_s_leak_abs")
    pi0_mode0 = require_col(cols, "m4d_pi0_mode_0")
    pix_mode0 = require_col(cols, "m4d_pix_mode_0")
    piy_mode0 = require_col(cols, "m4d_piy_mode_0")
    piz_mode0 = require_col(cols, "m4d_piz_mode_0")
    piw_mode0 = require_col(cols, "m4d_piw_mode_0")
    int_divpi0_mode0 = require_col(cols, "m4d_int_divpi0_mode0")
    int_divpix_mode0 = require_col(cols, "m4d_int_divpix_mode0")
    int_divpiy_mode0 = require_col(cols, "m4d_int_divpiy_mode0")
    int_divpiz_mode0 = require_col(cols, "m4d_int_divpiz_mode0")
    int_divpiw_mode0 = require_col(cols, "m4d_int_divpiw_mode0")
    int_srcpi0_mode0 = require_col(cols, "m4d_int_srcpi0_mode0")
    int_srcpix_mode0 = require_col(cols, "m4d_int_srcpix_mode0")
    int_srcpiy_mode0 = require_col(cols, "m4d_int_srcpiy_mode0")
    int_srcpiz_mode0 = require_col(cols, "m4d_int_srcpiz_mode0")
    int_srcpiw_mode0 = require_col(cols, "m4d_int_srcpiw_mode0")

    ensure_finite(times, "harris:time")
    ensure_finite(mixed_ew2, "harris:m4d_mixed_ew2")
    ensure_finite(mixed_c2, "harris:m4d_mixed_c2")
    ensure_finite(cont_mode0_max_abs, "harris:m4d_cont_mode0_max_abs")
    ensure_finite(int_jw_ew, "harris:m4d_int_jw_ew")
    ensure_finite(int_s_leak_abs, "harris:m4d_int_s_leak_abs")
    ensure_finite(pi0_mode0, "harris:m4d_pi0_mode_0")
    ensure_finite(pix_mode0, "harris:m4d_pix_mode_0")
    ensure_finite(piy_mode0, "harris:m4d_piy_mode_0")
    ensure_finite(piz_mode0, "harris:m4d_piz_mode_0")
    ensure_finite(piw_mode0, "harris:m4d_piw_mode_0")
    ensure_finite(int_divpi0_mode0, "harris:m4d_int_divpi0_mode0")
    ensure_finite(int_divpix_mode0, "harris:m4d_int_divpix_mode0")
    ensure_finite(int_divpiy_mode0, "harris:m4d_int_divpiy_mode0")
    ensure_finite(int_divpiz_mode0, "harris:m4d_int_divpiz_mode0")
    ensure_finite(int_divpiw_mode0, "harris:m4d_int_divpiw_mode0")
    ensure_finite(int_srcpi0_mode0, "harris:m4d_int_srcpi0_mode0")
    ensure_finite(int_srcpix_mode0, "harris:m4d_int_srcpix_mode0")
    ensure_finite(int_srcpiy_mode0, "harris:m4d_int_srcpiy_mode0")
    ensure_finite(int_srcpiz_mode0, "harris:m4d_int_srcpiz_mode0")
    ensure_finite(int_srcpiw_mode0, "harris:m4d_int_srcpiw_mode0")

    max_mixed_ew2 = max(mixed_ew2)
    max_mixed_c2 = max(mixed_c2)
    max_cont_mode0 = max(cont_mode0_max_abs)
    final_jw_ew = int_jw_ew[-1]
    final_s_leak_abs = int_s_leak_abs[-1]

    (
        pi0_transport_max_norm,
        _,
        pi0_transport_max_abs_rate,
        _,
    ) = transport_balance_metrics(
        times, pi0_mode0, int_divpi0_mode0, int_srcpi0_mode0, transport_abs_rate_tol
    )
    (
        pix_transport_max_norm,
        _,
        pix_transport_max_abs_rate,
        _,
    ) = transport_balance_metrics(
        times, pix_mode0, int_divpix_mode0, int_srcpix_mode0, transport_abs_rate_tol
    )
    (
        piy_transport_max_norm,
        _,
        piy_transport_max_abs_rate,
        _,
    ) = transport_balance_metrics(
        times, piy_mode0, int_divpiy_mode0, int_srcpiy_mode0, transport_abs_rate_tol
    )
    (
        piz_transport_max_norm,
        _,
        piz_transport_max_abs_rate,
        _,
    ) = transport_balance_metrics(
        times, piz_mode0, int_divpiz_mode0, int_srcpiz_mode0, transport_abs_rate_tol
    )
    (
        piw_transport_max_norm,
        _,
        piw_transport_max_abs_rate,
        _,
    ) = transport_balance_metrics(
        times, piw_mode0, int_divpiw_mode0, int_srcpiw_mode0, transport_abs_rate_tol
    )

    def rate_status(max_norm, max_abs):
        if math.isnan(max_norm):
            return "N/A"
        return (
            "PASS"
            if (max_norm <= transport_norm_tol or max_abs <= transport_abs_rate_tol)
            else "FAIL"
        )

    pi0_transport_status = rate_status(
        pi0_transport_max_norm, pi0_transport_max_abs_rate
    )
    pix_transport_status = rate_status(
        pix_transport_max_norm, pix_transport_max_abs_rate
    )
    piy_transport_status = rate_status(
        piy_transport_max_norm, piy_transport_max_abs_rate
    )
    piz_transport_status = rate_status(
        piz_transport_max_norm, piz_transport_max_abs_rate
    )
    piw_transport_status = rate_status(
        piw_transport_max_norm, piw_transport_max_abs_rate
    )
    em_transport_status = (
        "PASS"
        if all(
            s == "PASS"
            for s in (
                pi0_transport_status,
                pix_transport_status,
                piy_transport_status,
                piz_transport_status,
                piw_transport_status,
            )
        )
        else "FAIL"
    )

    failures = []
    if max_mixed_ew2 < min_mixed_ew2:
        failures.append("mixed_ew2")
    if max_mixed_c2 < min_mixed_c2:
        failures.append("mixed_c2")
    if pi0_transport_status != "PASS":
        failures.append("pi0_transport")
    if pix_transport_status != "PASS":
        failures.append("pix_transport")
    if piy_transport_status != "PASS":
        failures.append("piy_transport")
    if piz_transport_status != "PASS":
        failures.append("piz_transport")
    if piw_transport_status != "PASS":
        failures.append("piw_transport")

    return {
        "case": "harris_short_conservative",
        "max_mixed_ew2": max_mixed_ew2,
        "max_mixed_c2": max_mixed_c2,
        "max_cont_mode0_abs": max_cont_mode0,
        "final_int_jw_ew": final_jw_ew,
        "final_int_s_leak_abs": final_s_leak_abs,
        "pi0_transport_max_norm": pi0_transport_max_norm,
        "pi0_transport_max_abs_rate": pi0_transport_max_abs_rate,
        "pi0_transport_status": pi0_transport_status,
        "pix_transport_max_norm": pix_transport_max_norm,
        "pix_transport_max_abs_rate": pix_transport_max_abs_rate,
        "pix_transport_status": pix_transport_status,
        "piy_transport_max_norm": piy_transport_max_norm,
        "piy_transport_max_abs_rate": piy_transport_max_abs_rate,
        "piy_transport_status": piy_transport_status,
        "piz_transport_max_norm": piz_transport_max_norm,
        "piz_transport_max_abs_rate": piz_transport_max_abs_rate,
        "piz_transport_status": piz_transport_status,
        "piw_transport_max_norm": piw_transport_max_norm,
        "piw_transport_max_abs_rate": piw_transport_max_abs_rate,
        "piw_transport_status": piw_transport_status,
        "em_transport_status": em_transport_status,
        "status": "FAIL" if failures else "PASS",
        "failures": "|".join(failures) if failures else "none",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to athenaPK binary")
    parser.add_argument("--workdir", default=".", help="Run directory")
    parser.add_argument(
        "--output-dir",
        default="modes4d_em_conservative_outputs",
        help="Directory for copied history files (relative to workdir if not absolute)",
    )
    parser.add_argument(
        "--scalar-input",
        default=None,
        help="Scalar pulse input deck (default: inputs/em4d_scalar_pulse.in)",
    )
    parser.add_argument(
        "--harris-input",
        default=None,
        help="Harris full input deck (default: inputs/harris_4d_full.in)",
    )
    parser.add_argument(
        "--harris-tlim",
        type=float,
        default=2.0e-2,
        help="Short Harris runtime limit for this regression",
    )
    parser.add_argument(
        "--harris-nlim",
        type=int,
        default=200,
        help="Short Harris cycle limit for this regression",
    )
    parser.add_argument("--min-scalar-mixed-ew2", type=float, default=1.0e-8)
    parser.add_argument("--min-scalar-mixed-c2", type=float, default=1.0e-7)
    parser.add_argument("--min-scalar-em-a2-mode1", type=float, default=1.0e-8)
    parser.add_argument("--min-scalar-pulse-xc-shift", type=float, default=0.0)
    parser.add_argument("--min-harris-mixed-ew2", type=float, default=1.0e-6)
    parser.add_argument("--min-harris-mixed-c2", type=float, default=1.0e-7)
    parser.add_argument("--transport-norm-tol", type=float, default=5.0e-2)
    parser.add_argument("--transport-abs-rate-tol", type=float, default=1.0e-8)
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]

    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    if args.scalar_input is None:
        scalar_input = repo_root / "inputs" / "em4d_scalar_pulse.in"
    else:
        scalar_input = Path(args.scalar_input).resolve()
    if args.harris_input is None:
        harris_input = repo_root / "inputs" / "harris_4d_full.in"
    else:
        harris_input = Path(args.harris_input).resolve()

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = workdir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    scalar_cols = run_case(
        binary=binary,
        input_path=scalar_input,
        workdir=workdir,
        output_hst=output_dir / "scalar_conservative.out.hst",
        extra_args=["modes4d/em_conservative_transport=true"],
    )
    harris_cols = run_case(
        binary=binary,
        input_path=harris_input,
        workdir=workdir,
        output_hst=output_dir / "harris_short_conservative.out.hst",
        extra_args=[
            "modes4d/em_conservative_transport=true",
            f"parthenon/time/tlim={args.harris_tlim}",
            f"parthenon/time/nlim={args.harris_nlim}",
        ],
    )

    scalar_result = analyze_scalar(
        scalar_cols,
        args.min_scalar_mixed_ew2,
        args.min_scalar_mixed_c2,
        args.min_scalar_em_a2_mode1,
        args.min_scalar_pulse_xc_shift,
    )
    harris_result = analyze_harris(
        harris_cols,
        args.min_harris_mixed_ew2,
        args.min_harris_mixed_c2,
        args.transport_norm_tol,
        args.transport_abs_rate_tol,
    )

    print(
        "case,max_mixed_ew2,max_mixed_c2,max_em_a2_mode_1,pulse_xc_shift,max_cont_mode0_abs,final_int_jw_ew,final_int_s_leak_abs,pi0_transport_max_norm,pi0_transport_max_abs_rate,pi0_transport_status,pix_transport_max_norm,pix_transport_max_abs_rate,pix_transport_status,piy_transport_max_norm,piy_transport_max_abs_rate,piy_transport_status,piz_transport_max_norm,piz_transport_max_abs_rate,piz_transport_status,piw_transport_max_norm,piw_transport_max_abs_rate,piw_transport_status,em_transport_status,status,failures"
    )
    print(
        ",".join(
            [
                scalar_result["case"],
                fmt(scalar_result["max_mixed_ew2"]),
                fmt(scalar_result["max_mixed_c2"]),
                fmt(scalar_result["max_em_a2_mode_1"]),
                fmt(scalar_result["pulse_xc_shift"]),
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "N/A",
                "nan",
                "nan",
                "N/A",
                "nan",
                "nan",
                "N/A",
                "nan",
                "nan",
                "N/A",
                "nan",
                "nan",
                "N/A",
                "N/A",
                scalar_result["status"],
                scalar_result["failures"],
            ]
        )
    )
    print(
        ",".join(
            [
                harris_result["case"],
                fmt(harris_result["max_mixed_ew2"]),
                fmt(harris_result["max_mixed_c2"]),
                "nan",
                "nan",
                fmt(harris_result["max_cont_mode0_abs"]),
                fmt(harris_result["final_int_jw_ew"]),
                fmt(harris_result["final_int_s_leak_abs"]),
                fmt(harris_result["pi0_transport_max_norm"]),
                fmt(harris_result["pi0_transport_max_abs_rate"]),
                harris_result["pi0_transport_status"],
                fmt(harris_result["pix_transport_max_norm"]),
                fmt(harris_result["pix_transport_max_abs_rate"]),
                harris_result["pix_transport_status"],
                fmt(harris_result["piy_transport_max_norm"]),
                fmt(harris_result["piy_transport_max_abs_rate"]),
                harris_result["piy_transport_status"],
                fmt(harris_result["piz_transport_max_norm"]),
                fmt(harris_result["piz_transport_max_abs_rate"]),
                harris_result["piz_transport_status"],
                fmt(harris_result["piw_transport_max_norm"]),
                fmt(harris_result["piw_transport_max_abs_rate"]),
                harris_result["piw_transport_status"],
                harris_result["em_transport_status"],
                harris_result["status"],
                harris_result["failures"],
            ]
        )
    )

    failures = []
    if scalar_result["status"] != "PASS":
        failures.append(f"scalar:{scalar_result['failures']}")
    if harris_result["status"] != "PASS":
        failures.append(f"harris:{harris_result['failures']}")
    if failures:
        raise SystemExit("conservative EM regression failed: " + ", ".join(failures))


if __name__ == "__main__":
    main()
