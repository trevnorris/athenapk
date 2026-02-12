#!/usr/bin/env python3
"""Windowed residual attribution for precomputed Harris scan .hst outputs."""

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


def parse_window_fractions(raw: str):
    vals = [float(x.strip()) for x in raw.split(",") if x.strip()]
    if len(vals) < 2:
        raise ValueError("Need at least two fractions.")
    if abs(vals[0]) > 1.0e-12 or abs(vals[-1] - 1.0) > 1.0e-12:
        raise ValueError("Fractions must start at 0 and end at 1.")
    for i in range(1, len(vals)):
        if vals[i] <= vals[i - 1]:
            raise ValueError("Fractions must be strictly increasing.")
    return vals


def fmt(value):
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    max_abs = max(abs(v) for v in finite)
    rms_abs = math.sqrt(sum(v * v for v in finite) / len(finite))
    final = finite[-1]
    return (max_abs, rms_abs, final)


def residual_rate_series(
    times,
    em_u_bulk,
    int_ja_ea,
    int_jw_ew,
    int_bridge_power=None,
    bridge_sign=0.0,
):
    if (
        times is None
        or em_u_bulk is None
        or int_ja_ea is None
        or int_jw_ew is None
        or len(times) < 2
        or len(em_u_bulk) != len(times)
        or len(int_ja_ea) != len(times)
        or len(int_jw_ew) != len(times)
    ):
        return ([], [])
    if int_bridge_power is not None and len(int_bridge_power) != len(times):
        return ([], [])

    out_t = []
    out_r = []
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        du_dt = (em_u_bulk[i] - em_u_bulk[i - 1]) / dt
        dja_dt = (int_ja_ea[i] - int_ja_ea[i - 1]) / dt
        djw_dt = (int_jw_ew[i] - int_jw_ew[i - 1]) / dt
        dbridge_dt = 0.0
        if int_bridge_power is not None:
            dbridge_dt = (int_bridge_power[i] - int_bridge_power[i - 1]) / dt
        out_t.append(times[i])
        out_r.append(du_dt + dja_dt + djw_dt + bridge_sign * dbridge_dt)
    return (out_t, out_r)


def select_window(times, values, t_lo, t_hi, include_left=True, include_right=False):
    out = []
    for t, v in zip(times, values):
        ge = (t >= t_lo) if include_left else (t > t_lo)
        le = (t <= t_hi) if include_right else (t < t_hi)
        if ge and le:
            out.append(v)
    return out


def compute_window_rows(case_dir, case_name, cols, gain_tag, gain, fractions):
    times = cols.get("time", [])
    if len(times) < 2:
        return []
    t0 = times[0]
    t1 = times[-1]
    span = t1 - t0
    if span <= 0.0:
        return []

    cont_mode0_max_abs = cols.get("m4d_cont_mode0_max_abs")
    cont_times = times if cont_mode0_max_abs is not None else []

    int_bridge_mode0 = cols.get("m4d_int_bridge_power_mode0")
    int_bridge_mode2 = cols.get("m4d_int_bridge_power_mode2")
    int_bridge_sum = cols.get("m4d_int_bridge_power_sum")
    if (
        int_bridge_sum is None
        and int_bridge_mode0 is not None
        and int_bridge_mode2 is not None
        and len(int_bridge_mode0) == len(int_bridge_mode2)
    ):
        int_bridge_sum = [t0 + t2 for t0, t2 in zip(int_bridge_mode0, int_bridge_mode2)]
    led_t_baseline, led_r_baseline = residual_rate_series(
        times,
        cols.get("m4d_em_u_bulk"),
        cols.get("m4d_int_ja_ea"),
        cols.get("m4d_int_jw_ew"),
        int_bridge_power=None,
        bridge_sign=0.0,
    )
    led_t_plus, led_r_plus = residual_rate_series(
        times,
        cols.get("m4d_em_u_bulk"),
        cols.get("m4d_int_ja_ea"),
        cols.get("m4d_int_jw_ew"),
        int_bridge_power=int_bridge_mode0,
        bridge_sign=1.0,
    )
    led_t_minus, led_r_minus = residual_rate_series(
        times,
        cols.get("m4d_em_u_bulk"),
        cols.get("m4d_int_ja_ea"),
        cols.get("m4d_int_jw_ew"),
        int_bridge_power=int_bridge_mode0,
        bridge_sign=-1.0,
    )
    led_t_sum, led_r_sum = residual_rate_series(
        times,
        cols.get("m4d_em_u_bulk"),
        cols.get("m4d_int_ja_ea"),
        cols.get("m4d_int_jw_ew"),
        int_bridge_power=int_bridge_sum,
        bridge_sign=1.0,
    )

    rows = []
    for idx in range(len(fractions) - 1):
        f_lo = fractions[idx]
        f_hi = fractions[idx + 1]
        w_name = "early" if idx == 0 else ("late" if idx == len(fractions) - 2 else "mid")
        t_lo = t0 + f_lo * span
        t_hi = t0 + f_hi * span
        include_right = idx == (len(fractions) - 2)

        cont_vals = (
            select_window(cont_times, cont_mode0_max_abs, t_lo, t_hi, True, include_right)
            if cont_mode0_max_abs is not None
            else []
        )
        cont_max, cont_rms, cont_final = finite_stats(cont_vals)

        base_vals = select_window(led_t_baseline, led_r_baseline, t_lo, t_hi, True, include_right)
        plus_vals = select_window(led_t_plus, led_r_plus, t_lo, t_hi, True, include_right)
        minus_vals = select_window(led_t_minus, led_r_minus, t_lo, t_hi, True, include_right)
        sum_vals = select_window(led_t_sum, led_r_sum, t_lo, t_hi, True, include_right)
        base_max, base_rms, base_final = finite_stats(base_vals)
        plus_max, plus_rms, plus_final = finite_stats(plus_vals)
        minus_max, minus_rms, minus_final = finite_stats(minus_vals)
        sum_max, sum_rms, sum_final = finite_stats(sum_vals)

        bridge_mode0_abs_vals = (
            [abs(v) for v in select_window(times, int_bridge_mode0, t_lo, t_hi, True, include_right)]
            if int_bridge_mode0 is not None
            else []
        )
        bridge_mode2_abs_vals = (
            [abs(v) for v in select_window(times, int_bridge_mode2, t_lo, t_hi, True, include_right)]
            if int_bridge_mode2 is not None
            else []
        )
        bridge_sum_abs_vals = (
            [abs(v) for v in select_window(times, int_bridge_sum, t_lo, t_hi, True, include_right)]
            if int_bridge_sum is not None
            else []
        )
        bridge_mode0_abs_max, bridge_mode0_abs_rms, bridge_mode0_abs_final = finite_stats(
            bridge_mode0_abs_vals
        )
        bridge_mode2_abs_max, bridge_mode2_abs_rms, bridge_mode2_abs_final = finite_stats(
            bridge_mode2_abs_vals
        )
        bridge_sum_abs_max, bridge_sum_abs_rms, bridge_sum_abs_final = finite_stats(
            bridge_sum_abs_vals
        )
        bridge_cancel_ratio_vals = []
        if (
            int_bridge_mode0 is not None
            and int_bridge_mode2 is not None
            and len(int_bridge_mode0) == len(times)
            and len(int_bridge_mode2) == len(times)
        ):
            mode0_window = select_window(
                times, int_bridge_mode0, t_lo, t_hi, True, include_right
            )
            mode2_window = select_window(
                times, int_bridge_mode2, t_lo, t_hi, True, include_right
            )
            for t0, t2 in zip(mode0_window, mode2_window):
                denom = abs(t0) + abs(t2)
                bridge_cancel_ratio_vals.append(abs(t0 + t2) / denom if denom > 0.0 else math.nan)
        (
            bridge_cancel_ratio_max,
            bridge_cancel_ratio_rms,
            bridge_cancel_ratio_final,
        ) = finite_stats(bridge_cancel_ratio_vals)

        candidates = []
        if math.isfinite(base_max):
            candidates.append(("baseline", base_max, base_rms, base_final))
        if math.isfinite(plus_max):
            candidates.append(("plus", plus_max, plus_rms, plus_final))
        if math.isfinite(minus_max):
            candidates.append(("minus", minus_max, minus_rms, minus_final))
        if math.isfinite(sum_max):
            candidates.append(("sum", sum_max, sum_rms, sum_final))
        if candidates:
            best = min(candidates, key=lambda x: x[1])
            best_label, best_max, best_rms, best_final = best
        else:
            best_label, best_max, best_rms, best_final = ("N/A", math.nan, math.nan, math.nan)

        rows.append(
            {
                "case_dir": str(case_dir),
                "gain_tag": gain_tag,
                "gain": gain,
                "case": case_name,
                "window": w_name,
                "t_start": t_lo,
                "t_end": t_hi,
                "cont_mode0_max_abs_max": cont_max,
                "cont_mode0_max_abs_rms": cont_rms,
                "cont_mode0_max_abs_final": cont_final,
                "ledger_baseline_max_abs_rate": base_max,
                "ledger_baseline_rms_abs_rate": base_rms,
                "ledger_baseline_final_rate": base_final,
                "ledger_plus_max_abs_rate": plus_max,
                "ledger_plus_rms_abs_rate": plus_rms,
                "ledger_plus_final_rate": plus_final,
                "ledger_minus_max_abs_rate": minus_max,
                "ledger_minus_rms_abs_rate": minus_rms,
                "ledger_minus_final_rate": minus_final,
                "ledger_sum_max_abs_rate": sum_max,
                "ledger_sum_rms_abs_rate": sum_rms,
                "ledger_sum_final_rate": sum_final,
                "ledger_best_label": best_label,
                "ledger_best_max_abs_rate": best_max,
                "ledger_best_rms_abs_rate": best_rms,
                "ledger_best_final_rate": best_final,
                "bridge_mode0_abs_max": bridge_mode0_abs_max,
                "bridge_mode0_abs_rms": bridge_mode0_abs_rms,
                "bridge_mode0_abs_final": bridge_mode0_abs_final,
                "bridge_mode2_abs_max": bridge_mode2_abs_max,
                "bridge_mode2_abs_rms": bridge_mode2_abs_rms,
                "bridge_mode2_abs_final": bridge_mode2_abs_final,
                "bridge_sum_abs_max": bridge_sum_abs_max,
                "bridge_sum_abs_rms": bridge_sum_abs_rms,
                "bridge_sum_abs_final": bridge_sum_abs_final,
                "bridge_cancel_ratio_max": bridge_cancel_ratio_max,
                "bridge_cancel_ratio_rms": bridge_cancel_ratio_rms,
                "bridge_cancel_ratio_final": bridge_cancel_ratio_final,
                "full_over_controlled_psi0_ratio": math.nan,
                "full_over_controlled_psiproj_ratio": math.nan,
                "full_over_controlled_psiw0_ratio": math.nan,
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        default="/projects/fluid-engine/out/step64_window_g*",
        help="Glob for case directories containing controlled/full .hst outputs",
    )
    parser.add_argument(
        "--window-fractions",
        default="0,0.333333,0.666667,1",
        help="Comma-separated monotonically increasing fractional boundaries from 0 to 1",
    )
    parser.add_argument(
        "--out-csv",
        default="/projects/fluid-engine/out/step64_windowed_residuals.csv",
        help="Path to write windowed residual CSV",
    )
    parser.add_argument(
        "--out-md",
        default="/projects/fluid-engine/out/step64_windowed_residuals.md",
        help="Path to write compact markdown report",
    )
    args = parser.parse_args()

    fractions = parse_window_fractions(args.window_fractions)
    script_dir = Path(__file__).resolve().parent
    hsm = load_harris_scan_matrix_module(script_dir)

    dirs = [Path(p) for p in sorted(glob.glob(args.glob))]
    if not dirs:
        raise SystemExit(f"No directories matched: {args.glob}")

    rows = []
    full_final_by_dir = {}
    ctrl_final_by_dir = {}
    for d in dirs:
        controlled_hst = d / "harris_controlled.out1.hst"
        full_hst = d / "harris_full.out1.hst"
        if not controlled_hst.exists() or not full_hst.exists():
            continue
        gain_tag = d.name.split("_g", 1)[-1]
        gain = parse_gain_tag(gain_tag)

        controlled_cols = hsm.parse_hst(controlled_hst)
        full_cols = hsm.parse_hst(full_hst)

        ctrl_final_by_dir[str(d)] = {
            "psi0": controlled_cols["m4d_psi0_span"][-1],
            "psiproj": controlled_cols.get("m4d_psi_proj_span", controlled_cols["m4d_psi0_span"])[-1],
            "psiw0": controlled_cols.get("m4d_psi_w0_span", controlled_cols["m4d_psi0_span"])[-1],
        }
        full_final_by_dir[str(d)] = {
            "psi0": full_cols["m4d_psi0_span"][-1],
            "psiproj": full_cols.get("m4d_psi_proj_span", full_cols["m4d_psi0_span"])[-1],
            "psiw0": full_cols.get("m4d_psi_w0_span", full_cols["m4d_psi0_span"])[-1],
        }

        rows.extend(compute_window_rows(d, "controlled", controlled_cols, gain_tag, gain, fractions))
        rows.extend(compute_window_rows(d, "full", full_cols, gain_tag, gain, fractions))

    # Populate ratios on full rows.
    for r in rows:
        if r["case"] != "full":
            continue
        dkey = r["case_dir"]
        f = full_final_by_dir.get(dkey, {})
        c = ctrl_final_by_dir.get(dkey, {})
        psi0 = f.get("psi0", math.nan)
        p0 = c.get("psi0", math.nan)
        psip = f.get("psiproj", math.nan)
        pp = c.get("psiproj", math.nan)
        psiw = f.get("psiw0", math.nan)
        pw = c.get("psiw0", math.nan)
        r["full_over_controlled_psi0_ratio"] = (
            psi0 / p0 if math.isfinite(psi0) and math.isfinite(p0) and abs(p0) > 0.0 else math.nan
        )
        r["full_over_controlled_psiproj_ratio"] = (
            psip / pp if math.isfinite(psip) and math.isfinite(pp) and abs(pp) > 0.0 else math.nan
        )
        r["full_over_controlled_psiw0_ratio"] = (
            psiw / pw if math.isfinite(psiw) and math.isfinite(pw) and abs(pw) > 0.0 else math.nan
        )

    rows.sort(
        key=lambda r: (
            (math.inf if not math.isfinite(r["gain"]) else r["gain"]),
            r["case_dir"],
            0 if r["case"] == "controlled" else 1,
            {"early": 0, "mid": 1, "late": 2}.get(r["window"], 9),
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
        f.write("# Windowed Residual Attribution\n\n")
        f.write(f"Parsed `{len(rows)}` rows from glob `{args.glob}`.\n\n")
        full_rows = [r for r in rows if r["case"] == "full"]
        if not full_rows:
            f.write("No full-case rows found.\n")
        else:
            f.write(
                "| gain | window | best ledger | best max | bridge sum max | cancel ratio | cont max | psi0 ratio | psiw0 ratio |\n"
            )
            f.write("|---:|:---|:---|---:|---:|---:|---:|---:|---:|\n")
            for r in full_rows:
                f.write(
                    "| "
                    + " | ".join(
                        [
                            fmt(r["gain"]),
                            r["window"],
                            r["ledger_best_label"],
                            fmt(r["ledger_best_max_abs_rate"]),
                            fmt(r["bridge_sum_abs_max"]),
                            fmt(r["bridge_cancel_ratio_max"]),
                            fmt(r["cont_mode0_max_abs_max"]),
                            fmt(r["full_over_controlled_psi0_ratio"]),
                            fmt(r["full_over_controlled_psiw0_ratio"]),
                        ]
                    )
                    + " |\n"
                )
            f.write("\n")
            f.write("## Worst Early/Mid/Late Full-Case Ledger Residual by Gain\n\n")
            by_gain = {}
            for r in full_rows:
                g = r["gain"]
                by_gain.setdefault(g, []).append(r)
            for g in sorted(by_gain.keys()):
                grp = by_gain[g]
                worst = max(
                    grp,
                    key=lambda x: (
                        -math.inf
                        if not math.isfinite(x["ledger_best_max_abs_rate"])
                        else x["ledger_best_max_abs_rate"]
                    ),
                )
                f.write(
                    f"- gain={fmt(g)}: worst_window={worst['window']}, "
                    f"ledger_best={fmt(worst['ledger_best_max_abs_rate'])}, "
                    f"label={worst['ledger_best_label']}, "
                    f"psi0_ratio={fmt(worst['full_over_controlled_psi0_ratio'])}\n"
                )

    print(f"summary_csv,{out_csv}")
    print(f"summary_md,{out_md}")
    print(f"rows,{len(rows)}")


if __name__ == "__main__":
    main()
