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
    header_index = -1
    for idx, line in enumerate(lines):
        if line.startswith("# [1]="):
            header = line
            header_index = idx
            break
    if header is None:
        raise RuntimeError(f"Could not find history header in {path}")

    labels = []
    for token in header.replace("#", "").split():
        if token.startswith("[") and "]=" in token:
            labels.append(token.split("]=", 1)[1])

    rows = []
    for line in lines[header_index + 1 :]:
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != len(labels):
            continue
        rows.append([float(x) for x in parts])

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


def combine_abs_series(series_list):
    length = None
    for series in series_list:
        if series is None:
            continue
        if length is None:
            length = len(series)
        elif len(series) != length:
            return None
    if length is None:
        return None

    combined = []
    for i in range(length):
        total = 0.0
        any_finite = False
        for series in series_list:
            if series is None:
                continue
            value = series[i]
            if math.isfinite(value):
                total += abs(value)
                any_finite = True
        combined.append(total if any_finite else math.nan)
    return combined


def mode_activity_proxy(jw_mode_l2, jw_ew):
    candidates = []
    if isinstance(jw_mode_l2, float) and math.isfinite(jw_mode_l2):
        candidates.append(jw_mode_l2)
    if isinstance(jw_ew, float) and math.isfinite(jw_ew):
        candidates.append(abs(jw_ew))
    if not candidates:
        return math.nan
    return max(candidates)


def resolve_input_path(raw_value, repo_root, workdir):
    path = Path(raw_value)
    if path.is_absolute():
        return path.resolve()
    for base in (Path.cwd(), workdir, repo_root):
        candidate = (base / path).resolve()
        if candidate.exists():
            return candidate
    return (repo_root / path).resolve()


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
        or len(quantity_mode0) != len(times)
        or len(int_divquantity_mode0) != len(times)
    ):
        return (math.nan, math.nan, math.nan, math.nan)
    if int_srcquantity_mode0 is not None and len(int_srcquantity_mode0) != len(times):
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
        dsrc_dt = 0.0
        if int_srcquantity_mode0 is not None:
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


def em_bulk_ledger_abs_rate_metrics(
    times, em_u_bulk, int_ja_ea, int_jw_ew, skip_initial_interval=True
):
    if (
        len(times) < 2
        or em_u_bulk is None
        or int_ja_ea is None
        or int_jw_ew is None
        or len(em_u_bulk) != len(times)
        or len(int_ja_ea) != len(times)
        or len(int_jw_ew) != len(times)
    ):
        return (math.nan, math.nan, math.nan)

    abs_rates = []
    final_rate = math.nan
    for i in range(1, len(times)):
        if skip_initial_interval and i == 1:
            continue
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        du_dt = (em_u_bulk[i] - em_u_bulk[i - 1]) / dt
        dja_dt = (int_ja_ea[i] - int_ja_ea[i - 1]) / dt
        djw_dt = (int_jw_ew[i] - int_jw_ew[i - 1]) / dt
        rate = du_dt + dja_dt + djw_dt
        abs_rates.append(abs(rate))
        final_rate = rate

    if not abs_rates:
        return (math.nan, math.nan, final_rate)

    max_abs = max(abs_rates)
    rms_abs = math.sqrt(sum(x * x for x in abs_rates) / len(abs_rates))
    return (max_abs, rms_abs, final_rate)


def reconnection_rate_metrics(times, psi, skip_initial_interval=True):
    if len(times) < 2 or len(psi) != len(times):
        return (math.nan, math.nan, math.nan, math.nan, math.nan, math.nan)

    max_rate = -math.inf
    min_rate = math.inf
    max_abs_rate = -math.inf
    final_rate = math.nan
    time_at_max_rate = math.nan
    time_at_max_abs_rate = math.nan

    for i in range(1, len(times)):
        if skip_initial_interval and i == 1:
            continue
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        rate = (psi[i] - psi[i - 1]) / dt
        final_rate = rate
        if rate > max_rate:
            max_rate = rate
            time_at_max_rate = times[i]
        if rate < min_rate:
            min_rate = rate
        if abs(rate) > max_abs_rate:
            max_abs_rate = abs(rate)
            time_at_max_abs_rate = times[i]

    if max_rate == -math.inf:
        return (math.nan, math.nan, math.nan, final_rate, math.nan, math.nan)

    return (
        max_rate,
        min_rate,
        max_abs_rate,
        final_rate,
        time_at_max_rate,
        time_at_max_abs_rate,
    )


def max_abs_finite(series):
    values = [abs(v) for v in series if math.isfinite(v)]
    if not values:
        return math.nan
    return max(values)


def first_crossing_time(times, series, abs_threshold):
    if abs_threshold <= 0.0 or len(times) != len(series):
        return math.nan
    for t, value in zip(times, series):
        if math.isfinite(value) and abs(value) >= abs_threshold:
            return t
    return math.nan


def local_peak_count(
    times,
    series,
    min_abs_value,
    min_prominence,
    t_min=-math.inf,
    t_max=math.inf,
):
    if len(times) < 3 or len(series) != len(times):
        return 0
    count = 0
    for i in range(1, len(series) - 1):
        t = times[i]
        if t < t_min or t > t_max:
            continue
        left = series[i - 1]
        center = series[i]
        right = series[i + 1]
        if (
            not math.isfinite(left)
            or not math.isfinite(center)
            or not math.isfinite(right)
        ):
            continue
        if center > left and center >= right:
            prominence = center - max(left, right)
            if abs(center) >= min_abs_value and prominence >= min_prominence:
                count += 1
    return count


def analyze_case(
    case_name,
    cols,
    closure_norm_tol,
    closure_abs_rate_tol,
    local_mode0_abs_rate_tol,
    em_bulk_ledger_abs_rate_tol,
    cycle_onset_psi0_abs_threshold,
    cycle_peak_min_amplitude_rel,
    cycle_peak_min_amplitude_abs,
    cycle_peak_min_prominence_rel,
    cycle_peak_min_prominence_abs,
):
    times = cols["time"]
    psi0 = cols["m4d_psi0_span"]
    psi_proj = maybe_col(cols, "m4d_psi_proj_span")
    if psi_proj is None or len(psi_proj) != len(psi0):
        psi_proj = psi0
    psi_proj_matched = maybe_col(cols, "m4d_psi_proj_matched_span")
    if psi_proj_matched is None or len(psi_proj_matched) != len(psi0):
        psi_proj_matched = psi_proj
    psi_proj_point = maybe_col(cols, "m4d_psi_proj_point_span")
    if psi_proj_point is None or len(psi_proj_point) != len(psi0):
        psi_proj_point = psi_proj
    psi_proj_gaussian = maybe_col(cols, "m4d_psi_proj_gaussian_span")
    if psi_proj_gaussian is None or len(psi_proj_gaussian) != len(psi0):
        psi_proj_gaussian = psi_proj
    psi_w0 = maybe_col(cols, "m4d_psi_w0_span")
    if psi_w0 is None or len(psi_w0) != len(psi0):
        psi_w0 = psi_proj
    leak = cols["m4d_int_s_leak"]
    leak_abs = cols.get("m4d_int_s_leak_abs", leak)
    jw_ew = cols["m4d_int_jw_ew"]
    ja_ea = maybe_col(cols, "m4d_int_ja_ea")
    epar2 = cols["m4d_brane_epar2"]
    mixed_ew2 = cols["m4d_mixed_ew2"]
    mixed_c2 = cols["m4d_mixed_c2"]
    em_u_bulk = maybe_col(cols, "m4d_em_u_bulk")
    em_u_resolved = maybe_col(cols, "m4d_em_u_resolved")
    em_sw = maybe_col(cols, "m4d_em_sw")
    em_leak_w = maybe_col(cols, "m4d_em_leak_w")
    emf_vw_c_abs = maybe_col(cols, "m4d_emf_vw_c_abs")
    emf_vw_c_par_abs = maybe_col(cols, "m4d_emf_vw_c_par_abs")
    emf_cov_vxb_abs = maybe_col(cols, "m4d_emf_cov_vxb_abs")
    emf_cov_vxb_par_abs = maybe_col(cols, "m4d_emf_cov_vxb_par_abs")
    jw_l1 = maybe_col(cols, "m4d_jw_l1")
    jw_l2 = maybe_col(cols, "m4d_jw_l2")
    ew_l1 = maybe_col(cols, "m4d_ew_l1")
    ew_l2 = maybe_col(cols, "m4d_ew_l2")
    helicity_sub = maybe_col(cols, "m4d_helicity_sub")
    edotb_sub = maybe_col(cols, "m4d_edotb_sub")
    em_a2_mode1 = maybe_col(cols, "m4d_em_a2_mode_1")
    em_pi2_mode1 = maybe_col(cols, "m4d_em_pi2_mode_1")
    jw_mode1_l2 = maybe_col(cols, "m4d_jw_mode_l2_1")
    jw_mode1 = maybe_col(cols, "m4d_jw_mode_1")
    charge_mode0 = maybe_col(cols, "m4d_charge_mode_0")
    momx_mode0 = maybe_col(cols, "m4d_momx_mode_0")
    momy_mode0 = maybe_col(cols, "m4d_momy_mode_0")
    momz_mode0 = maybe_col(cols, "m4d_momz_mode_0")
    momw_mode0 = maybe_col(cols, "m4d_momw_mode_0")
    energy_mode0 = maybe_col(cols, "m4d_energy_mode_0")
    pi0_mode0 = maybe_col(cols, "m4d_pi0_mode_0")
    pix_mode0 = maybe_col(cols, "m4d_pix_mode_0")
    piy_mode0 = maybe_col(cols, "m4d_piy_mode_0")
    piz_mode0 = maybe_col(cols, "m4d_piz_mode_0")
    piw_mode0 = maybe_col(cols, "m4d_piw_mode_0")
    int_divj_mode0 = maybe_col(cols, "m4d_int_divj_mode0")
    int_divmomx_mode0 = maybe_col(cols, "m4d_int_divmomx_mode0")
    int_divmomy_mode0 = maybe_col(cols, "m4d_int_divmomy_mode0")
    int_divmomz_mode0 = maybe_col(cols, "m4d_int_divmomz_mode0")
    int_divmomw_mode0 = maybe_col(cols, "m4d_int_divmomw_mode0")
    int_divenergy_mode0 = maybe_col(cols, "m4d_int_divenergy_mode0")
    int_divpi0_mode0 = maybe_col(cols, "m4d_int_divpi0_mode0")
    int_divpix_mode0 = maybe_col(cols, "m4d_int_divpix_mode0")
    int_divpiy_mode0 = maybe_col(cols, "m4d_int_divpiy_mode0")
    int_divpiz_mode0 = maybe_col(cols, "m4d_int_divpiz_mode0")
    int_divpiw_mode0 = maybe_col(cols, "m4d_int_divpiw_mode0")
    int_srcmomx_mode0 = maybe_col(cols, "m4d_int_srcmomx_mode0")
    int_srcmomy_mode0 = maybe_col(cols, "m4d_int_srcmomy_mode0")
    int_srcmomz_mode0 = maybe_col(cols, "m4d_int_srcmomz_mode0")
    int_srcmomw_mode0 = maybe_col(cols, "m4d_int_srcmomw_mode0")
    int_srcenergy_mode0 = maybe_col(cols, "m4d_int_srcenergy_mode0")
    int_srcpi0_mode0 = maybe_col(cols, "m4d_int_srcpi0_mode0")
    int_srcpix_mode0 = maybe_col(cols, "m4d_int_srcpix_mode0")
    int_srcpiy_mode0 = maybe_col(cols, "m4d_int_srcpiy_mode0")
    int_srcpiz_mode0 = maybe_col(cols, "m4d_int_srcpiz_mode0")
    int_srcpiw_mode0 = maybe_col(cols, "m4d_int_srcpiw_mode0")
    int_src_em_laplacian_abs = maybe_col(cols, "m4d_int_src_em_laplacian_abs")
    int_src_em_mass_abs = maybe_col(cols, "m4d_int_src_em_mass_abs")
    int_src_em_current_abs = maybe_col(cols, "m4d_int_src_em_current_abs")
    int_src_em_damping_abs = maybe_col(cols, "m4d_int_src_em_damping_abs")
    int_src_em_spatial_mixed_abs = maybe_col(cols, "m4d_int_src_em_spatial_mixed_abs")
    int_src_em_timelike_abs = maybe_col(cols, "m4d_int_src_em_timelike_abs")
    int_src_em_gauge_abs = maybe_col(cols, "m4d_int_src_em_gauge_abs")
    int_div_mode0_plasma_abs = maybe_col(cols, "m4d_int_div_mode0_plasma_abs")
    int_div_mode0_em_abs = maybe_col(cols, "m4d_int_div_mode0_em_abs")
    int_div_mode0_total_abs = maybe_col(cols, "m4d_int_div_mode0_total_abs")
    int_src_mode0_plasma_abs = maybe_col(cols, "m4d_int_src_mode0_plasma_abs")
    int_src_mode0_em_abs = maybe_col(cols, "m4d_int_src_mode0_em_abs")
    int_src_mode0_timelike_abs = maybe_col(cols, "m4d_int_src_mode0_timelike_abs")
    int_src_mode0_total_abs = maybe_col(cols, "m4d_int_src_mode0_total_abs")
    cont_local_mode0_max_abs = maybe_col(cols, "m4d_cont_mode0_max_abs")
    cont_local_mode0_l1 = maybe_col(cols, "m4d_cont_mode0_l1")
    cont_local_mode0_l2 = maybe_col(cols, "m4d_cont_mode0_l2")
    gauge_l1 = maybe_col(cols, "m4d_gauge_l1")
    gauge_l2 = maybe_col(cols, "m4d_gauge_l2")
    gauge_max_abs = maybe_col(cols, "m4d_gauge_max_abs")
    gauge_mode0_l1 = maybe_col(cols, "m4d_gauge_mode0_l1")
    gauge_mode0_l2 = maybe_col(cols, "m4d_gauge_mode0_l2")
    gauge_mode0_max_abs = maybe_col(cols, "m4d_gauge_mode0_max_abs")

    psi_w0_minus_psi_proj = [w - p for w, p in zip(psi_w0, psi_proj)]
    psi_w0_minus_psi0 = [w - p0 for w, p0 in zip(psi_w0, psi0)]
    psi_w0_over_psi_proj = [
        (w / p) if math.isfinite(w) and math.isfinite(p) and abs(p) > 0.0 else math.nan
        for w, p in zip(psi_w0, psi_proj)
    ]
    psi_proj_over_psi0 = [
        (p / p0) if math.isfinite(p) and math.isfinite(p0) and abs(p0) > 0.0 else math.nan
        for p, p0 in zip(psi_proj, psi0)
    ]
    psi_w0_over_psi0 = [
        (w / p0) if math.isfinite(w) and math.isfinite(p0) and abs(p0) > 0.0 else math.nan
        for w, p0 in zip(psi_w0, psi0)
    ]

    if int_div_mode0_plasma_abs is None:
        int_div_mode0_plasma_abs = combine_abs_series(
            [
                int_divmomx_mode0,
                int_divmomy_mode0,
                int_divmomz_mode0,
                int_divmomw_mode0,
                int_divenergy_mode0,
            ]
        )
    if int_div_mode0_em_abs is None:
        int_div_mode0_em_abs = combine_abs_series(
            [
                int_divpi0_mode0,
                int_divpix_mode0,
                int_divpiy_mode0,
                int_divpiz_mode0,
                int_divpiw_mode0,
            ]
        )
    if int_div_mode0_total_abs is None:
        int_div_mode0_total_abs = combine_abs_series(
            [
                int_divj_mode0,
                int_divmomx_mode0,
                int_divmomy_mode0,
                int_divmomz_mode0,
                int_divmomw_mode0,
                int_divenergy_mode0,
                int_divpi0_mode0,
                int_divpix_mode0,
                int_divpiy_mode0,
                int_divpiz_mode0,
                int_divpiw_mode0,
            ]
        )
    if int_src_mode0_plasma_abs is None:
        int_src_mode0_plasma_abs = combine_abs_series(
            [
                int_srcmomx_mode0,
                int_srcmomy_mode0,
                int_srcmomz_mode0,
                int_srcmomw_mode0,
                int_srcenergy_mode0,
            ]
        )
    if int_src_mode0_em_abs is None:
        int_src_mode0_em_abs = combine_abs_series(
            [
                int_srcpi0_mode0,
                int_srcpix_mode0,
                int_srcpiy_mode0,
                int_srcpiz_mode0,
                int_srcpiw_mode0,
            ]
        )
    if int_src_mode0_timelike_abs is None:
        int_src_mode0_timelike_abs = combine_abs_series(
            [maybe_col(cols, "m4d_int_src_timelike_a0_from_piw"), maybe_col(cols, "m4d_int_src_timelike_aw_from_pi0")]
        )
    if int_src_mode0_total_abs is None:
        int_src_mode0_total_abs = combine_abs_series(
            [
                int_srcmomx_mode0,
                int_srcmomy_mode0,
                int_srcmomz_mode0,
                int_srcmomw_mode0,
                int_srcenergy_mode0,
                int_srcpi0_mode0,
                int_srcpix_mode0,
                int_srcpiy_mode0,
                int_srcpiz_mode0,
                int_srcpiw_mode0,
                maybe_col(cols, "m4d_int_src_timelike_a0_from_piw"),
                maybe_col(cols, "m4d_int_src_timelike_aw_from_pi0"),
            ]
        )

    em_u_sub = None
    if (
        em_u_bulk is not None
        and em_u_resolved is not None
        and len(em_u_bulk) == len(em_u_resolved)
    ):
        em_u_sub = [ub - ur for ub, ur in zip(em_u_bulk, em_u_resolved)]

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

    (
        em_bulk_ledger_max_abs_rate,
        em_bulk_ledger_rms_abs_rate,
        em_bulk_ledger_final_rate,
    ) = em_bulk_ledger_abs_rate_metrics(times, em_u_bulk, ja_ea, jw_ew)
    em_bulk_ledger_status = (
        "PASS"
        if (
            not math.isnan(em_bulk_ledger_max_abs_rate)
            and em_bulk_ledger_max_abs_rate <= em_bulk_ledger_abs_rate_tol
        )
        else ("N/A" if math.isnan(em_bulk_ledger_max_abs_rate) else "FAIL")
    )
    (
        max_dpsi0_dt,
        min_dpsi0_dt,
        max_abs_dpsi0_dt,
        final_dpsi0_dt,
        time_at_max_dpsi0_dt,
        time_at_max_abs_dpsi0_dt,
    ) = reconnection_rate_metrics(times, psi0)
    (
        max_dpsi_proj_dt,
        min_dpsi_proj_dt,
        max_abs_dpsi_proj_dt,
        final_dpsi_proj_dt,
        time_at_max_dpsi_proj_dt,
        time_at_max_abs_dpsi_proj_dt,
    ) = reconnection_rate_metrics(times, psi_proj)
    (
        max_dpsi_w0_dt,
        min_dpsi_w0_dt,
        max_abs_dpsi_w0_dt,
        final_dpsi_w0_dt,
        time_at_max_dpsi_w0_dt,
        time_at_max_abs_dpsi_w0_dt,
    ) = reconnection_rate_metrics(times, psi_w0)
    time_first_abs_psi0_onset = first_crossing_time(
        times, psi0, cycle_onset_psi0_abs_threshold
    )
    psi_w0_peak_scale = max_abs_finite(psi_w0)
    psi_w0_peak_min_abs = cycle_peak_min_amplitude_abs
    psi_w0_peak_min_prom = cycle_peak_min_prominence_abs
    if math.isfinite(psi_w0_peak_scale):
        psi_w0_peak_min_abs = max(
            psi_w0_peak_min_abs, cycle_peak_min_amplitude_rel * psi_w0_peak_scale
        )
        psi_w0_peak_min_prom = max(
            psi_w0_peak_min_prom, cycle_peak_min_prominence_rel * psi_w0_peak_scale
        )
    src_mode0_total_peak_scale = max_abs_finite(int_src_mode0_total_abs or [])
    src_mode0_total_peak_min_abs = cycle_peak_min_amplitude_abs
    src_mode0_total_peak_min_prom = cycle_peak_min_prominence_abs
    if math.isfinite(src_mode0_total_peak_scale):
        src_mode0_total_peak_min_abs = max(
            src_mode0_total_peak_min_abs,
            cycle_peak_min_amplitude_rel * src_mode0_total_peak_scale,
        )
        src_mode0_total_peak_min_prom = max(
            src_mode0_total_peak_min_prom,
            cycle_peak_min_prominence_rel * src_mode0_total_peak_scale,
        )
    psi_w0_peak_count_total = local_peak_count(
        times, psi_w0, psi_w0_peak_min_abs, psi_w0_peak_min_prom
    )
    src_mode0_total_peak_count_total = (
        local_peak_count(
            times,
            int_src_mode0_total_abs,
            src_mode0_total_peak_min_abs,
            src_mode0_total_peak_min_prom,
        )
        if int_src_mode0_total_abs is not None
        else math.nan
    )
    if math.isfinite(time_first_abs_psi0_onset):
        pre_t_max = math.nextafter(time_first_abs_psi0_onset, -math.inf)
        psi_w0_peak_count_pre_onset = local_peak_count(
            times,
            psi_w0,
            psi_w0_peak_min_abs,
            psi_w0_peak_min_prom,
            t_max=pre_t_max,
        )
        psi_w0_peak_count_post_onset = local_peak_count(
            times,
            psi_w0,
            psi_w0_peak_min_abs,
            psi_w0_peak_min_prom,
            t_min=time_first_abs_psi0_onset,
        )
        if int_src_mode0_total_abs is not None:
            src_mode0_total_peak_count_pre_onset = local_peak_count(
                times,
                int_src_mode0_total_abs,
                src_mode0_total_peak_min_abs,
                src_mode0_total_peak_min_prom,
                t_max=pre_t_max,
            )
            src_mode0_total_peak_count_post_onset = local_peak_count(
                times,
                int_src_mode0_total_abs,
                src_mode0_total_peak_min_abs,
                src_mode0_total_peak_min_prom,
                t_min=time_first_abs_psi0_onset,
            )
        else:
            src_mode0_total_peak_count_pre_onset = math.nan
            src_mode0_total_peak_count_post_onset = math.nan
    else:
        psi_w0_peak_count_pre_onset = psi_w0_peak_count_total
        psi_w0_peak_count_post_onset = 0
        if int_src_mode0_total_abs is not None:
            src_mode0_total_peak_count_pre_onset = src_mode0_total_peak_count_total
            src_mode0_total_peak_count_post_onset = 0
        else:
            src_mode0_total_peak_count_pre_onset = math.nan
            src_mode0_total_peak_count_post_onset = math.nan

    (
        momx_transport_max_norm,
        momx_transport_rms_norm,
        momx_transport_max_abs_rate,
        momx_transport_final_rate,
    ) = transport_balance_metrics(
        times, momx_mode0, int_divmomx_mode0, int_srcmomx_mode0, closure_abs_rate_tol
    )
    (
        momy_transport_max_norm,
        momy_transport_rms_norm,
        momy_transport_max_abs_rate,
        momy_transport_final_rate,
    ) = transport_balance_metrics(
        times, momy_mode0, int_divmomy_mode0, int_srcmomy_mode0, closure_abs_rate_tol
    )
    (
        momz_transport_max_norm,
        momz_transport_rms_norm,
        momz_transport_max_abs_rate,
        momz_transport_final_rate,
    ) = transport_balance_metrics(
        times, momz_mode0, int_divmomz_mode0, int_srcmomz_mode0, closure_abs_rate_tol
    )
    (
        momw_transport_max_norm,
        momw_transport_rms_norm,
        momw_transport_max_abs_rate,
        momw_transport_final_rate,
    ) = transport_balance_metrics(
        times, momw_mode0, int_divmomw_mode0, int_srcmomw_mode0, closure_abs_rate_tol
    )
    (
        energy_transport_max_norm,
        energy_transport_rms_norm,
        energy_transport_max_abs_rate,
        energy_transport_final_rate,
    ) = transport_balance_metrics(
        times,
        energy_mode0,
        int_divenergy_mode0,
        int_srcenergy_mode0,
        closure_abs_rate_tol,
    )
    (
        pi0_transport_max_norm,
        pi0_transport_rms_norm,
        pi0_transport_max_abs_rate,
        pi0_transport_final_rate,
    ) = transport_balance_metrics(
        times, pi0_mode0, int_divpi0_mode0, int_srcpi0_mode0, closure_abs_rate_tol
    )
    (
        pix_transport_max_norm,
        pix_transport_rms_norm,
        pix_transport_max_abs_rate,
        pix_transport_final_rate,
    ) = transport_balance_metrics(
        times, pix_mode0, int_divpix_mode0, int_srcpix_mode0, closure_abs_rate_tol
    )
    (
        piy_transport_max_norm,
        piy_transport_rms_norm,
        piy_transport_max_abs_rate,
        piy_transport_final_rate,
    ) = transport_balance_metrics(
        times, piy_mode0, int_divpiy_mode0, int_srcpiy_mode0, closure_abs_rate_tol
    )
    (
        piz_transport_max_norm,
        piz_transport_rms_norm,
        piz_transport_max_abs_rate,
        piz_transport_final_rate,
    ) = transport_balance_metrics(
        times, piz_mode0, int_divpiz_mode0, int_srcpiz_mode0, closure_abs_rate_tol
    )
    (
        piw_transport_max_norm,
        piw_transport_rms_norm,
        piw_transport_max_abs_rate,
        piw_transport_final_rate,
    ) = transport_balance_metrics(
        times, piw_mode0, int_divpiw_mode0, int_srcpiw_mode0, closure_abs_rate_tol
    )

    def rate_status(max_norm, max_abs):
        if math.isnan(max_norm):
            return "N/A"
        return "PASS" if (max_norm <= closure_norm_tol or max_abs <= closure_abs_rate_tol) else "FAIL"

    momx_transport_status = rate_status(
        momx_transport_max_norm, momx_transport_max_abs_rate
    )
    momy_transport_status = rate_status(
        momy_transport_max_norm, momy_transport_max_abs_rate
    )
    momz_transport_status = rate_status(
        momz_transport_max_norm, momz_transport_max_abs_rate
    )
    momw_transport_status = rate_status(
        momw_transport_max_norm, momw_transport_max_abs_rate
    )
    energy_transport_status = rate_status(
        energy_transport_max_norm, energy_transport_max_abs_rate
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
    plasma_transport_statuses = [
        momx_transport_status,
        momy_transport_status,
        momz_transport_status,
        momw_transport_status,
        energy_transport_status,
    ]
    em_transport_statuses = [
        pi0_transport_status,
        pix_transport_status,
        piy_transport_status,
        piz_transport_status,
        piw_transport_status,
    ]
    em_transport_closure_status = (
        "PASS" if all(s in ("PASS", "N/A") for s in em_transport_statuses) else "FAIL"
    )
    transport_statuses = plasma_transport_statuses + em_transport_statuses
    transport_closure_status = (
        "PASS"
        if all(s in ("PASS", "N/A") for s in transport_statuses)
        else "FAIL"
    )

    result = {
        "case": case_name,
        "final_psi0_span": psi0[-1],
        "final_psi_proj_span": psi_proj[-1],
        "final_psi_proj_matched_span": psi_proj_matched[-1],
        "final_psi_proj_point_span": psi_proj_point[-1],
        "final_psi_proj_gaussian_span": psi_proj_gaussian[-1],
        "final_psi_w0_span": psi_w0[-1],
        "final_psi_w0_minus_psi_proj_span": psi_w0_minus_psi_proj[-1],
        "final_psi_w0_minus_psi0_span": psi_w0_minus_psi0[-1],
        "final_psi_w0_over_psi_proj_span": psi_w0_over_psi_proj[-1],
        "final_psi_proj_over_psi0_span": psi_proj_over_psi0[-1],
        "final_psi_w0_over_psi0_span": psi_w0_over_psi0[-1],
        "max_abs_psi_w0_minus_psi_proj_span": max(
            (abs(v) for v in psi_w0_minus_psi_proj if math.isfinite(v)),
            default=math.nan,
        ),
        "final_s_leak": leak[-1],
        "final_s_leak_abs": leak_abs[-1],
        "final_jw_ew": jw_ew[-1],
        "final_ja_ea": ja_ea[-1] if ja_ea is not None else math.nan,
        "final_brane_epar2": epar2[-1],
        "final_mixed_ew2": mixed_ew2[-1],
        "final_mixed_c2": mixed_c2[-1],
        "final_em_u_bulk": em_u_bulk[-1] if em_u_bulk is not None else math.nan,
        "final_em_u_resolved": em_u_resolved[-1]
        if em_u_resolved is not None
        else math.nan,
        "final_em_u_sub": em_u_sub[-1] if em_u_sub is not None else math.nan,
        "final_em_sw": em_sw[-1] if em_sw is not None else math.nan,
        "final_em_leak_w": em_leak_w[-1] if em_leak_w is not None else math.nan,
        "final_emf_vw_c_abs": emf_vw_c_abs[-1] if emf_vw_c_abs is not None else math.nan,
        "final_emf_vw_c_par_abs": emf_vw_c_par_abs[-1]
        if emf_vw_c_par_abs is not None
        else math.nan,
        "final_emf_cov_vxb_abs": emf_cov_vxb_abs[-1]
        if emf_cov_vxb_abs is not None
        else math.nan,
        "final_emf_cov_vxb_par_abs": emf_cov_vxb_par_abs[-1]
        if emf_cov_vxb_par_abs is not None
        else math.nan,
        "final_jw_l1": jw_l1[-1] if jw_l1 is not None else math.nan,
        "final_jw_l2": jw_l2[-1] if jw_l2 is not None else math.nan,
        "final_ew_l1": ew_l1[-1] if ew_l1 is not None else math.nan,
        "final_ew_l2": ew_l2[-1] if ew_l2 is not None else math.nan,
        "final_helicity_sub": helicity_sub[-1] if helicity_sub is not None else math.nan,
        "final_edotb_sub": edotb_sub[-1] if edotb_sub is not None else math.nan,
        "corr_psi0_s_leak": pearson(psi0, leak),
        "corr_psi0_s_leak_abs": pearson(psi0, leak_abs),
        "corr_psi0_jw_ew": pearson(psi0, jw_ew),
        "corr_psi0_brane_epar2": pearson(psi0, epar2),
        "corr_psi0_mixed_ew2": pearson(psi0, mixed_ew2),
        "corr_psi0_mixed_c2": pearson(psi0, mixed_c2),
        "corr_psi0_em_u_sub": pearson(psi0, em_u_sub) if em_u_sub is not None else math.nan,
        "corr_psi0_em_sw": pearson(psi0, em_sw) if em_sw is not None else math.nan,
        "corr_psi0_em_leak_w": pearson(psi0, em_leak_w)
        if em_leak_w is not None
        else math.nan,
        "corr_psi0_emf_vw_c_abs": pearson(psi0, emf_vw_c_abs)
        if emf_vw_c_abs is not None
        else math.nan,
        "corr_psi0_emf_cov_vxb_abs": pearson(psi0, emf_cov_vxb_abs)
        if emf_cov_vxb_abs is not None
        else math.nan,
        "corr_psi0_helicity_sub": pearson(psi0, helicity_sub)
        if helicity_sub is not None
        else math.nan,
        "corr_psi0_edotb_sub": pearson(psi0, edotb_sub)
        if edotb_sub is not None
        else math.nan,
        "corr_psiproj_s_leak_abs": pearson(psi_proj, leak_abs),
        "corr_psiproj_jw_ew": pearson(psi_proj, jw_ew),
        "corr_psiproj_mixed_ew2": pearson(psi_proj, mixed_ew2),
        "corr_psiw0_s_leak_abs": pearson(psi_w0, leak_abs),
        "corr_psiw0_jw_ew": pearson(psi_w0, jw_ew),
        "corr_psiw0_mixed_ew2": pearson(psi_w0, mixed_ew2),
        "corr_psiw0_emf_vw_c_abs": pearson(psi_w0, emf_vw_c_abs)
        if emf_vw_c_abs is not None
        else math.nan,
        "corr_psiw0_emf_cov_vxb_abs": pearson(psi_w0, emf_cov_vxb_abs)
        if emf_cov_vxb_abs is not None
        else math.nan,
        "final_em_a2_mode_1": em_a2_mode1[-1] if em_a2_mode1 is not None else math.nan,
        "final_em_pi2_mode_1": em_pi2_mode1[-1] if em_pi2_mode1 is not None else math.nan,
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
        "final_int_divpi0_mode0": int_divpi0_mode0[-1]
        if int_divpi0_mode0 is not None
        else math.nan,
        "final_int_divpix_mode0": int_divpix_mode0[-1]
        if int_divpix_mode0 is not None
        else math.nan,
        "final_int_divpiy_mode0": int_divpiy_mode0[-1]
        if int_divpiy_mode0 is not None
        else math.nan,
        "final_int_divpiz_mode0": int_divpiz_mode0[-1]
        if int_divpiz_mode0 is not None
        else math.nan,
        "final_int_divpiw_mode0": int_divpiw_mode0[-1]
        if int_divpiw_mode0 is not None
        else math.nan,
        "final_int_srcmomx_mode0": int_srcmomx_mode0[-1]
        if int_srcmomx_mode0 is not None
        else math.nan,
        "final_int_srcmomy_mode0": int_srcmomy_mode0[-1]
        if int_srcmomy_mode0 is not None
        else math.nan,
        "final_int_srcmomz_mode0": int_srcmomz_mode0[-1]
        if int_srcmomz_mode0 is not None
        else math.nan,
        "final_int_srcmomw_mode0": int_srcmomw_mode0[-1]
        if int_srcmomw_mode0 is not None
        else math.nan,
        "final_int_srcenergy_mode0": int_srcenergy_mode0[-1]
        if int_srcenergy_mode0 is not None
        else math.nan,
        "final_int_srcpi0_mode0": int_srcpi0_mode0[-1]
        if int_srcpi0_mode0 is not None
        else math.nan,
        "final_int_srcpix_mode0": int_srcpix_mode0[-1]
        if int_srcpix_mode0 is not None
        else math.nan,
        "final_int_srcpiy_mode0": int_srcpiy_mode0[-1]
        if int_srcpiy_mode0 is not None
        else math.nan,
        "final_int_srcpiz_mode0": int_srcpiz_mode0[-1]
        if int_srcpiz_mode0 is not None
        else math.nan,
        "final_int_srcpiw_mode0": int_srcpiw_mode0[-1]
        if int_srcpiw_mode0 is not None
        else math.nan,
        "final_int_src_em_laplacian_abs": int_src_em_laplacian_abs[-1]
        if int_src_em_laplacian_abs is not None
        else math.nan,
        "final_int_src_em_mass_abs": int_src_em_mass_abs[-1]
        if int_src_em_mass_abs is not None
        else math.nan,
        "final_int_src_em_current_abs": int_src_em_current_abs[-1]
        if int_src_em_current_abs is not None
        else math.nan,
        "final_int_src_em_damping_abs": int_src_em_damping_abs[-1]
        if int_src_em_damping_abs is not None
        else math.nan,
        "final_int_src_em_spatial_mixed_abs": int_src_em_spatial_mixed_abs[-1]
        if int_src_em_spatial_mixed_abs is not None
        else math.nan,
        "final_int_src_em_timelike_abs": int_src_em_timelike_abs[-1]
        if int_src_em_timelike_abs is not None
        else math.nan,
        "final_int_src_em_gauge_abs": int_src_em_gauge_abs[-1]
        if int_src_em_gauge_abs is not None
        else math.nan,
        "final_int_div_mode0_plasma_abs": int_div_mode0_plasma_abs[-1]
        if int_div_mode0_plasma_abs is not None
        else math.nan,
        "final_int_div_mode0_em_abs": int_div_mode0_em_abs[-1]
        if int_div_mode0_em_abs is not None
        else math.nan,
        "final_int_div_mode0_total_abs": int_div_mode0_total_abs[-1]
        if int_div_mode0_total_abs is not None
        else math.nan,
        "final_int_src_mode0_plasma_abs": int_src_mode0_plasma_abs[-1]
        if int_src_mode0_plasma_abs is not None
        else math.nan,
        "final_int_src_mode0_em_abs": int_src_mode0_em_abs[-1]
        if int_src_mode0_em_abs is not None
        else math.nan,
        "final_int_src_mode0_timelike_abs": int_src_mode0_timelike_abs[-1]
        if int_src_mode0_timelike_abs is not None
        else math.nan,
        "final_int_src_mode0_total_abs": int_src_mode0_total_abs[-1]
        if int_src_mode0_total_abs is not None
        else math.nan,
        "final_gauge_l1": gauge_l1[-1] if gauge_l1 is not None else math.nan,
        "final_gauge_l2": gauge_l2[-1] if gauge_l2 is not None else math.nan,
        "final_gauge_max_abs": gauge_max_abs[-1] if gauge_max_abs is not None else math.nan,
        "final_gauge_mode0_l1": gauge_mode0_l1[-1]
        if gauge_mode0_l1 is not None
        else math.nan,
        "final_gauge_mode0_l2": gauge_mode0_l2[-1]
        if gauge_mode0_l2 is not None
        else math.nan,
        "final_gauge_mode0_max_abs": gauge_mode0_max_abs[-1]
        if gauge_mode0_max_abs is not None
        else math.nan,
        "corr_psi0_jw_mode_l2_1": pearson(psi0, jw_mode1_l2)
        if jw_mode1_l2 is not None
        else math.nan,
        "corr_psiw0_minus_psiproj_mixed_ew2": pearson(psi_w0_minus_psi_proj, mixed_ew2),
        "corr_psiw0_minus_psiproj_emf_vw_c_abs": pearson(
            psi_w0_minus_psi_proj, emf_vw_c_abs
        )
        if emf_vw_c_abs is not None
        else math.nan,
        "corr_psiw0_minus_psiproj_emf_cov_vxb_abs": pearson(
            psi_w0_minus_psi_proj, emf_cov_vxb_abs
        )
        if emf_cov_vxb_abs is not None
        else math.nan,
        "corr_psiw0_minus_psiproj_src_mode0_total_abs": pearson(
            psi_w0_minus_psi_proj, int_src_mode0_total_abs
        )
        if int_src_mode0_total_abs is not None
        else math.nan,
        "corr_psiw0_src_mode0_total_abs": pearson(psi_w0, int_src_mode0_total_abs)
        if int_src_mode0_total_abs is not None
        else math.nan,
        "corr_psiw0_div_mode0_total_abs": pearson(psi_w0, int_div_mode0_total_abs)
        if int_div_mode0_total_abs is not None
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
        "max_dpsi0_dt": max_dpsi0_dt,
        "min_dpsi0_dt": min_dpsi0_dt,
        "max_abs_dpsi0_dt": max_abs_dpsi0_dt,
        "final_dpsi0_dt": final_dpsi0_dt,
        "time_at_max_dpsi0_dt": time_at_max_dpsi0_dt,
        "time_at_max_abs_dpsi0_dt": time_at_max_abs_dpsi0_dt,
        "max_dpsi_proj_dt": max_dpsi_proj_dt,
        "min_dpsi_proj_dt": min_dpsi_proj_dt,
        "max_abs_dpsi_proj_dt": max_abs_dpsi_proj_dt,
        "final_dpsi_proj_dt": final_dpsi_proj_dt,
        "time_at_max_dpsi_proj_dt": time_at_max_dpsi_proj_dt,
        "time_at_max_abs_dpsi_proj_dt": time_at_max_abs_dpsi_proj_dt,
        "max_dpsi_w0_dt": max_dpsi_w0_dt,
        "min_dpsi_w0_dt": min_dpsi_w0_dt,
        "max_abs_dpsi_w0_dt": max_abs_dpsi_w0_dt,
        "final_dpsi_w0_dt": final_dpsi_w0_dt,
        "time_at_max_dpsi_w0_dt": time_at_max_dpsi_w0_dt,
        "time_at_max_abs_dpsi_w0_dt": time_at_max_abs_dpsi_w0_dt,
        "time_first_abs_psi0_onset": time_first_abs_psi0_onset,
        "psi_w0_peak_min_abs": psi_w0_peak_min_abs,
        "psi_w0_peak_min_prominence": psi_w0_peak_min_prom,
        "psi_w0_peak_count_total": psi_w0_peak_count_total,
        "psi_w0_peak_count_pre_onset": psi_w0_peak_count_pre_onset,
        "psi_w0_peak_count_post_onset": psi_w0_peak_count_post_onset,
        "src_mode0_total_peak_min_abs": src_mode0_total_peak_min_abs,
        "src_mode0_total_peak_min_prominence": src_mode0_total_peak_min_prom,
        "src_mode0_total_peak_count_total": src_mode0_total_peak_count_total,
        "src_mode0_total_peak_count_pre_onset": src_mode0_total_peak_count_pre_onset,
        "src_mode0_total_peak_count_post_onset": src_mode0_total_peak_count_post_onset,
        "em_bulk_ledger_max_abs_rate": em_bulk_ledger_max_abs_rate,
        "em_bulk_ledger_rms_abs_rate": em_bulk_ledger_rms_abs_rate,
        "em_bulk_ledger_final_rate": em_bulk_ledger_final_rate,
        "em_bulk_ledger_status": em_bulk_ledger_status,
        "momx_transport_max_norm": momx_transport_max_norm,
        "momx_transport_rms_norm": momx_transport_rms_norm,
        "momx_transport_max_abs_rate": momx_transport_max_abs_rate,
        "momx_transport_final_rate": momx_transport_final_rate,
        "momy_transport_max_norm": momy_transport_max_norm,
        "momy_transport_rms_norm": momy_transport_rms_norm,
        "momy_transport_max_abs_rate": momy_transport_max_abs_rate,
        "momy_transport_final_rate": momy_transport_final_rate,
        "momz_transport_max_norm": momz_transport_max_norm,
        "momz_transport_rms_norm": momz_transport_rms_norm,
        "momz_transport_max_abs_rate": momz_transport_max_abs_rate,
        "momz_transport_final_rate": momz_transport_final_rate,
        "momw_transport_max_norm": momw_transport_max_norm,
        "momw_transport_rms_norm": momw_transport_rms_norm,
        "momw_transport_max_abs_rate": momw_transport_max_abs_rate,
        "momw_transport_final_rate": momw_transport_final_rate,
        "energy_transport_max_norm": energy_transport_max_norm,
        "energy_transport_rms_norm": energy_transport_rms_norm,
        "energy_transport_max_abs_rate": energy_transport_max_abs_rate,
        "energy_transport_final_rate": energy_transport_final_rate,
        "pi0_transport_max_norm": pi0_transport_max_norm,
        "pi0_transport_rms_norm": pi0_transport_rms_norm,
        "pi0_transport_max_abs_rate": pi0_transport_max_abs_rate,
        "pi0_transport_final_rate": pi0_transport_final_rate,
        "pix_transport_max_norm": pix_transport_max_norm,
        "pix_transport_rms_norm": pix_transport_rms_norm,
        "pix_transport_max_abs_rate": pix_transport_max_abs_rate,
        "pix_transport_final_rate": pix_transport_final_rate,
        "piy_transport_max_norm": piy_transport_max_norm,
        "piy_transport_rms_norm": piy_transport_rms_norm,
        "piy_transport_max_abs_rate": piy_transport_max_abs_rate,
        "piy_transport_final_rate": piy_transport_final_rate,
        "piz_transport_max_norm": piz_transport_max_norm,
        "piz_transport_rms_norm": piz_transport_rms_norm,
        "piz_transport_max_abs_rate": piz_transport_max_abs_rate,
        "piz_transport_final_rate": piz_transport_final_rate,
        "piw_transport_max_norm": piw_transport_max_norm,
        "piw_transport_rms_norm": piw_transport_rms_norm,
        "piw_transport_max_abs_rate": piw_transport_max_abs_rate,
        "piw_transport_final_rate": piw_transport_final_rate,
        "momx_transport_status": momx_transport_status,
        "momy_transport_status": momy_transport_status,
        "momz_transport_status": momz_transport_status,
        "momw_transport_status": momw_transport_status,
        "energy_transport_status": energy_transport_status,
        "pi0_transport_status": pi0_transport_status,
        "pix_transport_status": pix_transport_status,
        "piy_transport_status": piy_transport_status,
        "piz_transport_status": piz_transport_status,
        "piw_transport_status": piw_transport_status,
        "em_transport_closure_status": em_transport_closure_status,
        "transport_closure_status": transport_closure_status,
    }
    return result


def run_case(binary, input_path, workdir, output_hst, extra_args=None):
    for filename in ("parthenon.out0.hst", "parthenon.out1.hst"):
        try:
            (workdir / filename).unlink()
        except FileNotFoundError:
            pass

    cmd = [str(binary), "-i", str(input_path)]
    if extra_args:
        cmd.extend(extra_args)
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
        "final_psi_proj_span",
        "final_psi_proj_matched_span",
        "final_psi_proj_point_span",
        "final_psi_proj_gaussian_span",
        "final_psi_w0_span",
        "final_psi_w0_minus_psi_proj_span",
        "final_psi_w0_minus_psi0_span",
        "final_psi_w0_over_psi_proj_span",
        "final_psi_proj_over_psi0_span",
        "final_psi_w0_over_psi0_span",
        "max_abs_psi_w0_minus_psi_proj_span",
        "final_s_leak",
        "final_s_leak_abs",
        "final_jw_ew",
        "final_ja_ea",
        "final_brane_epar2",
        "final_mixed_ew2",
        "final_mixed_c2",
        "final_em_u_bulk",
        "final_em_u_resolved",
        "final_em_u_sub",
        "final_em_sw",
        "final_em_leak_w",
        "final_emf_vw_c_abs",
        "final_emf_vw_c_par_abs",
        "final_emf_cov_vxb_abs",
        "final_emf_cov_vxb_par_abs",
        "final_jw_l1",
        "final_jw_l2",
        "final_ew_l1",
        "final_ew_l2",
        "final_helicity_sub",
        "final_edotb_sub",
        "final_em_a2_mode_1",
        "final_em_pi2_mode_1",
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
        "final_int_divpi0_mode0",
        "final_int_divpix_mode0",
        "final_int_divpiy_mode0",
        "final_int_divpiz_mode0",
        "final_int_divpiw_mode0",
        "final_int_srcmomx_mode0",
        "final_int_srcmomy_mode0",
        "final_int_srcmomz_mode0",
        "final_int_srcmomw_mode0",
        "final_int_srcenergy_mode0",
        "final_int_srcpi0_mode0",
        "final_int_srcpix_mode0",
        "final_int_srcpiy_mode0",
        "final_int_srcpiz_mode0",
        "final_int_srcpiw_mode0",
        "final_int_src_em_laplacian_abs",
        "final_int_src_em_mass_abs",
        "final_int_src_em_current_abs",
        "final_int_src_em_damping_abs",
        "final_int_src_em_spatial_mixed_abs",
        "final_int_src_em_timelike_abs",
        "final_int_src_em_gauge_abs",
        "final_int_div_mode0_plasma_abs",
        "final_int_div_mode0_em_abs",
        "final_int_div_mode0_total_abs",
        "final_int_src_mode0_plasma_abs",
        "final_int_src_mode0_em_abs",
        "final_int_src_mode0_timelike_abs",
        "final_int_src_mode0_total_abs",
        "final_gauge_l1",
        "final_gauge_l2",
        "final_gauge_max_abs",
        "final_gauge_mode0_l1",
        "final_gauge_mode0_l2",
        "final_gauge_mode0_max_abs",
        "corr_psi0_s_leak",
        "corr_psi0_s_leak_abs",
        "corr_psi0_jw_ew",
        "corr_psi0_brane_epar2",
        "corr_psi0_mixed_ew2",
        "corr_psi0_mixed_c2",
        "corr_psi0_em_u_sub",
        "corr_psi0_em_sw",
        "corr_psi0_em_leak_w",
        "corr_psi0_emf_vw_c_abs",
        "corr_psi0_emf_cov_vxb_abs",
        "corr_psi0_helicity_sub",
        "corr_psi0_edotb_sub",
        "corr_psi0_jw_mode_l2_1",
        "corr_psiproj_s_leak_abs",
        "corr_psiproj_jw_ew",
        "corr_psiproj_mixed_ew2",
        "corr_psiw0_s_leak_abs",
        "corr_psiw0_jw_ew",
        "corr_psiw0_mixed_ew2",
        "corr_psiw0_emf_vw_c_abs",
        "corr_psiw0_emf_cov_vxb_abs",
        "corr_psiw0_minus_psiproj_mixed_ew2",
        "corr_psiw0_minus_psiproj_emf_vw_c_abs",
        "corr_psiw0_minus_psiproj_emf_cov_vxb_abs",
        "corr_psiw0_minus_psiproj_src_mode0_total_abs",
        "corr_psiw0_src_mode0_total_abs",
        "corr_psiw0_div_mode0_total_abs",
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
        "max_dpsi0_dt",
        "min_dpsi0_dt",
        "max_abs_dpsi0_dt",
        "final_dpsi0_dt",
        "time_at_max_dpsi0_dt",
        "time_at_max_abs_dpsi0_dt",
        "max_dpsi_proj_dt",
        "min_dpsi_proj_dt",
        "max_abs_dpsi_proj_dt",
        "final_dpsi_proj_dt",
        "time_at_max_dpsi_proj_dt",
        "time_at_max_abs_dpsi_proj_dt",
        "max_dpsi_w0_dt",
        "min_dpsi_w0_dt",
        "max_abs_dpsi_w0_dt",
        "final_dpsi_w0_dt",
        "time_at_max_dpsi_w0_dt",
        "time_at_max_abs_dpsi_w0_dt",
        "time_first_abs_psi0_onset",
        "psi_w0_peak_min_abs",
        "psi_w0_peak_min_prominence",
        "psi_w0_peak_count_total",
        "psi_w0_peak_count_pre_onset",
        "psi_w0_peak_count_post_onset",
        "src_mode0_total_peak_min_abs",
        "src_mode0_total_peak_min_prominence",
        "src_mode0_total_peak_count_total",
        "src_mode0_total_peak_count_pre_onset",
        "src_mode0_total_peak_count_post_onset",
        "em_bulk_ledger_max_abs_rate",
        "em_bulk_ledger_rms_abs_rate",
        "em_bulk_ledger_final_rate",
        "em_bulk_ledger_status",
        "momx_transport_max_norm",
        "momx_transport_rms_norm",
        "momx_transport_max_abs_rate",
        "momx_transport_final_rate",
        "momx_transport_status",
        "momy_transport_max_norm",
        "momy_transport_rms_norm",
        "momy_transport_max_abs_rate",
        "momy_transport_final_rate",
        "momy_transport_status",
        "momz_transport_max_norm",
        "momz_transport_rms_norm",
        "momz_transport_max_abs_rate",
        "momz_transport_final_rate",
        "momz_transport_status",
        "momw_transport_max_norm",
        "momw_transport_rms_norm",
        "momw_transport_max_abs_rate",
        "momw_transport_final_rate",
        "momw_transport_status",
        "energy_transport_max_norm",
        "energy_transport_rms_norm",
        "energy_transport_max_abs_rate",
        "energy_transport_final_rate",
        "energy_transport_status",
        "pi0_transport_max_norm",
        "pi0_transport_rms_norm",
        "pi0_transport_max_abs_rate",
        "pi0_transport_final_rate",
        "pi0_transport_status",
        "pix_transport_max_norm",
        "pix_transport_rms_norm",
        "pix_transport_max_abs_rate",
        "pix_transport_final_rate",
        "pix_transport_status",
        "piy_transport_max_norm",
        "piy_transport_rms_norm",
        "piy_transport_max_abs_rate",
        "piy_transport_final_rate",
        "piy_transport_status",
        "piz_transport_max_norm",
        "piz_transport_rms_norm",
        "piz_transport_max_abs_rate",
        "piz_transport_final_rate",
        "piz_transport_status",
        "piw_transport_max_norm",
        "piw_transport_rms_norm",
        "piw_transport_max_abs_rate",
        "piw_transport_final_rate",
        "piw_transport_status",
        "em_transport_closure_status",
        "transport_closure_status",
        "correlation_status",
        "correlation_failures",
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
    min_em_a2_mode_1,
    min_em_pi2_mode_1,
    min_jw_mode_l2_1,
    min_jw_mode_activity_proxy,
    min_em_leak_w_abs,
    min_helicity_sub_abs,
    min_edotb_sub_abs,
    min_src_em_laplacian_abs,
    min_src_em_mass_abs,
    min_src_em_current_abs,
    min_src_em_spatial_mixed_abs,
    min_src_em_timelike_abs,
    min_emf_vw_c_abs,
    min_emf_cov_vxb_abs,
    min_jw_l1,
    min_ew_l1,
):
    failures = []
    if min_jw_ew_abs > 0.0 and abs(row["final_jw_ew"]) < min_jw_ew_abs:
        failures.append("jw_ew")
    if min_s_leak_abs > 0.0 and row["final_s_leak_abs"] < min_s_leak_abs:
        failures.append("s_leak_abs")
    if min_mixed_ew2 > 0.0 and row["final_mixed_ew2"] < min_mixed_ew2:
        failures.append("mixed_ew2")
    if min_em_a2_mode_1 > 0.0 and row["final_em_a2_mode_1"] < min_em_a2_mode_1:
        failures.append("em_a2_mode_1")
    if min_em_pi2_mode_1 > 0.0 and row["final_em_pi2_mode_1"] < min_em_pi2_mode_1:
        failures.append("em_pi2_mode_1")
    if min_jw_mode_l2_1 > 0.0 and row["final_jw_mode_l2_1"] < min_jw_mode_l2_1:
        failures.append("jw_mode_l2_1")
    jw_mode_activity_proxy = mode_activity_proxy(
        row["final_jw_mode_l2_1"], row["final_jw_ew"]
    )
    if (
        min_jw_mode_activity_proxy > 0.0
        and jw_mode_activity_proxy < min_jw_mode_activity_proxy
    ):
        failures.append("jw_mode_activity_proxy")
    if min_em_leak_w_abs > 0.0 and abs(row["final_em_leak_w"]) < min_em_leak_w_abs:
        failures.append("em_leak_w")
    if min_helicity_sub_abs > 0.0 and abs(row["final_helicity_sub"]) < min_helicity_sub_abs:
        failures.append("helicity_sub")
    if min_edotb_sub_abs > 0.0 and abs(row["final_edotb_sub"]) < min_edotb_sub_abs:
        failures.append("edotb_sub")
    if (
        min_src_em_laplacian_abs > 0.0
        and abs(row["final_int_src_em_laplacian_abs"]) < min_src_em_laplacian_abs
    ):
        failures.append("src_em_laplacian_abs")
    if min_src_em_mass_abs > 0.0 and abs(row["final_int_src_em_mass_abs"]) < min_src_em_mass_abs:
        failures.append("src_em_mass_abs")
    if (
        min_src_em_current_abs > 0.0
        and abs(row["final_int_src_em_current_abs"]) < min_src_em_current_abs
    ):
        failures.append("src_em_current_abs")
    if (
        min_src_em_spatial_mixed_abs > 0.0
        and abs(row["final_int_src_em_spatial_mixed_abs"]) < min_src_em_spatial_mixed_abs
    ):
        failures.append("src_em_spatial_mixed_abs")
    if (
        min_src_em_timelike_abs > 0.0
        and abs(row["final_int_src_em_timelike_abs"]) < min_src_em_timelike_abs
    ):
        failures.append("src_em_timelike_abs")
    if min_emf_vw_c_abs > 0.0 and abs(row["final_emf_vw_c_abs"]) < min_emf_vw_c_abs:
        failures.append("emf_vw_c_abs")
    if (
        min_emf_cov_vxb_abs > 0.0
        and abs(row["final_emf_cov_vxb_abs"]) < min_emf_cov_vxb_abs
    ):
        failures.append("emf_cov_vxb_abs")
    if min_jw_l1 > 0.0 and abs(row["final_jw_l1"]) < min_jw_l1:
        failures.append("jw_l1")
    if min_ew_l1 > 0.0 and abs(row["final_ew_l1"]) < min_ew_l1:
        failures.append("ew_l1")

    if failures:
        return ("FAIL", "|".join(failures))
    return ("PASS", "none")


def evaluate_full_correlation(
    row,
    min_abs_corr_psi0_s_leak_abs,
    min_abs_corr_psi0_jw_ew,
    min_abs_corr_psi0_mixed_ew2,
    min_abs_corr_psi0_mixed_c2,
    min_abs_corr_psi0_jw_mode_l2_1,
    min_abs_corr_psi0_em_leak_w,
    min_abs_corr_psi0_emf_vw_c_abs,
    min_abs_corr_psi0_emf_cov_vxb_abs,
    min_abs_corr_psi0_helicity_sub,
    min_abs_corr_psi0_edotb_sub,
):
    failures = []

    def check(key, threshold, label):
        if threshold <= 0.0:
            return
        value = row.get(key, math.nan)
        if math.isnan(value) or abs(value) < threshold:
            failures.append(label)

    check("corr_psi0_s_leak_abs", min_abs_corr_psi0_s_leak_abs, "corr_psi0_s_leak_abs")
    check("corr_psi0_jw_ew", min_abs_corr_psi0_jw_ew, "corr_psi0_jw_ew")
    check("corr_psi0_mixed_ew2", min_abs_corr_psi0_mixed_ew2, "corr_psi0_mixed_ew2")
    check("corr_psi0_mixed_c2", min_abs_corr_psi0_mixed_c2, "corr_psi0_mixed_c2")
    check(
        "corr_psi0_jw_mode_l2_1",
        min_abs_corr_psi0_jw_mode_l2_1,
        "corr_psi0_jw_mode_l2_1",
    )
    check("corr_psi0_em_leak_w", min_abs_corr_psi0_em_leak_w, "corr_psi0_em_leak_w")
    check(
        "corr_psi0_emf_vw_c_abs",
        min_abs_corr_psi0_emf_vw_c_abs,
        "corr_psi0_emf_vw_c_abs",
    )
    check(
        "corr_psi0_emf_cov_vxb_abs",
        min_abs_corr_psi0_emf_cov_vxb_abs,
        "corr_psi0_emf_cov_vxb_abs",
    )
    check(
        "corr_psi0_helicity_sub",
        min_abs_corr_psi0_helicity_sub,
        "corr_psi0_helicity_sub",
    )
    check("corr_psi0_edotb_sub", min_abs_corr_psi0_edotb_sub, "corr_psi0_edotb_sub")

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
        "--athena-arg",
        action="append",
        default=[],
        help="Additional athenaPK runtime override passed to both controlled/full runs",
    )
    parser.add_argument(
        "--controlled-arg",
        action="append",
        default=[],
        help="Additional athenaPK runtime override passed only to controlled run",
    )
    parser.add_argument(
        "--full-arg",
        action="append",
        default=[],
        help="Additional athenaPK runtime override passed only to full run",
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
        "--check-transport-closure",
        action="store_true",
        help="Include plasma+EM transport closure status in fail-on-check logic",
    )
    parser.add_argument(
        "--check-em-bulk-ledger",
        action="store_true",
        help="Include EM bulk-ledger closure status in fail-on-check logic",
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
        "--full-min-em-a2-mode-1",
        type=float,
        default=0.0,
        help="Minimum final_em_a2_mode_1 required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-em-pi2-mode-1",
        type=float,
        default=0.0,
        help="Minimum final_em_pi2_mode_1 required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-jw-mode-l2-1",
        type=float,
        default=0.0,
        help="Minimum final_jw_mode_l2_1 required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-jw-mode-activity-proxy",
        type=float,
        default=0.0,
        help="Minimum max(final_jw_mode_l2_1, |final_jw_ew|) required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-em-leak-w",
        type=float,
        default=0.0,
        help="Minimum |final_em_leak_w| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-helicity-sub",
        type=float,
        default=0.0,
        help="Minimum |final_helicity_sub| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-edotb-sub",
        type=float,
        default=0.0,
        help="Minimum |final_edotb_sub| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-src-em-laplacian-abs",
        type=float,
        default=0.0,
        help="Minimum |final_int_src_em_laplacian_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-src-em-mass-abs",
        type=float,
        default=0.0,
        help="Minimum |final_int_src_em_mass_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-src-em-current-abs",
        type=float,
        default=0.0,
        help="Minimum |final_int_src_em_current_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-src-em-spatial-mixed-abs",
        type=float,
        default=0.0,
        help="Minimum |final_int_src_em_spatial_mixed_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-src-em-timelike-abs",
        type=float,
        default=0.0,
        help="Minimum |final_int_src_em_timelike_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-emf-vw-c-abs",
        type=float,
        default=0.0,
        help="Minimum |final_emf_vw_c_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-emf-cov-vxb-abs",
        type=float,
        default=0.0,
        help="Minimum |final_emf_cov_vxb_abs| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-jw-l1",
        type=float,
        default=0.0,
        help="Minimum |final_jw_l1| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-ew-l1",
        type=float,
        default=0.0,
        help="Minimum |final_ew_l1| required for full-case activity PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-s-leak-abs",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, int_s_leak_abs)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-jw-ew",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, int_jw_ew)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-mixed-ew2",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, mixed_ew2)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-mixed-c2",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, mixed_c2)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-jw-mode-l2-1",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, jw_mode_l2_1)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-em-leak-w",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, em_leak_w)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-emf-vw-c-abs",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, emf_vw_c_abs)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-emf-cov-vxb-abs",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, emf_cov_vxb_abs)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-helicity-sub",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, helicity_sub)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--full-min-abs-corr-psi0-edotb-sub",
        type=float,
        default=0.0,
        help="Minimum |corr(psi0_span, edotb_sub)| for full-case correlation PASS (disabled when 0)",
    )
    parser.add_argument(
        "--em-bulk-ledger-abs-rate-tol",
        type=float,
        default=1.0,
        help="Maximum absolute EM bulk-ledger residual rate for PASS",
    )
    parser.add_argument(
        "--cycle-onset-psi0-abs-threshold",
        type=float,
        default=1.0,
        help="Absolute psi0 threshold defining onset time for pre/post peak diagnostics",
    )
    parser.add_argument(
        "--cycle-peak-min-amplitude-rel",
        type=float,
        default=1.0e-6,
        help="Minimum local-peak absolute amplitude as a fraction of series max|value|",
    )
    parser.add_argument(
        "--cycle-peak-min-amplitude-abs",
        type=float,
        default=1.0e-12,
        help="Minimum local-peak absolute amplitude floor",
    )
    parser.add_argument(
        "--cycle-peak-min-prominence-rel",
        type=float,
        default=1.0e-3,
        help="Minimum local-peak prominence as a fraction of series max|value|",
    )
    parser.add_argument(
        "--cycle-peak-min-prominence-abs",
        type=float,
        default=1.0e-12,
        help="Minimum local-peak prominence floor",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    controlled_input = resolve_input_path(args.controlled_input, repo_root, workdir)
    full_input = resolve_input_path(args.full_input, repo_root, workdir)
    output_dir_arg = Path(args.output_dir)
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (workdir / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    controlled_args = list(args.athena_arg) + list(args.controlled_arg)
    full_args = list(args.athena_arg) + list(args.full_arg)
    controlled_cols = run_case(
        binary,
        controlled_input,
        workdir,
        output_dir / "harris_controlled.out1.hst",
        controlled_args,
    )
    full_cols = run_case(
        binary,
        full_input,
        workdir,
        output_dir / "harris_full.out1.hst",
        full_args,
    )

    results = [
        analyze_case(
            "controlled",
            controlled_cols,
            args.closure_norm_tol,
            args.closure_abs_rate_tol,
            args.closure_local_mode0_abs_rate_tol,
            args.em_bulk_ledger_abs_rate_tol,
            args.cycle_onset_psi0_abs_threshold,
            args.cycle_peak_min_amplitude_rel,
            args.cycle_peak_min_amplitude_abs,
            args.cycle_peak_min_prominence_rel,
            args.cycle_peak_min_prominence_abs,
        ),
        analyze_case(
            "full",
            full_cols,
            args.closure_norm_tol,
            args.closure_abs_rate_tol,
            args.closure_local_mode0_abs_rate_tol,
            args.em_bulk_ledger_abs_rate_tol,
            args.cycle_onset_psi0_abs_threshold,
            args.cycle_peak_min_amplitude_rel,
            args.cycle_peak_min_amplitude_abs,
            args.cycle_peak_min_prominence_rel,
            args.cycle_peak_min_prominence_abs,
        ),
    ]
    for row in results:
        row["correlation_status"] = "N/A"
        row["correlation_failures"] = "n/a"
        row["activity_status"] = "N/A"
        row["activity_failures"] = "n/a"
        if row["case"] == "full":
            corr_status, corr_failures = evaluate_full_correlation(
                row,
                args.full_min_abs_corr_psi0_s_leak_abs,
                args.full_min_abs_corr_psi0_jw_ew,
                args.full_min_abs_corr_psi0_mixed_ew2,
                args.full_min_abs_corr_psi0_mixed_c2,
                args.full_min_abs_corr_psi0_jw_mode_l2_1,
                args.full_min_abs_corr_psi0_em_leak_w,
                args.full_min_abs_corr_psi0_emf_vw_c_abs,
                args.full_min_abs_corr_psi0_emf_cov_vxb_abs,
                args.full_min_abs_corr_psi0_helicity_sub,
                args.full_min_abs_corr_psi0_edotb_sub,
            )
            row["correlation_status"] = corr_status
            row["correlation_failures"] = corr_failures
            status, failures = evaluate_full_activity(
                row,
                args.full_min_jw_ew_abs,
                args.full_min_s_leak_abs,
                args.full_min_mixed_ew2,
                args.full_min_em_a2_mode_1,
                args.full_min_em_pi2_mode_1,
                args.full_min_jw_mode_l2_1,
                args.full_min_jw_mode_activity_proxy,
                args.full_min_abs_em_leak_w,
                args.full_min_abs_helicity_sub,
                args.full_min_abs_edotb_sub,
                args.full_min_src_em_laplacian_abs,
                args.full_min_src_em_mass_abs,
                args.full_min_src_em_current_abs,
                args.full_min_src_em_spatial_mixed_abs,
                args.full_min_src_em_timelike_abs,
                args.full_min_emf_vw_c_abs,
                args.full_min_emf_cov_vxb_abs,
                args.full_min_jw_l1,
                args.full_min_ew_l1,
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
            or r["correlation_status"] == "FAIL"
            or r["activity_status"] == "FAIL"
            or (
                args.check_transport_closure
                and r["transport_closure_status"] == "FAIL"
            )
            or (
                args.check_em_bulk_ledger
                and r["em_bulk_ledger_status"] == "FAIL"
            )
        ]
        if failing:
            raise SystemExit(f"scan checks failed for: {', '.join(failing)}")


if __name__ == "__main__":
    main()
