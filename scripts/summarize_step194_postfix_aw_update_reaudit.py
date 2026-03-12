#!/usr/bin/env python3
"""Post-fix A_w update-power re-audit on the corrected baseline."""

from __future__ import annotations

import argparse, csv, importlib.util, math
from pathlib import Path


def load_hsm(script_dir: Path):
    module_path = script_dir / 'harris_scan_matrix.py'
    spec = importlib.util.spec_from_file_location('harris_scan_matrix', module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Could not load module spec from {module_path}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def derivative_series(times, values):
    if values is None or len(times) < 2 or len(values) != len(times):
        return [], []
    out_t, out_v = [], []
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0.0:
            continue
        out_t.append(times[i])
        out_v.append((values[i] - values[i - 1]) / dt)
    return out_t, out_v


def pearson(xs, ys):
    if len(xs) != len(ys) or len(xs) < 2:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def fit_single(y, x):
    xx = sum(v * v for v in x)
    if xx <= 0.0:
        return math.nan
    return sum(a * b for a, b in zip(y, x)) / xx


def finite_stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return (math.nan, math.nan, math.nan)
    return (
        max(abs(v) for v in finite),
        math.sqrt(sum(v * v for v in finite) / len(finite)),
        finite[-1],
    )


def window(times, values):
    start, end = times[0], times[-1]
    lo = start + (end - start) * 2.0 / 3.0
    return [v for t, v in zip(times, values) if t >= lo]


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return 'nan'
    if v == math.inf:
        return 'inf'
    return f'{float(v):.6e}'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-root', required=True)
    ap.add_argument('--out-csv', required=True)
    ap.add_argument('--out-md', required=True)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    hsm = load_hsm(script_dir)
    hst = Path(args.input_root) / 'baseline_quadrature' / 'harris_full.out1.hst'
    if not hst.exists():
        raise SystemExit(f'Missing HST: {hst}')
    cols = hsm.parse_hst(hst)

    times = cols['time']
    t_r, du = derivative_series(times, cols['m4d_em_u_bulk'])
    _, dja = derivative_series(times, cols['m4d_int_ja_ea'])
    _, djw = derivative_series(times, cols['m4d_int_jw_ew'])
    residual = [u + a + w for u, a, w in zip(du, dja, djw)]
    y = window(t_r, residual)
    _, residual_rms, _ = finite_stats(y)

    def ds(name):
        return derivative_series(times, cols[name])[1]

    def ws(name):
        return window(t_r, ds(name))

    x_gradx = [sum(vals) for vals in zip(ws('m4d_int_aw_mode0_gradx_power_da_dt_discrete'), ws('m4d_int_aw_mode1_gradx_power_da_dt_discrete'), ws('m4d_int_aw_mode2_gradx_power_da_dt_discrete'), ws('m4d_int_aw_mode3_gradx_power_da_dt_discrete'))]
    x_gradz = [sum(vals) for vals in zip(ws('m4d_int_aw_mode0_gradz_power_da_dt_discrete'), ws('m4d_int_aw_mode1_gradz_power_da_dt_discrete'), ws('m4d_int_aw_mode2_gradz_power_da_dt_discrete'), ws('m4d_int_aw_mode3_gradz_power_da_dt_discrete'))]
    x_gradxz = [a + b for a, b in zip(x_gradx, x_gradz)]
    x_pi_drive_x = ws('m4d_int_aw_gradx_power_pi_drive_discrete_sum')
    x_pi_drive_z = ws('m4d_int_aw_gradz_power_pi_drive_discrete_sum')
    x_pi_drive_xz = [a + b for a, b in zip(x_pi_drive_x, x_pi_drive_z)]
    mode1_target = ws('m4d_mixed_cz_daw_dz2_mode1')
    mode2_target = ws('m4d_mixed_cz_daw_dz2_mode2')
    mode1_da_dt = ws('m4d_int_aw_mode1_gradz_power_da_dt_discrete')
    mode2_da_dt = ws('m4d_int_aw_mode2_gradz_power_da_dt_discrete')
    mode1_pi = ws('m4d_int_aw_mode1_gradz_power_pi_drive_discrete')
    mode2_pi = ws('m4d_int_aw_mode2_gradz_power_pi_drive_discrete')

    rows = []
    def add(component, vals, target=y):
        _, rms_abs, final = finite_stats(vals)
        alpha = fit_single(target, vals)
        corrected = [yy - alpha * vv for yy, vv in zip(target, vals)] if math.isfinite(alpha) else []
        _, corrected_rms, _ = finite_stats(corrected)
        rows.append({
            'component': component,
            'target_rms_abs_rate': finite_stats(target)[1],
            'rms_abs_rate': rms_abs,
            'corr_with_target': pearson(vals, target),
            'fit_coeff': alpha,
            'corrected_rms_abs_rate': corrected_rms,
            'reduction_factor': (finite_stats(target)[1] / corrected_rms) if corrected_rms > 0 else math.inf,
            'final_rate': final,
        })

    add('bulk_residual', y, y)
    add('aw_da_dt_gradx_sum', x_gradx, y)
    add('aw_da_dt_gradz_sum', x_gradz, y)
    add('aw_da_dt_gradxz_sum', x_gradxz, y)
    add('aw_pi_drive_gradx_sum', x_pi_drive_x, y)
    add('aw_pi_drive_gradz_sum', x_pi_drive_z, y)
    add('aw_pi_drive_gradxz_sum', x_pi_drive_xz, y)
    add('mode1_target_vs_da_dt', mode1_da_dt, mode1_target)
    add('mode1_target_vs_pi_drive', mode1_pi, mode1_target)
    add('mode2_target_vs_da_dt', mode2_da_dt, mode2_target)
    add('mode2_target_vs_pi_drive', mode2_pi, mode2_target)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    out_md = Path(args.out_md)
    with out_md.open('w', encoding='utf-8') as f:
        f.write('# Step-194 Post-fix A_w Update Re-audit\n\n')
        f.write(f'- residual_rms_abs_rate: `{fmt(residual_rms)}`\n\n')
        f.write('| component | target_rms_abs_rate | corr_with_target | fit_coeff | corrected_rms_abs_rate | reduction_factor |\n')
        f.write('|:---|---:|---:|---:|---:|---:|\n')
        for row in rows:
            f.write(
                f"| {row['component']} | {fmt(row['target_rms_abs_rate'])} | {fmt(row['corr_with_target'])} | {fmt(row['fit_coeff'])} | {fmt(row['corrected_rms_abs_rate'])} | {fmt(row['reduction_factor'])} |\n"
            )

if __name__ == '__main__':
    main()
