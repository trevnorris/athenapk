#!/usr/bin/env python3
"""Audit projection-kernel mode coefficients against the Modes4D basis."""

import argparse
import csv
import math
import re
from pathlib import Path


PI = 3.141592653589793238462643383279502884


def hermite_physicists(n, x):
    if n < 0:
        raise ValueError("n must be non-negative")
    if n == 0:
        return (1.0, 0.0)
    hm2 = 1.0
    hm1 = 2.0 * x
    if n == 1:
        return (hm1, hm2)
    for k in range(2, n + 1):
        h = (2.0 * x * hm1) - (2.0 * float(k - 1) * hm2)
        hm2 = hm1
        hm1 = h
    return (hm1, hm2)


def compute_gauss_hermite_standard(n):
    if n <= 0:
        raise ValueError("Gauss-Hermite order must be positive")

    # Match src/4d_modes/mode_tables.cpp Jacobi-matrix eigen solve.
    a = [[0.0 for _ in range(n)] for _ in range(n)]
    v = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        v[i][i] = 1.0
    for i in range(n - 1):
        beta = math.sqrt(float(i + 1) / 2.0)
        a[i][i + 1] = beta
        a[i + 1][i] = beta

    tol = 1.0e-15
    max_iters = 128 * n * n
    for _ in range(max_iters):
        p = 0
        q = 1
        max_offdiag = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                val = abs(a[i][j])
                if val > max_offdiag:
                    max_offdiag = val
                    p = i
                    q = j
        if max_offdiag < tol:
            break

        app = a[p][p]
        aqq = a[q][q]
        apq = a[p][q]
        tau = (aqq - app) / (2.0 * apq)
        t = (1.0 if tau >= 0.0 else -1.0) / (abs(tau) + math.sqrt(1.0 + tau * tau))
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = t * c

        a[p][p] = app - (t * apq)
        a[q][q] = aqq + (t * apq)
        a[p][q] = 0.0
        a[q][p] = 0.0

        for k in range(n):
            if k == p or k == q:
                continue
            akp = a[k][p]
            akq = a[k][q]
            a[k][p] = (c * akp) - (s * akq)
            a[p][k] = a[k][p]
            a[k][q] = (s * akp) + (c * akq)
            a[q][k] = a[k][q]

        for k in range(n):
            vkp = v[k][p]
            vkq = v[k][q]
            v[k][p] = (c * vkp) - (s * vkq)
            v[k][q] = (s * vkp) + (c * vkq)

    node_weight_pairs = []
    for i in range(n):
        node = a[i][i]
        weight = math.sqrt(PI) * v[0][i] * v[0][i]
        node_weight_pairs.append((node, weight))
    node_weight_pairs.sort(key=lambda pair: pair[0])

    nodes = [p[0] for p in node_weight_pairs]
    weights = [p[1] for p in node_weight_pairs]
    return (nodes, weights)


def build_mode_tables(n_modes, n_quadrature, lambda_value):
    std_nodes, std_weights = compute_gauss_hermite_standard(n_quadrature)
    nodes = [lambda_value * x for x in std_nodes]
    weights = [lambda_value * w for w in std_weights]

    phi = [[0.0 for _ in range(n_quadrature)] for _ in range(n_modes)]
    normalizer_factor = lambda_value * math.sqrt(PI)
    for n in range(n_modes):
        two_to_n = math.ldexp(1.0, n)
        n_factorial = math.gamma(float(n) + 1.0)
        normalization = math.sqrt(normalizer_factor * two_to_n * n_factorial)
        for q in range(n_quadrature):
            x = nodes[q] / lambda_value
            hn, _ = hermite_physicists(n, x)
            phi[n][q] = hn / normalization
    return (nodes, weights, phi)


def brane_point_coefficients(n_modes, lambda_value):
    hermite_zero = [0.0 for _ in range(n_modes)]
    if n_modes > 0:
        hermite_zero[0] = 1.0
    if n_modes > 1:
        hermite_zero[1] = 0.0
    for n in range(1, n_modes - 1):
        hermite_zero[n + 1] = -2.0 * float(n) * hermite_zero[n - 1]

    normalizer_factor = lambda_value * math.sqrt(PI)
    coeff = [0.0 for _ in range(n_modes)]
    for n in range(n_modes):
        two_to_n = math.ldexp(1.0, n)
        n_factorial = math.gamma(float(n) + 1.0)
        normalization = math.sqrt(normalizer_factor * two_to_n * n_factorial)
        coeff[n] = hermite_zero[n] / normalization
    return coeff


def projection_coefficients(kernel, sigma_factor, n_modes, lambda_value, nodes, weights, phi):
    if kernel == "point":
        return brane_point_coefficients(n_modes, lambda_value)

    z_int = lambda_value * math.sqrt(PI)
    coeff = [0.0 for _ in range(n_modes)]
    n_quad = len(nodes)

    if kernel == "gaussian":
        lambda_sq = lambda_value * lambda_value
        sigma_sq = sigma_factor * sigma_factor
        w_norm = sigma_factor * z_int
        for n in range(n_modes):
            mode_coeff = 0.0
            for q in range(n_quad):
                w = nodes[q]
                z = math.exp(-(w * w) / lambda_sq)
                if z <= 0.0:
                    continue
                w_kernel = math.exp(-(w * w) / (sigma_sq * lambda_sq)) / w_norm
                mode_coeff += weights[q] * (w_kernel / z) * phi[n][q]
            coeff[n] = mode_coeff
        return coeff

    if kernel == "matched":
        for n in range(n_modes):
            coeff[n] = sum(weights[q] * phi[n][q] for q in range(n_quad)) / z_int
        return coeff

    raise ValueError(f"Unsupported kernel: {kernel}")


def parse_hst(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = None
    header_idx = -1
    for idx, line in enumerate(lines):
        if line.startswith("# [1]="):
            header = line
            header_idx = idx
            break
    if header is None:
        raise RuntimeError(f"Could not find history header in {path}")

    labels = []
    for token in header.replace("#", "").split():
        if token.startswith("[") and "]=" in token:
            labels.append(token.split("]=", 1)[1])

    data_rows = []
    for line in lines[header_idx + 1 :]:
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != len(labels):
            continue
        data_rows.append([float(x) for x in parts])
    if not data_rows:
        raise RuntimeError(f"No data rows in {path}")

    cols = {name: [] for name in labels}
    for row in data_rows:
        for i, name in enumerate(labels):
            cols[name].append(row[i])
    return cols


def find_mode_keys(prefix, cols):
    re_key = re.compile(rf"^{re.escape(prefix)}(\d+)$")
    keys = []
    for key in cols.keys():
        m = re_key.match(key)
        if m:
            keys.append((int(m.group(1)), key))
    keys.sort(key=lambda pair: pair[0])
    return keys


def parity_metrics(values_by_mode):
    total = 0.0
    even = 0.0
    odd = 0.0
    for n, value in values_by_mode.items():
        if not math.isfinite(value):
            continue
        mag = abs(value)
        total += mag
        if n % 2 == 0:
            even += mag
        else:
            odd += mag
    odd_frac = odd / total if total > 0.0 else math.nan
    even_frac = even / total if total > 0.0 else math.nan
    return (total, even, odd, even_frac, odd_frac)


def fmt(value):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "nan"
    if isinstance(value, str):
        return value
    return f"{value:.6e}"


def parse_csv_floats(raw):
    values = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(float(token))
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-modes", type=int, default=8, help="Number of modes to audit")
    parser.add_argument(
        "--n-quadrature",
        type=int,
        default=None,
        help="Quadrature order (default: n_modes+2)",
    )
    parser.add_argument("--lambda-value", type=float, default=1.0, help="Modes4D lambda")
    parser.add_argument(
        "--sigma-factors",
        default="0.5,0.75,1.0,1.5,2.0",
        help="Comma-separated gaussian sigma_factor values",
    )
    parser.add_argument(
        "--full-hst",
        default=None,
        help="Optional harris_full.out1.hst path for parity weighting diagnostics",
    )
    parser.add_argument(
        "--output-dir",
        default="projection_coeff_audit_outputs",
        help="Directory for coefficient + summary CSV outputs",
    )
    args = parser.parse_args()

    if args.n_modes <= 0:
        raise SystemExit("--n-modes must be > 0")
    n_quadrature = args.n_quadrature if args.n_quadrature is not None else (args.n_modes + 2)
    if n_quadrature < args.n_modes:
        raise SystemExit("--n-quadrature must be >= --n-modes")
    if args.lambda_value <= 0.0:
        raise SystemExit("--lambda-value must be > 0")

    sigma_factors = parse_csv_floats(args.sigma_factors)
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    nodes, weights, phi = build_mode_tables(args.n_modes, n_quadrature, args.lambda_value)

    kernel_specs = [("matched", math.nan), ("point", math.nan)]
    kernel_specs.extend(("gaussian", sigma) for sigma in sigma_factors)

    coeff_rows = []
    summary_rows = []

    hst_cols = None
    if args.full_hst:
        hst_cols = parse_hst(args.full_hst)

    for kernel, sigma in kernel_specs:
        sigma_factor = sigma if math.isfinite(sigma) else 1.0
        coeff = projection_coefficients(
            kernel=kernel,
            sigma_factor=sigma_factor,
            n_modes=args.n_modes,
            lambda_value=args.lambda_value,
            nodes=nodes,
            weights=weights,
            phi=phi,
        )

        coeff_values = {n: coeff[n] for n in range(args.n_modes)}
        coeff_total, coeff_even, coeff_odd, coeff_even_frac, coeff_odd_frac = parity_metrics(
            coeff_values
        )

        # Optional hst-weighted parity proxy: |M_n| * sqrt(E_n)
        weighted_total = math.nan
        weighted_even_frac = math.nan
        weighted_odd_frac = math.nan
        em_a2_even_frac = math.nan
        em_a2_odd_frac = math.nan
        jw_l2_even_frac = math.nan
        jw_l2_odd_frac = math.nan
        if hst_cols is not None:
            em_keys = find_mode_keys("m4d_em_a2_mode_", hst_cols)
            if em_keys:
                em_vals = {n: hst_cols[key][-1] for n, key in em_keys}
                _, _, _, em_a2_even_frac, em_a2_odd_frac = parity_metrics(em_vals)
                weighted = {}
                for n, value in em_vals.items():
                    amp = math.sqrt(max(value, 0.0)) if math.isfinite(value) else math.nan
                    weighted[n] = abs(coeff_values.get(n, 0.0)) * amp if math.isfinite(amp) else math.nan
                weighted_total, _, _, weighted_even_frac, weighted_odd_frac = parity_metrics(
                    weighted
                )

            jw_keys = find_mode_keys("m4d_jw_mode_l2_", hst_cols)
            if jw_keys:
                jw_vals = {n: hst_cols[key][-1] for n, key in jw_keys}
                _, _, _, jw_l2_even_frac, jw_l2_odd_frac = parity_metrics(jw_vals)

        row = {
            "kernel": kernel,
            "sigma_factor": sigma if math.isfinite(sigma) else math.nan,
            "coeff_l1_total": coeff_total,
            "coeff_l1_even": coeff_even,
            "coeff_l1_odd": coeff_odd,
            "coeff_even_fraction": coeff_even_frac,
            "coeff_odd_fraction": coeff_odd_frac,
            "weighted_coeff_energy_proxy_total": weighted_total,
            "weighted_coeff_energy_proxy_even_fraction": weighted_even_frac,
            "weighted_coeff_energy_proxy_odd_fraction": weighted_odd_frac,
            "full_hst_em_a2_even_fraction": em_a2_even_frac,
            "full_hst_em_a2_odd_fraction": em_a2_odd_frac,
            "full_hst_jw_l2_even_fraction": jw_l2_even_frac,
            "full_hst_jw_l2_odd_fraction": jw_l2_odd_frac,
        }
        for n in range(args.n_modes):
            row[f"coeff_mode_{n}"] = coeff[n]
        summary_rows.append(row)

        for n in range(args.n_modes):
            coeff_rows.append(
                {
                    "kernel": kernel,
                    "sigma_factor": sigma if math.isfinite(sigma) else math.nan,
                    "mode": n,
                    "coeff": coeff[n],
                }
            )

    coeff_csv = output_dir / "projection_coefficients_by_mode.csv"
    with coeff_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=["kernel", "sigma_factor", "mode", "coeff"])
        writer.writeheader()
        writer.writerows(coeff_rows)

    summary_fieldnames = list(summary_rows[0].keys())
    summary_csv = output_dir / "projection_coefficients_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=summary_fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(",".join(summary_fieldnames))
    for row in summary_rows:
        print(",".join(fmt(row.get(k, "")) for k in summary_fieldnames))
    print(f"coeff_csv,{coeff_csv}")
    print(f"summary_csv,{summary_csv}")


if __name__ == "__main__":
    main()
