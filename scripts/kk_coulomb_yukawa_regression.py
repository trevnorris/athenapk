#!/usr/bin/env python3
"""Validate Gaussian KK Coulomb+Yukawa tower behavior."""

import argparse
import math


def coupling_ratio_even(m):
    return math.comb(2 * m, m) / (4**m)


def c_n(n, lam):
    c0 = 1.0 / (lam * math.sqrt(math.pi))
    if n % 2 == 1:
        return 0.0
    return c0 * coupling_ratio_even(n // 2)


def m_n(n, lam):
    return math.sqrt(2.0 * n) / lam


def normalized_potential(r, lam, n_terms):
    c0 = c_n(0, lam)
    tower = 0.0
    for n in range(n_terms):
        tower += c_n(n, lam) * math.exp(-m_n(n, lam) * r)
    return tower / c0


def validate_lambda(
    lam,
    n_terms,
    r1_scale,
    r2_scale,
    mass_rel_tol,
    amp_rel_tol,
):
    expected_mass = 2.0 / lam
    expected_amp = 0.5
    c0 = c_n(0, lam)

    # Check documented Gaussian KK pattern.
    mode_pattern_max_rel_err = 0.0
    odd_max_abs = 0.0
    for n in range(n_terms):
        cn = c_n(n, lam)
        if n % 2 == 1:
            odd_max_abs = max(odd_max_abs, abs(cn))
            continue
        m = n // 2
        expected = c0 * coupling_ratio_even(m)
        if expected > 0.0:
            rel_err = abs(cn - expected) / expected
            mode_pattern_max_rel_err = max(mode_pattern_max_rel_err, rel_err)

    r1 = r1_scale * lam
    r2 = r2_scale * lam
    delta1 = normalized_potential(r1, lam, n_terms) - 1.0
    delta2 = normalized_potential(r2, lam, n_terms) - 1.0
    if delta1 <= 0.0 or delta2 <= 0.0:
        return {
            "lambda": lam,
            "status": "FAIL",
            "failures": "nonpositive_delta",
        }

    mass_est = math.log(delta1 / delta2) / (r2 - r1)
    amp_est = delta1 * math.exp(mass_est * r1)
    mass_rel_err = abs(mass_est - expected_mass) / expected_mass
    amp_rel_err = abs(amp_est - expected_amp) / expected_amp

    # Require tail suppression at large distance in r/lambda units.
    far_delta = normalized_potential(8.0 * lam, lam, n_terms) - 1.0

    failures = []
    if odd_max_abs > 0.0:
        failures.append("odd_couplings")
    if mode_pattern_max_rel_err > 0.0:
        failures.append("mode_pattern")
    if mass_rel_err > mass_rel_tol:
        failures.append("yukawa_mass")
    if amp_rel_err > amp_rel_tol:
        failures.append("yukawa_amplitude")
    if far_delta > 1.0e-6:
        failures.append("far_tail")

    status = "FAIL" if failures else "PASS"
    return {
        "lambda": lam,
        "mass_est": mass_est,
        "mass_expected": expected_mass,
        "mass_rel_err": mass_rel_err,
        "amp_est": amp_est,
        "amp_expected": expected_amp,
        "amp_rel_err": amp_rel_err,
        "delta_r1": delta1,
        "delta_r2": delta2,
        "delta_far": far_delta,
        "odd_max_abs": odd_max_abs,
        "mode_pattern_max_rel_err": mode_pattern_max_rel_err,
        "status": status,
        "failures": "|".join(failures) if failures else "none",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lambdas",
        default="0.5,1.0,2.0",
        help="Comma-separated lambda values",
    )
    parser.add_argument(
        "--n-terms",
        type=int,
        default=64,
        help="Number of KK terms in the truncated tower",
    )
    parser.add_argument(
        "--r1-scale",
        type=float,
        default=4.0,
        help="First distance probe in units of lambda",
    )
    parser.add_argument(
        "--r2-scale",
        type=float,
        default=5.0,
        help="Second distance probe in units of lambda",
    )
    parser.add_argument(
        "--mass-rel-tol",
        type=float,
        default=0.05,
        help="Relative tolerance for inferred leading Yukawa mass",
    )
    parser.add_argument(
        "--amp-rel-tol",
        type=float,
        default=0.10,
        help="Relative tolerance for inferred leading Yukawa amplitude",
    )
    args = parser.parse_args()

    lam_values = [float(x.strip()) for x in args.lambdas.split(",") if x.strip()]
    if args.n_terms < 3:
        raise SystemExit("--n-terms must be >= 3")
    if args.r2_scale <= args.r1_scale:
        raise SystemExit("--r2-scale must be > --r1-scale")

    rows = []
    for lam in lam_values:
        if lam <= 0.0:
            raise SystemExit("All lambda values must be > 0")
        rows.append(
            validate_lambda(
                lam,
                args.n_terms,
                args.r1_scale,
                args.r2_scale,
                args.mass_rel_tol,
                args.amp_rel_tol,
            )
        )

    keys = [
        "lambda",
        "mass_est",
        "mass_expected",
        "mass_rel_err",
        "amp_est",
        "amp_expected",
        "amp_rel_err",
        "delta_r1",
        "delta_r2",
        "delta_far",
        "odd_max_abs",
        "mode_pattern_max_rel_err",
        "status",
        "failures",
    ]
    print(",".join(keys))
    for row in rows:
        vals = []
        for k in keys:
            v = row.get(k)
            if isinstance(v, float):
                vals.append(f"{v:.6e}")
            else:
                vals.append(str(v))
        print(",".join(vals))

    failing = [r for r in rows if r["status"] == "FAIL"]
    if failing:
        raise SystemExit(
            "kk coulomb+yukawa regression failed for lambda: "
            + ", ".join(f"{r['lambda']:.6g}" for r in failing)
        )


if __name__ == "__main__":
    main()
