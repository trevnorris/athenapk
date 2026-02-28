#!/usr/bin/env python3
"""Deterministic projection-operator verification for Modes4D diagnostics.

This validates that diagnostic projection kernels match their intended observer
definitions:
- matched: integral with W = Z / Z_int
- point:   evaluation at w = 0
- gaussian: integral with normalized gaussian kernel (sigma_factor)
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss, hermval


@dataclass
class CheckResult:
    check_id: str
    status: str
    detail: str


def append(results: list[CheckResult], check_id: str, ok: bool, detail: str) -> None:
    results.append(CheckResult(check_id=check_id, status="PASS" if ok else "FAIL", detail=detail))


def format_results(results: list[CheckResult]) -> str:
    lines = ["check_id,status,detail"]
    for r in results:
        lines.append(f"{r.check_id},{r.status},{r.detail.replace(',', ';')}")
    return "\n".join(lines)


def hermite_physicists_values(max_n: int, x: np.ndarray) -> np.ndarray:
    out = np.zeros((max_n + 1, x.size), dtype=np.float64)
    out[0, :] = 1.0
    if max_n >= 1:
        out[1, :] = 2.0 * x
    for n in range(1, max_n):
        out[n + 1, :] = (2.0 * x * out[n, :]) - (2.0 * n * out[n - 1, :])
    return out


def build_mode_tables(n_modes: int, n_quad: int, lam: float):
    x_std, w_std = hermgauss(n_quad)
    nodes = lam * x_std
    weights = lam * w_std
    phi = np.zeros((n_modes, n_quad), dtype=np.float64)
    normalizer_factor = lam * math.sqrt(math.pi)
    x = nodes / lam
    for n in range(n_modes):
        coeff = np.zeros(n + 1, dtype=np.float64)
        coeff[n] = 1.0
        hn = hermval(x, coeff)
        normalization = math.sqrt(normalizer_factor * (2.0**n) * math.factorial(n))
        phi[n, :] = hn / normalization
    return nodes, weights, phi


def mode_basis_values_at_w(n_modes: int, lam: float, w: np.ndarray) -> np.ndarray:
    x = w / lam
    hn = hermite_physicists_values(n_modes - 1, x)
    z_int = lam * math.sqrt(math.pi)
    phi = np.zeros((n_modes, w.size), dtype=np.float64)
    for n in range(n_modes):
        norm = math.sqrt(z_int * (2.0**n) * math.factorial(n))
        phi[n, :] = hn[n, :] / norm
    return phi


def profile_from_modes(coeff: np.ndarray, phi_vals: np.ndarray) -> np.ndarray:
    return np.einsum("n,nq->q", coeff, phi_vals)


def projection_coefficients(
    kernel: str,
    sigma_factor: float,
    lam: float,
    nodes: np.ndarray,
    weights: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    n_modes = phi.shape[0]
    z_int = lam * math.sqrt(math.pi)

    if kernel == "point":
        # phi_n(0) via recurrence-compatible Hermite evaluation
        coeff = np.zeros(n_modes, dtype=np.float64)
        for n in range(n_modes):
            p = np.zeros(n + 1, dtype=np.float64)
            p[n] = 1.0
            hn0 = float(hermval(0.0, p))
            norm = math.sqrt(z_int * (2.0**n) * math.factorial(n))
            coeff[n] = hn0 / norm
        return coeff

    if kernel == "matched":
        return np.einsum("q,nq->n", weights, phi) / z_int

    if kernel == "gaussian":
        lam_sq = lam * lam
        sig_sq = sigma_factor * sigma_factor
        w_norm = sigma_factor * z_int
        z = np.exp(-(nodes * nodes) / lam_sq)
        w_kernel = np.exp(-(nodes * nodes) / (sig_sq * lam_sq)) / w_norm
        weight_factor = np.where(z > 0.0, w_kernel / z, 0.0)
        return np.einsum("q,nq,q->n", weights, phi, weight_factor)

    raise ValueError(f"Unsupported kernel: {kernel}")


def matched_kernel(w: np.ndarray, lam: float) -> np.ndarray:
    z_int = lam * math.sqrt(math.pi)
    return np.exp(-(w * w) / (lam * lam)) / z_int


def gaussian_kernel(w: np.ndarray, lam: float, sigma: float) -> np.ndarray:
    z_int = lam * math.sqrt(math.pi)
    return np.exp(-(w * w) / ((sigma * lam) * (sigma * lam))) / (sigma * z_int)


def trapz_integral(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


def source_has(path: Path, pattern: str) -> bool:
    text = path.read_text(encoding="utf-8")
    return re.search(pattern, text, re.MULTILINE | re.DOTALL) is not None


def run_tests(source_root: Path, n_modes: int, n_quad: int, lam: float) -> int:
    if n_modes <= 0:
        raise ValueError("n_modes must be > 0")
    if n_quad < n_modes:
        raise ValueError("n_quadrature must be >= n_modes")
    if lam <= 0.0:
        raise ValueError("lambda must be > 0")

    results: list[CheckResult] = []
    nodes, weights, phi = build_mode_tables(n_modes=n_modes, n_quad=n_quad, lam=lam)

    # Dense grid reference for observer-kernel integrals.
    w_ref = np.linspace(-10.0 * lam, 10.0 * lam, 20001, dtype=np.float64)
    phi_ref = mode_basis_values_at_w(n_modes=n_modes, lam=lam, w=w_ref)
    phi0 = mode_basis_values_at_w(n_modes=n_modes, lam=lam, w=np.array([0.0], dtype=np.float64))[
        :, 0
    ]

    profiles = {
        "profile_even_mix": np.array(
            [0.83, 0.0, -0.47, 0.0, 0.19, 0.0, -0.08, 0.0], dtype=np.float64
        )[:n_modes],
        "profile_even_odd_mix": np.array(
            [0.71, -0.23, 0.31, -0.15, 0.07, -0.03, 0.02, -0.01], dtype=np.float64
        )[:n_modes],
        "profile_highmode_mix": np.array(
            [0.0, 0.0, 0.0, 0.0, 0.42, -0.21, 0.13, -0.06], dtype=np.float64
        )[:n_modes],
    }

    # Matched and point checks for each synthetic profile.
    tol_integral = 7.0e-6
    tol_point = 1.0e-11
    for profile_name, mode_coeff in profiles.items():
        q_ref = profile_from_modes(mode_coeff, phi_ref)

        m_matched = projection_coefficients(
            kernel="matched", sigma_factor=1.0, lam=lam, nodes=nodes, weights=weights, phi=phi
        )
        got_matched = float(np.dot(m_matched, mode_coeff))
        exp_matched = trapz_integral(matched_kernel(w_ref, lam) * q_ref, w_ref)
        err_matched = abs(got_matched - exp_matched)
        append(
            results,
            f"projection_matched_{profile_name}",
            err_matched <= tol_integral,
            f"abs_error={err_matched:.3e}; got={got_matched:.6e}; expected={exp_matched:.6e}",
        )

        m_point = projection_coefficients(
            kernel="point", sigma_factor=1.0, lam=lam, nodes=nodes, weights=weights, phi=phi
        )
        got_point = float(np.dot(m_point, mode_coeff))
        exp_point = float(np.dot(phi0, mode_coeff))
        err_point = abs(got_point - exp_point)
        append(
            results,
            f"projection_point_{profile_name}",
            err_point <= tol_point,
            f"abs_error={err_point:.3e}; got={got_point:.6e}; expected={exp_point:.6e}",
        )

    # Gaussian checks for representative sigma values.
    gaussian_sigmas = (0.5, 0.75, 1.0, 1.5)
    gaussian_profile = profiles["profile_even_odd_mix"]
    q_ref = profile_from_modes(gaussian_profile, phi_ref)
    for sigma in gaussian_sigmas:
        m_gauss = projection_coefficients(
            kernel="gaussian",
            sigma_factor=sigma,
            lam=lam,
            nodes=nodes,
            weights=weights,
            phi=phi,
        )
        got_gauss = float(np.dot(m_gauss, gaussian_profile))
        exp_gauss = trapz_integral(gaussian_kernel(w_ref, lam, sigma) * q_ref, w_ref)
        err_gauss = abs(got_gauss - exp_gauss)
        append(
            results,
            f"projection_gaussian_sigma_{str(sigma).replace('.', 'p')}",
            err_gauss <= tol_integral,
            f"abs_error={err_gauss:.3e}; got={got_gauss:.6e}; expected={exp_gauss:.6e}",
        )

    # Code-level contracts for point-projection equivalence semantics.
    diagnostics_cpp = source_root / "src/4d_modes/diagnostics4d.cpp"
    has_point_projection_coeffs = source_has(
        diagnostics_cpp,
        r'if \(kernel == "point"\)\s*\{\s*return BranePointCoefficients\(tables,\s*n_modes,\s*center_w\);\s*\}',
    )
    has_psiw0_point_coeffs = source_has(
        diagnostics_cpp,
        r"Real AyBranePointSpanWithParityHst\(.*?BranePointCoefficients\(tables,\s*n_modes,\s*center_w\)",
    )
    append(
        results,
        "code_contract_point_kernel_equals_w0_coeff_basis",
        has_point_projection_coeffs and has_psiw0_point_coeffs,
        "point projection and psi_w0 both use BranePointCoefficients",
    )

    has_matched_definition = source_has(
        diagnostics_cpp,
        r"Default:\s*matched kernel W = Z\(w-center_w\) / Z_int\.",
    )
    append(
        results,
        "code_contract_matched_kernel_definition",
        has_matched_definition,
        "matched observer kernel documented as W=Z/Z_int in diagnostics code",
    )

    print(format_results(results))
    return 1 if any(r.status == "FAIL" for r in results) else 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root", type=Path, default=default_root)
    p.add_argument("--n-modes", type=int, default=8)
    p.add_argument("--n-quadrature", type=int, default=24)
    p.add_argument("--lambda", dest="lam", type=float, default=1.0)
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    return run_tests(
        source_root=args.source_root.resolve(),
        n_modes=args.n_modes,
        n_quad=args.n_quadrature,
        lam=args.lam,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
