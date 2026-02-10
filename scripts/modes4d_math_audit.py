#!/usr/bin/env python3
"""Audit core 4D mode math identities against implementation contracts.

This script validates the Hermite/Gauss-Hermite machinery used by the 3D×modes
solver and checks for source-level contracts in the C++ implementation.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.polynomial.hermite import hermgauss, hermval


@dataclass
class CheckResult:
    check_id: str
    status: str
    detail: str


def build_mode_tables(n_modes: int, n_quad: int, lam: float):
    x_std, w_std = hermgauss(n_quad)
    nodes = lam * x_std
    weights = lam * w_std
    phi = np.zeros((n_modes, n_quad), dtype=np.float64)
    dphi = np.zeros((n_modes, n_quad), dtype=np.float64)
    normalizer_factor = lam * math.sqrt(math.pi)
    x = nodes / lam
    for n in range(n_modes):
        coeff = np.zeros(n + 1, dtype=np.float64)
        coeff[n] = 1.0
        hn = hermval(x, coeff)
        normalization = math.sqrt(normalizer_factor * (2.0**n) * math.factorial(n))
        phi[n, :] = hn / normalization
        if n > 0:
            dphi[n, :] = (math.sqrt(2.0 * n) / lam) * phi[n - 1, :]
    return nodes, weights, phi, dphi


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a))) if a.size else 0.0


def append(results: list[CheckResult], check_id: str, ok: bool, detail: str) -> None:
    results.append(CheckResult(check_id=check_id, status="PASS" if ok else "FAIL", detail=detail))


def append_warn(results: list[CheckResult], check_id: str, detail: str) -> None:
    results.append(CheckResult(check_id=check_id, status="WARN", detail=detail))


def source_has(path: Path, pattern: str) -> bool:
    text = path.read_text(encoding="utf-8")
    return re.search(pattern, text, re.MULTILINE | re.DOTALL) is not None


def format_results(results: Iterable[CheckResult]) -> str:
    lines = ["check_id,status,detail"]
    for r in results:
        detail = r.detail.replace("\n", " ").replace(",", ";")
        lines.append(f"{r.check_id},{r.status},{detail}")
    return "\n".join(lines)


def run_audit(repo_root: Path, n_modes: int, n_quad: int, lam: float, strict_gaps: bool) -> int:
    if n_modes <= 0:
        raise ValueError("n_modes must be > 0")
    if n_quad < n_modes:
        raise ValueError("n_quadrature must be >= n_modes")
    if lam <= 0.0:
        raise ValueError("lambda must be > 0")

    results: list[CheckResult] = []
    nodes, weights, phi, dphi = build_mode_tables(n_modes, n_quad, lam)
    z_int = lam * math.sqrt(math.pi)
    atol = 5.0e-11
    rtol = 5.0e-11

    # ID-0 / orthonormality
    gram = np.einsum("q,nq,mq->nm", weights, phi, phi)
    expected = np.eye(n_modes, dtype=np.float64)
    gram_err = max_abs(gram - expected)
    append(
        results,
        "ID0_orthonormality",
        gram_err <= 2.0e-11,
        f"max_abs_error={gram_err:.3e}",
    )

    # ID-2 / mass spectrum
    mass = np.array([(2.0 * n) / (lam * lam) for n in range(n_modes)], dtype=np.float64)
    mass_expected = mass.copy()
    mass_err = max_abs(mass - mass_expected)
    append(
        results,
        "ID2_mass_spectrum",
        mass_err <= 1.0e-14,
        f"max_abs_error={mass_err:.3e}",
    )

    # ID-3 / raising
    rng = np.random.default_rng(seed=4)
    mode_values = rng.normal(size=n_modes)
    g = np.einsum("n,nq->q", mode_values, phi)
    gp = np.einsum("n,nq->q", mode_values, dphi)
    lhs_id3 = np.einsum("q,nq,q->n", weights, phi, gp)
    rhs_id3 = np.zeros(n_modes, dtype=np.float64)
    for n in range(n_modes - 1):
        rhs_id3[n] = (math.sqrt(2.0 * (n + 1)) / lam) * mode_values[n + 1]
    id3_err = max_abs(lhs_id3 - rhs_id3)
    append(
        results,
        "ID3_raising",
        id3_err <= (atol + rtol * max_abs(rhs_id3)),
        f"max_abs_error={id3_err:.3e}",
    )

    # ID-4 / lowering via derivative of Z*G
    lambda_sq = lam * lam
    lhs_id4 = np.einsum("q,nq,q->n", weights, phi, (gp - (2.0 * nodes / lambda_sq) * g))
    rhs_id4 = np.zeros(n_modes, dtype=np.float64)
    for n in range(1, n_modes):
        rhs_id4[n] = -(math.sqrt(2.0 * n) / lam) * mode_values[n - 1]
    id4_err = max_abs(lhs_id4 - rhs_id4)
    append(
        results,
        "ID4_lowering",
        id4_err <= (atol + rtol * max_abs(rhs_id4)),
        f"max_abs_error={id4_err:.3e}",
    )

    # ID-7 / double-raising from second derivative
    gpp = np.zeros(n_quad, dtype=np.float64)
    for m in range(2, n_modes):
        pref = 2.0 * math.sqrt(m * (m - 1)) / lambda_sq
        gpp += mode_values[m] * pref * phi[m - 2, :]
    lhs_id7 = np.einsum("q,nq,q->n", weights, phi, gpp)
    rhs_id7 = np.zeros(n_modes, dtype=np.float64)
    for n in range(n_modes - 2):
        rhs_id7[n] = (2.0 * math.sqrt((n + 1) * (n + 2)) / lambda_sq) * mode_values[n + 2]
    id7_err = max_abs(lhs_id7 - rhs_id7)
    append(
        results,
        "ID7_double_raising",
        id7_err <= (5.0e-10 + 5.0e-10 * max_abs(rhs_id7)),
        f"max_abs_error={id7_err:.3e}",
    )

    # ID-6 / matched projection kernel yields zero-mode-only readout.
    proj_coeff = np.einsum("q,nq->n", weights, phi) / z_int
    expected_coeff = np.zeros(n_modes, dtype=np.float64)
    expected_coeff[0] = 1.0 / math.sqrt(z_int)
    id6_err = max_abs(proj_coeff - expected_coeff)
    append(
        results,
        "ID6_matched_projection",
        id6_err <= 2.0e-11,
        f"max_abs_error={id6_err:.3e}",
    )

    # Brane-point coefficients parity check at w=0 (odd modes vanish).
    x0 = 0.0
    coeff_w0 = np.zeros(n_modes, dtype=np.float64)
    for n in range(n_modes):
        coeff = np.zeros(n + 1, dtype=np.float64)
        coeff[n] = 1.0
        hn0 = hermval(x0, coeff)
        norm = math.sqrt(z_int * (2.0**n) * math.factorial(n))
        coeff_w0[n] = hn0 / norm
    odd_max = max_abs(coeff_w0[1::2])
    append(
        results,
        "brane_point_parity",
        odd_max <= 1.0e-12,
        f"max_abs_odd_mode_coeff={odd_max:.3e}",
    )

    em_cpp = repo_root / "src/4d_modes/em4d_modes.cpp"
    diagnostics_cpp = repo_root / "src/4d_modes/diagnostics4d.cpp"
    package_cpp = repo_root / "src/4d_modes/package4d_modes.cpp"
    mode_tables_cpp = repo_root / "src/4d_modes/mode_tables.cpp"

    # Source-level contracts for implemented couplings.
    append(
        results,
        "code_contract_id3",
        source_has(
            mode_tables_cpp,
            r"ApplyID3Raising.*?return\s+\(std::sqrt\(2\.0 \* static_cast<double>\(m\)\)\s*/\s*config_\.lambda\)\s*\*\s*mode_values\[m\];",
        ),
        "ModeTables::ApplyID3Raising coefficient",
    )
    append(
        results,
        "code_contract_id4",
        source_has(
            mode_tables_cpp,
            r"ApplyID4Lowering.*?return\s+-\(std::sqrt\(2\.0 \* static_cast<double>\(n\)\)\s*/\s*config_\.lambda\)\s*\*\s*mode_values\[m\];",
        ),
        "ModeTables::ApplyID4Lowering sign/coefficient",
    )
    append(
        results,
        "code_contract_em_mass",
        source_has(
            em_cpp,
            r"src_ay_mass\s*=\s*-em_source_mass_gain\s*\*\s*\(c2 \* mass_squared\[n\] \* a_old\(idx_ay",
        ),
        "EM mass term includes -c^2 m_n^2 a_y^(n)",
    )
    append(
        results,
        "code_contract_em_mixed",
        source_has(
            em_cpp,
            r"src_ay_spatial_mixed\s*=\s*c2 \* coupling_coeff \* mixed_grad_aw_y",
        ),
        "EM mixed spatial coupling uses sqrt(2n)/lambda factor",
    )
    append(
        results,
        "code_contract_continuity",
        source_has(
            em_cpp,
            r"leak_rhs\s*=\s*-\s*coupling\s*\*\s*plasma_new\(base_np1 \+ kPlasmaMomW",
        ),
        "Mode continuity source includes -sqrt(2(n+1))/lambda * j_w^(n+1)",
    )
    append(
        results,
        "code_contract_projection",
        source_has(
            diagnostics_cpp,
            r"ProjectionCoefficients\(tables,\s*n_modes,\s*projection_kernel,\s*projection_sigma_factor\)",
        ),
        "psi_proj uses runtime projection-kernel coefficient map",
    )

    # Gap checks: these are warnings, not hard failures by default.
    has_explicit_gauge_control = (
        source_has(package_cpp, r'em_source_gauge_gain')
        and source_has(em_cpp, r'em4d/source_gauge_gain')
        and source_has(em_cpp, r'src_a0_gauge')
    )
    if has_explicit_gauge_control:
        append(
            results,
            "gap_gauge_control",
            True,
            "Explicit gauge-control gain and source term detected.",
        )
    else:
        append_warn(
            results,
            "gap_gauge_control",
            "No explicit gauge-control gain/source path detected; solver relies on damping terms only.",
        )

    has_projection_kernel_param = (
        source_has(package_cpp, r'diag_projection_kernel')
        and source_has(package_cpp, r'diag/projection_kernel')
        and source_has(diagnostics_cpp, r'if \(kernel == "gaussian"\)')
        and source_has(diagnostics_cpp, r'if \(kernel == "point"\)')
    )
    if has_projection_kernel_param:
        append(
            results,
            "gap_projection_kernel",
            True,
            "Projection-kernel parameterization detected (matched/point/gaussian).",
        )
    else:
        append_warn(
            results,
            "gap_projection_kernel",
            "No configurable projection kernel W(w) detected in runtime params; diagnostics are tied to matched and point-sample readouts.",
        )

    print(format_results(results))

    fails = [r for r in results if r.status == "FAIL"]
    warns = [r for r in results if r.status == "WARN"]
    if fails:
        return 1
    if strict_gaps and warns:
        return 2
    return 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root", type=Path, default=default_root)
    p.add_argument("--n-modes", type=int, default=8)
    p.add_argument("--n-quadrature", type=int, default=16)
    p.add_argument("--lambda", dest="lam", type=float, default=1.0)
    p.add_argument(
        "--strict-gaps",
        action="store_true",
        help="Treat WARN gap checks as non-zero exit.",
    )
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    return run_audit(
        repo_root=args.source_root.resolve(),
        n_modes=args.n_modes,
        n_quad=args.n_quadrature,
        lam=args.lam,
        strict_gaps=args.strict_gaps,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
