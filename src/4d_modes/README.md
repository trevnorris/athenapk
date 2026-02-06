# 4D Modes Scaffold (`src/4d_modes`)

This directory contains the initial AthenaPK integration scaffold for the 3D x modes 4+1D solver described in `docs/plan.md`.

## Current scope

Implemented now:
- `ModeTables` for Hermite-Gaussian mode math (`mode_tables.cpp`)
- scaled Gauss-Hermite quadrature (`w_q = lambda * x_q`, `omega_q = lambda * omega_hat_q`)
- mode/node transforms and identity helpers (ID-3, ID-4, ID-7)
- package registration skeleton (`modes4d`) with EM/plasma field registration
- unsplit EM mode source update using the semi-discrete equations in
  `docs/mode_equations_verified.md`:
  - brane components include Laplacian, `m_n^2` term, nearest-neighbor `A_w` coupling,
    and `J^nu` source deposition
  - scalar component includes Laplacian and `J^w` source deposition
  - scalar-photon work (`J^w E_w`) and leakage (`S_leak`) accumulators are updated
- split plasma transverse-momentum (`momw`) source update:
  - reconstruct `E_w` and `C_a` at quadrature nodes from EM modes
  - apply species force `q/m * rho * (E_w - v^a C_a)` in node space
  - project updated `momw` back to mode coefficients
- split plasma brane-momentum + energy source update:
  - apply brane Lorentz terms `q/m * rho * (E + v x B + v_w C)` for `momx/momy/momz`
  - apply species work `q/m * rho * (v·E + v_w E_w)` for mode energy
  - project updated node values back to mode coefficients
- split plasma mode-continuity coupling update:
  - source-step term applies leakage coupling
    `d_t rho^(n) = -(sqrt(2(n+1))/lambda) j_w^(n+1)`
  - brane transport `div(j^a,(n))` is handled in the conservative flux pipeline by
    registering `plasma4d_cons` with `WithFluxes` and filling stage fluxes for
    `rho/momx/momy/momz/momw/energy` each stage
  - runtime gain `plasma_rho_divj_gain` scales that conservative brane-transport path
  - write updated `plasma4d_cons` back each source step
- `harris_4d` problem hook and zero-mode Harris initialization scaffold
- `harris_4d` now uses the same unsplit equation-based EM source update during runtime
- standalone unit-test executable for mode math (`modes4d_unit_tests`)
- initial history diagnostics:
  - `m4d_em_a2`
  - `m4d_em_pi2`
  - `m4d_plasma_cons2`
  - `m4d_int_jw_ew`
  - `m4d_int_s_leak`
  - `m4d_int_s_leak_abs`
  - `m4d_int_divj_mode0`
  - `m4d_brane_e2`
  - `m4d_brane_b2`
  - `m4d_brane_epar2`
  - `m4d_mixed_ew2`
  - `m4d_mixed_c2`
  - `m4d_psi0_span`
  - `m4d_pulse_xc`
  - `m4d_em_a2_mode_<n>`
  - `m4d_em_pi2_mode_<n>`
  - `m4d_jw_mode_l2_<n>`
  - `m4d_jw_mode_<n>`
  - `m4d_charge_mode_<n>`
  - local continuity closure channels (source-step residuals):
    - `m4d_cont_local_l1`
    - `m4d_cont_local_l2`
    - `m4d_cont_local_max_abs`
    - `m4d_cont_mode0_l1`
    - `m4d_cont_mode0_l2`
    - `m4d_cont_mode0_max_abs`

Not implemented yet:
- full two-fluid mode dynamics beyond the current source + advective-transport bring-up
- full conservative EM update in the AthenaPK flux pipeline (current path is source-step)
- full reconnection workflow/analysis

Additional bring-up path now available:
- `inputs/em4d_pulse.in` exercises EM-only source evolution for `em4d_a`/`em4d_pi`
- `inputs/harris_4d_controlled.in` and `inputs/harris_4d_full.in` provide a
  controlled-vs-full scan pair (`n_modes=1` vs `n_modes=4`), with full-case
  mixed-channel activation driven by plasma response (no direct `J^w` mode seed)

## Files

- `mode_tables.hpp`, `mode_tables.cpp`: Hermite basis, quadrature, transforms, sparse operator helpers
- `mode_transform.hpp`, `mode_transform.cpp`: forward/inverse transform wrappers
- `package4d_modes.hpp`, `package4d_modes.cpp`: `modes4d` package setup and parameters
- `em4d_modes.*`: registered EM mode fields (`em4d_a`, `em4d_pi`)
- `plasma4d_modes.*`: registered plasma mode fields (`plasma4d_cons`)
  plus conservative plasma-transport flux helpers
- `diagnostics4d.*`: placeholder diagnostics params
- `pgen_harris4d.*`: Harris-sheet initialization for the `harris_4d` problem path
- `mode_tables_tests.cpp`: standalone mode-math tests

## Build options

At AthenaPK configure time:

- `-DAthenaPK_ENABLE_4D_MODES=ON|OFF`
- `-DAthenaPK_ENABLE_4D_MODES_UNIT_TESTS=ON|OFF`

Defaults are `ON` for both.

## Quick start

### 1) Configure and build

```bash
cmake -S /projects/fluid-engine/athenapk -B /projects/fluid-engine/athenapk/build-baseline \
  -DPARTHENON_DISABLE_MPI=ON \
  -DPARTHENON_DISABLE_HDF5=ON \
  -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF \
  -DPARTHENON_ENABLE_TESTING=OFF \
  -DAthenaPK_ENABLE_4D_MODES=ON \
  -DAthenaPK_ENABLE_4D_MODES_UNIT_TESTS=ON

cmake --build /projects/fluid-engine/athenapk/build-baseline -j$(nproc)
```

### 2) Run mode-math unit tests

```bash
/projects/fluid-engine/athenapk/build-baseline/bin/modes4d_unit_tests
```

Expected output:

```text
modes4d_unit_tests: PASS
```

### 3) Run Harris scaffold input

```bash
/projects/fluid-engine/athenapk/build-baseline/bin/athenaPK \
  -i /projects/fluid-engine/athenapk/inputs/harris_4d.in
```

This input deck includes a history output stream (`file_type = hst`) so the
`m4d_*` diagnostics channels are emitted during the run.
For reconnection-oriented analysis, this includes `m4d_brane_epar2` (parallel
electric proxy), mixed-sector proxies (`m4d_mixed_ew2`, `m4d_mixed_c2`), and
`m4d_psi0_span` as a zero-mode flux-function span proxy.

### 4) Run EM-only pulse bring-up

```bash
/projects/fluid-engine/athenapk/build-baseline/bin/athenaPK \
  -i /projects/fluid-engine/athenapk/inputs/em4d_pulse.in
```

This path enables a minimal unsplit source evolution:
- equation-based mode update for `em4d_a`/`em4d_pi` (semi-discrete Maxwell form),
  with current sources derived from `plasma4d_cons`

with `c_wave` and `damping` configured in `<modes4d>` as:
- `em_c_wave`
- `em_damping`
- `em_mu0`
- `plasma_qom_ion`
- `plasma_qom_electron`
- `plasma_force_source_gain`
- `plasma_momw_source_gain`
- `plasma_momw_damping`
- `plasma_rho_divj_gain`
- `plasma_rho_floor`
- `plasma_energy_source_gain`
- `plasma_energy_floor`

### 5) Run controlled-vs-full Harris scan summary

```bash
python3 /projects/fluid-engine/athenapk/scripts/harris_scan_matrix.py \
  --binary /projects/fluid-engine/athenapk/build-baseline/bin/athenaPK \
  --workdir /projects/fluid-engine \
  --controlled-input /projects/fluid-engine/athenapk/inputs/harris_4d_controlled.in \
  --full-input /projects/fluid-engine/athenapk/inputs/harris_4d_full.in \
  --output-dir /projects/fluid-engine/athenapk/build-baseline/harris_scan_outputs
```

This prints a CSV table with final values and Pearson correlations between
`m4d_psi0_span` and key ledger/reconnection proxy channels.
It also reports continuity-closure metrics based on
`m4d_charge_mode_0`, `m4d_int_s_leak`, and `m4d_int_divj_mode0`, with `closure_status`
computed from a normalized residual threshold (`--closure-norm-tol`, default `5e-2`).
Low-signal cases use `--closure-abs-rate-tol` (default `1e-8`) as an absolute-rate gate.
The summary also includes local mode-0 closure residual channels from
`m4d_cont_mode0_*`, with pass/fail controlled by
`--closure-local-mode0-abs-rate-tol` (default `1e-8`).
`m4d_int_divj_mode0` is accumulated from the corrected conservative rho-transport fluxes.
Use `--fail-on-check` to force a nonzero exit code on failed checks.

## Input requirements for `harris_4d`

Required right now:
- `<job> problem_id = harris_4d`
- `<hydro> fluid = glmmhd`
- `<modes4d> enabled = true`

Key parameters currently used from `<problem/harris_4d>`:
- `b0`
- `guide_bz`
- `n_bg`
- `n_sheet`
- `p_bg`
- `sheet_half_width`
- `perturbation_amp`
- `drift_current_scale`
- `aw_mode1_amp`
- `piw_mode1_amp`

## Notes

- The `harris_4d` initializer currently sets:
  - hydro/MHD conservative state (`cons`) for a Harris-like zero-mode profile
  - mode-0 `A_y` perturbation in `em4d_a`
  - optional mode-1 `A_w` perturbation in `em4d_a` for mixed-sector triggering
  - species mode-0 plasma state in `plasma4d_cons` using q/m-aware splits that
    preserve mode-0 charge neutrality and target drift-current amplitude
- This is a bring-up scaffold and not yet the final physics model from the docs.
