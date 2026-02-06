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
- `harris_4d` problem hook and zero-mode Harris initialization scaffold
- standalone unit-test executable for mode math (`modes4d_unit_tests`)
- initial history diagnostics:
  - `m4d_em_a2`
  - `m4d_em_pi2`
  - `m4d_plasma_cons2`
  - `m4d_int_jw_ew`
  - `m4d_int_s_leak`
  - `m4d_pulse_xc`
  - `m4d_em_a2_mode_<n>`
  - `m4d_em_pi2_mode_<n>`

Not implemented yet:
- full two-fluid mode dynamics and `J^w` source evolution
- full conservative EM update in the AthenaPK flux pipeline (current path is source-step)
- full reconnection workflow/analysis

Additional bring-up path now available:
- `inputs/em4d_pulse.in` exercises EM-only source evolution for `em4d_a`/`em4d_pi`

## Files

- `mode_tables.hpp`, `mode_tables.cpp`: Hermite basis, quadrature, transforms, sparse operator helpers
- `mode_transform.hpp`, `mode_transform.cpp`: forward/inverse transform wrappers
- `package4d_modes.hpp`, `package4d_modes.cpp`: `modes4d` package setup and parameters
- `em4d_modes.*`: registered EM mode fields (`em4d_a`, `em4d_pi`)
- `plasma4d_modes.*`: registered plasma mode fields (`plasma4d_cons`)
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

## Notes

- The `harris_4d` initializer currently sets:
  - hydro/MHD conservative state (`cons`) for a Harris-like zero-mode profile
  - mode-0 `A_y` perturbation in `em4d_a`
  - placeholder species zero-mode values in `plasma4d_cons`
- This is a bring-up scaffold and not yet the final physics model from the docs.
