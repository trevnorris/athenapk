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
- optional conservative EM transport path in the hydro flux pipeline:
  - `modes4d/em_conservative_transport = true` registers `em4d_pi` with `WithFluxes`
  - stage fluxes add `-c_wave^2 * grad(em4d_a)` transport for `em4d_pi`
  - conservative spatial mixed-sector couplings are now included in flux form
    following `docs/mode_equations_verified.md`:
    - `pi_ax/pi_ay/pi_az` fluxes include
      `-(c_wave^2 * sqrt(2n)/lambda) * a_w^(n-1)`
    - `pi_aw` fluxes include
      `+(c_wave^2 * sqrt(2(n+1))/lambda) * a_{x,y,z}^(n+1)` contributions
  - source-step Laplacian terms are disabled when this path is enabled to avoid
    double counting
  - source-step spatial mixed couplings are disabled in this mode to avoid
    flux/source double counting; non-conservative mode still applies them as
    source terms
  - time-like mixed couplings from `docs/mode_equations_verified.md` are now
    explicitly included in source form for all transport modes:
    - `nu=0` brane component includes the
      `+sqrt(2n)/lambda * d^0 a_w^(n-1)` coupling term
    - scalar component includes the time-like piece of
      `-sqrt(2(n+1))/lambda * d_mu a^{mu,(n+1)}`
- split plasma transverse-momentum (`momw`) source update:
  - reconstruct `E_w` and `C_a` at quadrature nodes from EM modes
  - apply species force `q/m * rho * (E_w - v^a C_a)` in node space
  - apply transverse pressure-gradient source `-<phi_n, d_w p>_Z` via ID-3 mode coupling
    (runtime gain: `modes4d/plasma_momw_pressure_source_gain`)
  - project updated `momw` back to mode coefficients
- split plasma brane-momentum + energy source update:
  - apply brane Lorentz terms `q/m * rho * (E + v x B + v_w C)` for `momx/momy/momz`
  - apply species work `q/m * rho * (v·E + v_w E_w)` for mode energy
  - project updated node values back to mode coefficients
- split plasma mode-continuity coupling update:
  - source-step term applies leakage coupling
    `d_t rho^(n) = -(sqrt(2(n+1))/lambda) j_w^(n+1)`
  - optional source-step `w`-flux coupling for species momentum/energy modes:
    - applies projected `-<phi_n, d_w F_w(U)>_Z` leakage-style updates for
      `momx/momy/momz/momw/energy`
    - runtime gain: `modes4d/plasma_w_flux_source_gain`
  - brane transport `div(j^a,(n))` is handled in the conservative flux pipeline by
    registering `plasma4d_cons` with `WithFluxes` and filling stage fluxes for
    `rho/momx/momy/momz/momw/energy` each stage
  - runtime gain `plasma_rho_divj_gain` scales that conservative brane-transport path
  - write updated `plasma4d_cons` back each source step
- conservative plasma transport now supports a pseudospectral mode-coupling path:
  - reconstruct left/right face states at Gauss-Hermite nodes
  - compute node-local directional fluxes
  - project node fluxes back to modal flux coefficients
  - runtime toggle: `modes4d/plasma_pseudospectral_transport`
- `harris_4d` problem hook and zero-mode Harris initialization scaffold
- `harris_4d` now uses the same unsplit equation-based EM source update during runtime
- standalone unit-test executable for mode math (`modes4d_unit_tests`)
- initial history diagnostics:
  - `m4d_em_a2`
  - `m4d_em_pi2`
  - `m4d_plasma_cons2`
  - `m4d_int_jw_ew`
  - `m4d_int_ja_ea`
  - `m4d_int_s_leak`
  - `m4d_int_s_leak_abs`
  - `m4d_int_divj_mode0`
  - `m4d_int_divmomx_mode0`
  - `m4d_int_divmomy_mode0`
  - `m4d_int_divmomz_mode0`
  - `m4d_int_divmomw_mode0`
  - `m4d_int_divenergy_mode0`
  - `m4d_int_divpi0_mode0`
  - `m4d_int_divpix_mode0`
  - `m4d_int_divpiy_mode0`
  - `m4d_int_divpiz_mode0`
  - `m4d_int_divpiw_mode0`
  - `m4d_int_srcmomx_mode0`
  - `m4d_int_srcmomy_mode0`
  - `m4d_int_srcmomz_mode0`
  - `m4d_int_srcmomw_mode0`
  - `m4d_int_srcenergy_mode0`
  - `m4d_int_srcpi0_mode0`
  - `m4d_int_srcpix_mode0`
  - `m4d_int_srcpiy_mode0`
  - `m4d_int_srcpiz_mode0`
  - `m4d_int_srcpiw_mode0`
  - `m4d_int_src_timelike_a0_from_piw`
  - `m4d_int_src_timelike_a0_from_piw_abs`
  - `m4d_int_src_timelike_aw_from_pi0`
  - `m4d_int_src_timelike_aw_from_pi0_abs`
  - `m4d_brane_e2`
  - `m4d_brane_b2`
  - `m4d_brane_epar2`
  - `m4d_mixed_ew2`
  - `m4d_mixed_c2`
  - `m4d_em_u_bulk`
  - `m4d_em_u_resolved`
  - `m4d_helicity_sub`
  - `m4d_edotb_sub`
  - `m4d_em_sw`
  - `m4d_em_leak_w`
  - `m4d_psi0_span`
  - `m4d_pulse_xc`
  - `m4d_em_a2_mode_<n>`
  - `m4d_em_pi2_mode_<n>`
  - `m4d_pi0_mode_0`
  - `m4d_pix_mode_0`
  - `m4d_piy_mode_0`
  - `m4d_piz_mode_0`
  - `m4d_piw_mode_0`
  - `m4d_jw_mode_l2_<n>`
  - `m4d_jw_mode_<n>`
  - `m4d_charge_mode_<n>`
  - `m4d_momx_mode_<n>`
  - `m4d_momy_mode_<n>`
  - `m4d_momz_mode_<n>`
  - `m4d_momw_mode_<n>`
  - `m4d_energy_mode_<n>`
  - local continuity closure channels (source-step residuals):
    - `m4d_cont_local_l1`
    - `m4d_cont_local_l2`
    - `m4d_cont_local_max_abs`
    - `m4d_cont_mode0_l1`
    - `m4d_cont_mode0_l2`
    - `m4d_cont_mode0_max_abs`

Not implemented yet:
- full two-fluid closure beyond the current Lorentz + transverse pressure-gradient
  source path plus conservative transport bring-up
- pressure-coupled transport physics calibration beyond the current stabilized
  experimental path (the numerics are now bounded, but this path remains
  outside baseline CI)
- full conservative EM closure beyond the current spatial mixed-coupling transport path
  (`em4d_pi` transport now covers brane Laplacian + nearest-neighbor spatial
  mixed couplings, with time-like mixed terms wired in source form; mass/current/
  damping/gauge-control terms remain source-step)
- full reconnection workflow/analysis

Additional bring-up path now available:
- `inputs/em4d_pulse.in` exercises EM-only source evolution for `em4d_a`/`em4d_pi`
- `inputs/em4d_scalar_pulse.in` seeds a scalar-channel `A_w` pulse (`aw_mode1`)
  for mixed-sector activity checks
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
- `-DAthenaPK_ENABLE_4D_MODES_REGRESSION_TARGET=ON|OFF`

Defaults are `ON` for all three.

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
- `em_conservative_transport`
- `em_source_mass_gain`
- `em_source_current_gain`
- `em_source_damping_gain`
- `em_source_timelike_gain`
- `plasma_qom_ion`
- `plasma_qom_electron`
- `plasma_force_source_gain`
- `plasma_momw_source_gain`
- `plasma_momw_pressure_source_gain`
- `plasma_momw_damping`
- `plasma_w_flux_source_gain`
- `plasma_rho_divj_gain`
- `plasma_rho_floor`
- `plasma_energy_source_gain`
- `plasma_energy_floor`
- `plasma_gamma`
- `plasma_pressure_transport_gain`
- `plasma_pressure_rusanov_gain`
- `plasma_pressure_signal_speed_cap`
- `plasma_pressure_flux_relative_cap`
- `plasma_pressure_transport_max_mode`
- `plasma_pressure_floor`
- `plasma_pseudospectral_transport`

`<problem/em4d_pulse>` also supports:
- `component` (`ay_mode0`, `aw_mode0`, `a0_mode1`, `aw_mode1`)
- `pi_pulse_scale` (initial `pi` pulse multiplier)

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
The momentum/energy transport companions (`m4d_int_divmomx_mode0`,
`m4d_int_divmomy_mode0`, `m4d_int_divmomz_mode0`, `m4d_int_divmomw_mode0`,
`m4d_int_divenergy_mode0`) are accumulated the same way from corrected
`plasma4d_cons` fluxes.
Full-case activity gates can be enabled with:
- `--full-min-jw-ew-abs`
- `--full-min-s-leak-abs`
- `--full-min-mixed-ew2`
- `--full-min-jw-mode-l2-1`
- `--full-min-abs-em-leak-w`
- `--full-min-abs-helicity-sub`
- `--full-min-abs-edotb-sub`
- `--full-min-src-em-laplacian-abs`
- `--full-min-src-em-mass-abs`
- `--full-min-src-em-current-abs`
- `--full-min-src-em-spatial-mixed-abs`
- `--full-min-src-em-timelike-abs`
When set, these enforce nontrivial leakage/work/mixed-channel signal levels in the
full case and report `activity_status` / `activity_failures` in the CSV output.
The scan CSV also reports mode-0 transport-balance rates for momentum/energy:
- `momx_transport_*`
- `momy_transport_*`
- `momz_transport_*`
- `momw_transport_*`
- `energy_transport_*`
where each rate is computed as:
`d/dt(mode0_quantity) + d/dt(m4d_int_div*_mode0) - d/dt(m4d_int_src*_mode0)`.
The source-ledger channels (`m4d_int_srcmomx_mode0`, `m4d_int_srcmomy_mode0`,
`m4d_int_srcmomz_mode0`, `m4d_int_srcmomw_mode0`, `m4d_int_srcenergy_mode0`)
track source-step contributions, so the transport rate isolates unresolved residual.
The scan CSV also reports an EM bulk-ledger residual-rate metric:
`d/dt(m4d_em_u_bulk) + d/dt(m4d_int_ja_ea) + d/dt(m4d_int_jw_ew)`,
with channels:
- `em_bulk_ledger_max_abs_rate`
- `em_bulk_ledger_rms_abs_rate`
- `em_bulk_ledger_final_rate`
- `em_bulk_ledger_status`
and optional gating:
- `--check-em-bulk-ledger`
- `--em-bulk-ledger-abs-rate-tol`
EM mode-0 transport-balance rates are also reported:
- `pi0_transport_*`
- `pix_transport_*`
- `piy_transport_*`
- `piz_transport_*`
- `piw_transport_*`
using
`d/dt(m4d_pi*_mode_0) + d/dt(m4d_int_divpi*_mode0) - d/dt(m4d_int_srcpi*_mode0)`.
Use `--check-transport-closure --fail-on-check` to include these transport
status checks in nonzero-exit gating (`transport_closure_status` now covers both
plasma and EM channels).

### 6) Run controlled-limit regression smoke target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_controlled_regression
```

This target runs `scripts/modes4d_controlled_regression.py`, which:
- runs baseline AthenaPK smoke cases
  (`linear_wave3d`, `brio_wu`, `sod`, `orszag_tang`)
- applies short smoke-runtime overrides (`tlim=0.05`, `nlim=200`)
- auto-disables HDF5/phdf output blocks in copied inputs (`dt=-1`) for
  non-HDF5 builds
- runs `inputs/harris_4d_controlled.in`
- enforces controlled-limit inactivity gates on:
  - `m4d_int_jw_ew`
  - `m4d_int_s_leak_abs`
  - `m4d_mixed_ew2`
  - `m4d_mixed_c2`

### 7) Run controlled-limit parity regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_controlled_parity_regression
```

This target runs `scripts/modes4d_controlled_parity_regression.py`, which:
- runs `linear_wave3d`, `brio_wu`, `sod`, and `orszag_tang` in three variants:
  - baseline
  - forced `<modes4d> enabled=true, n_modes=1`
  - forced `<modes4d> enabled=true, n_modes=4`
- applies short smoke runtime overrides (`tlim=0.05`, `nlim=200`)
- compares final `time`, `mass`, `1-mom`, `2-mom`, `3-mom`, `tot-E` parity
  between baseline and `n_modes=1`
- for `linear_wave3d`, compares `linearwave-errors.dat` parity between baseline
  and `n_modes=1` (RMS/L1/max columns)
- for `n_modes=4`, enforces controlled-limit inactivity of:
  - all `m4d_*_mode_<n>` channels with `n>=1`
  - mixed/leak channels (`m4d_mixed_ew2`, `m4d_mixed_c2`,
    `m4d_int_jw_ew`, `m4d_int_s_leak`, `m4d_int_s_leak_abs`)
- fails if any of those parity/inactivity checks exceed configured tolerances

### 8) Run Gaussian KK Coulomb+Yukawa regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_kk_coulomb_regression
```

This target runs `scripts/kk_coulomb_yukawa_regression.py`, validating the
Gaussian tower formulas used in `docs/4d_em_fields_summary.md`:
- `m_n^2 = 2n/lambda^2`
- odd-mode decoupling (`c_{2m+1}=0`)
- even-mode couplings (`c_{2m}/c_0 = (1/4^m) * binom(2m,m)`)
- leading static correction
  `A_0(r) = A_0^(Coulomb) * [1 + 1/2 * exp(-2r/lambda) + ...]`

### 9) Run scalar-channel pulse regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_scalar_pulse_regression
```

This target runs `scripts/em4d_scalar_pulse_regression.py` against
`inputs/em4d_scalar_pulse.in` and enforces mixed-sector activity checks:
- `max(m4d_mixed_ew2)`
- `max(m4d_mixed_c2)`
- `max(m4d_em_a2_mode_1)`
- `|m4d_pulse_xc(final)-m4d_pulse_xc(initial)|`

### 10) Run conservative-EM transport regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_em_conservative_regression
```

This target runs `scripts/modes4d_em_conservative_regression.py` and checks
the optional conservative EM transport path
(`modes4d/em_conservative_transport=true`) with:
- scalar-channel pulse activity gates
- short Harris full-case smoke with finite-history checks and mixed-channel activity
- mode-0 EM transport-rate closure gates for
  `pi0/pix/piy/piz/piw`, using
  `d/dt(m4d_pi*_mode_0) + d/dt(m4d_int_divpi*_mode0) - d/dt(m4d_int_srcpi*_mode0)`
  with normalized and absolute-rate tolerances
- short-window EM bulk-ledger absolute-rate gate:
  `d/dt(m4d_em_u_bulk) + d/dt(m4d_int_ja_ea) + d/dt(m4d_int_jw_ew)`
  with `--em-bulk-ledger-abs-rate-tol` (default `5e-2`)

Companion targeted time-like coupling check:

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_em_timelike_regression
```

This target runs `scripts/modes4d_em_timelike_regression.py` and executes two
short pulse cases tuned to exercise the source-step time-like mixed couplings:
- `aw_mode0` seeded with `pi` checks the `pi_w^(n-1) -> rhs_a0^(n)` path
- `a0_mode1` seeded with `pi` checks the `pi_0^(n+1) -> rhs_aw^(n)` path

Companion conservative-path diagnostics check:

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_em_conservative_timelike_diag_regression
```

This target runs `scripts/modes4d_em_conservative_timelike_diag_regression.py`
and verifies conservative-path history channels for explicit time-like source
terms are present, finite, and active where expected:
- `m4d_int_src_timelike_a0_from_piw`
- `m4d_int_src_timelike_a0_from_piw_abs`
- `m4d_int_src_timelike_aw_from_pi0`
- `m4d_int_src_timelike_aw_from_pi0_abs`

Default gates now enforce:
- `aw_mode0` inactivity ceilings on conservative-path time-like diagnostics
- `a0_mode1` activation floors for both absolute diagnostic deltas

### 11) Run the dedicated tuned Harris regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_regression
```

This target runs `scripts/harris_regression_tuned.py`, which executes the
controlled/full scan with fixed closure + activity thresholds for the tuned
`inputs/harris_4d_full.in` baseline, and includes transport-closure gating
(`--check-transport-closure`, plasma + EM mode-0 channels), EM bulk-ledger
gating (`--check-em-bulk-ledger`) plus full-case correlation gates tying
`m4d_psi0_span` to leakage/mixed-channel activity:
- `|corr_psi0_s_leak_abs|`
- `|corr_psi0_jw_ew|`
- `|corr_psi0_mixed_ew2|`
- `|corr_psi0_mixed_c2|`
- `|corr_psi0_jw_mode_l2_1|`
- `|corr_psi0_em_leak_w|`
- `|corr_psi0_helicity_sub|`
- `|corr_psi0_edotb_sub|`

It also enforces full-case activity floors on:
- `final_em_a2_mode_1`
- `final_em_pi2_mode_1`
- `|final_em_leak_w|`
- `|final_helicity_sub|`
- `|final_edotb_sub|`
- `|final_int_src_em_laplacian_abs|`
- `|final_int_src_em_mass_abs|`
- `|final_int_src_em_current_abs|`
- `|final_int_src_em_spatial_mixed_abs|`
- `|final_int_src_em_timelike_abs|`

### 12) Run the tuned Harris pseudospectral-transport regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_pseudospectral_regression
```

This target runs the tuned Harris gate set with explicit transport-mode
overrides so the coupling mode is always validated regardless of deck defaults:
- controlled override:
  - `modes4d/plasma_pseudospectral_transport=false`
- full override:
  - `modes4d/plasma_pseudospectral_transport=true`

It enforces the same closure/transport/activity/correlation gates as
`modes4d_harris_regression`.

### 13) Run the tuned Harris `N_w`-convergence regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_nw_convergence_regression
```

This target runs `scripts/harris_nw_convergence_regression.py` and enforces:
- per-run closure gates across `N_w = 2,3,4`:
  - `closure_status = PASS`
  - `closure_local_mode0_status = PASS`
  - `transport_closure_status = PASS`
- highest-mode activity floors at `N_w=4`:
  - `|final_jw_ew| >= 1e-15`
  - `final_s_leak_abs >= 1e-14`
  - `final_mixed_ew2 >= 1e-3`
  - `final_em_a2_mode_1 >= 1e-6`
  - `final_em_pi2_mode_1 >= 1e-3`
  - `|final_em_leak_w| >= 1e-4`
  - `|final_helicity_sub| >= 1e-12`
  - `|final_edotb_sub| >= 1e-11`
- high-mode stabilization checks between the two highest `N_w` runs (`3 -> 4`)
  using relative-delta gating (`<= 0.7`) for:
  - `final_psi0_span`
  - `final_s_leak_abs`
  - `final_jw_ew`
  - `final_mixed_ew2`
  - `final_mixed_c2`
  - `final_em_a2_mode_1`
  - `final_em_pi2_mode_1`
  - `final_em_leak_w`
  - `final_helicity_sub`
  - `final_edotb_sub`

### 14) Run the experimental pressure-transport Harris regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_pressure_experimental_regression
```

This target runs the same tuned Harris gate set as `modes4d_harris_regression`,
but adds a full-case override:
- `modes4d/plasma_pressure_transport_gain=1.0e-3`

It is an experimental check for the stabilized nonzero pressure-transport path
and is not part of `modes4d_phase10_regression`.

### 15) Run the pressure-window Harris regression target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_pressure_window_regression
```

This target scans finite nonzero pressure gains using tuned Harris gates:
- `1e-8`
- `1e-7`
- `1e-6`
- `1e-3`

It enforces:
- tuned closure/correlation/activity gates (same as `modes4d_harris_regression`)
- transport-closure gating
- explicit non-finite (`nan`/`inf`) checks on full-case ledger channels

It is a dedicated pressure-path regression and is intentionally not part of
`modes4d_phase10_regression`.

### 16) Run the Phase-10 validation bundle target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_phase10_regression
```

This target runs the full non-MPI Phase-10 validation stack in one command,
with fail-fast behavior on the first failing check:
- `modes4d_controlled_regression`
- `modes4d_controlled_parity_regression`
- `modes4d_cpaw_strict1d_regression`
- `modes4d_kk_coulomb_regression`
- `modes4d_scalar_pulse_regression`
- `modes4d_em_timelike_regression`
- `modes4d_em_conservative_timelike_diag_regression`
- `modes4d_em_conservative_regression`
- `modes4d_em_conservative_transport_ledger_regression`
- `modes4d_em_conservative_source_channels_regression`
- `modes4d_harris_regression`
- `modes4d_harris_pseudospectral_regression`
- `modes4d_harris_nw_convergence_regression`
- `modes4d_harris_lundquist_scan`

### 17) Run the Harris Lundquist scan target

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_lundquist_scan
```

This target runs `scripts/harris_lundquist_scan.py`, which performs a controlled
vs full Harris scan over multiple Lundquist numbers and writes:
- per-`S` run logs under
  `build-baseline/modes4d_harris_lundquist_scan/outputs/S_<value>/`
- aggregate CSV summary:
  `build-baseline/modes4d_harris_lundquist_scan/outputs/lundquist_scan_summary.csv`

Default scan settings:
- `S = 250,500,1000,2000`
- `eta = (L * V_A) / S`
- inferred from `inputs/harris_4d_full.in` unless overridden:
  - `L = x1max - x1min`
  - `V_A = |b0| / sqrt(mu0 * n_bg)`

Each scan point forwards resistive runtime overrides to AthenaPK:
- `diffusion/resistivity=ohmic`
- `diffusion/resistivity_coeff=fixed`
- `diffusion/integrator=rkl2`
- `diffusion/rkl2_max_dt_ratio=100.0`
- `diffusion/ohm_diff_coeff_code=<eta>`

It reports per-`S` closure/transport/correlation/activity statuses from
`harris_scan_matrix.py` and scan-level correlations versus `log10(S)`.

This CMake target enforces scan-level acceptance gates and exits nonzero on
failure:
- minimum absolute correlations:
  - `|corr_logS_full_final_psi0_span| >= 0.7`
  - `|corr_logS_full_final_s_leak_abs| >= 0.5`
  - `|corr_logS_full_final_jw_ew_abs| >= 0.5`
  - `|corr_logS_full_final_mixed_ew2| >= 0.7`
  - `|corr_logS_full_final_em_leak_w_abs| >= 0.2`
  - `|corr_logS_full_final_helicity_sub_abs| >= 0.5`
  - `|corr_logS_full_final_edotb_sub_abs| >= 0.5`
- signed trend checks:
  - `corr_logS_full_final_helicity_sub_abs <= -0.2`
  - `corr_logS_full_final_edotb_sub_abs <= -0.2`

The scan script can also be used directly with custom gates:
- `--scan-min-abs-corr metric=threshold`
- `--scan-min-corr metric=threshold`
- `--scan-max-corr metric=threshold`

### 18) Run MPI + HDF5 Harris full-case helper

```bash
./scripts/run_harris_full_mpi_hdf5.sh --ranks 2 --tlim 0.02 --nlim 200
```

This helper runs `inputs/harris_4d_full.in` with MPI and forces an HDF5 output
cadence via `parthenon/output0/dt`, writing outputs under:
- `build-mpi-hdf5/runs/harris_full_mpi/`

### 19) Generate a `phdf` quicklook slice

```bash
python3 scripts/plot_phdf_slice.py \
  --phdf build-mpi-hdf5/runs/harris_full_mpi/parthenon.out0.final.phdf \
  --dataset em4d_a \
  --component 0 \
  --axis z \
  --output build-mpi-hdf5/runs/harris_full_mpi/quicklook_em4d_a0_z.png
```

This produces a 2D PNG slice from a cell-centered Parthenon dataset in `phdf`
output (single-level runs).

### 20) Run projected two-fluid `w`-flux regression

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_harris_wflux_regression
```

This target runs `scripts/harris_wflux_regression.py`, which executes two tuned
Harris scan runs:
- baseline full-case `modes4d/plasma_w_flux_source_gain=0.0`
- enabled full-case `modes4d/plasma_w_flux_source_gain=0.1`

and enforces quantitative nonzero-effect gates on `final_jw_mode_1` while
requiring closure/correlation/activity/transport statuses to remain `PASS`.

### 21) Run strict conservative-EM transport/source ledger gate

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_em_conservative_transport_ledger_regression
```

This target runs `scripts/modes4d_em_conservative_regression.py` with stricter
EM transport closure settings:
- longer short-Harris window (`tlim=0.05`, `nlim=500`)
- tighter transport thresholds
  (`transport_norm_tol=1e-2`, `transport_abs_rate_tol=1e-10`)

It keeps scalar/mixed-channel conservative checks and enforces mode-0
`em4d_pi` transport/source ledger closure in conservative mode.

### 22) Run conservative EM source-channel decomposition/gain gate

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_em_conservative_source_channels_regression
```

This target runs `scripts/modes4d_em_conservative_source_channels_regression.py`
and checks the term-resolved conservative EM source diagnostics:
- `m4d_int_src_em_laplacian_abs`
- `m4d_int_src_em_mass_abs`
- `m4d_int_src_em_current_abs`
- `m4d_int_src_em_damping_abs`
- `m4d_int_src_em_spatial_mixed_abs`
- `m4d_int_src_em_timelike_abs`

It enforces, in short Harris runs:
- conservative path suppresses Laplacian and spatial-mixed source channels
- non-conservative path activates Laplacian and spatial-mixed source channels
- mass/current/timelike channels remain active in conservative baseline
- `em_source_mass_gain=0`, `em_source_current_gain=0`,
  `em_source_timelike_gain=0` each suppress their corresponding channel
- damping channel is exercised with `em_damping>0` and then suppressed by
  `em_source_damping_gain=0`

### 23) Run strict-1D CPAW crash-regression gate

```bash
cmake --build /projects/fluid-engine/athenapk/build-baseline \
  --target modes4d_cpaw_strict1d_regression
```

This target runs `scripts/cpaw_strict1d_regression.py`, which:
- executes `inputs/cpaw.in` with strict-1D overrides
  (`nx2=nx3=1`, with matching meshblock overrides)
- disables HDF5 output in non-HDF5 builds (`parthenon/output0/dt=-1`)
- applies short crash-gating runtime limits (`tlim=0.01`, `nlim=50`)
- verifies run completion:
  - zero exit status
  - `"Driver completed."` in run logs
  - at least one completed cycle (`final_cycle >= 1`)
- verifies quantitative CPAW parity diagnostics from `cpaw-errors.dat`:
  - file exists and reflects strict-1D mesh overrides (`Nx2=1`, `Nx3=1`)
  - final `RMS-Error` is finite and below a configured ceiling
    (`--max-rms-error`, default `1e-2`)
  - final transverse component errors are active
    (`|B2c|`, `|B3c| >= --min-transverse-error`, default `1e-6`)
  - final transverse component errors remain symmetric
    (`|B2c-B3c| / max(|B2c|,|B3c|) <= --max-transverse-rel-diff`,
    default `1e-2`)

The same gate is included in `modes4d_phase10_regression`.

### 24) Run strict-1D CPAW MPI+HDF5 quantitative gate

```bash
cmake --build /projects/fluid-engine/athenapk/build-mpi-hdf5 \
  --target modes4d_cpaw_strict1d_mpi_hdf5_regression
```

This target is available only in MPI+HDF5-enabled builds
(`PARTHENON_DISABLE_MPI=OFF`, `PARTHENON_DISABLE_HDF5=OFF`).

It runs `scripts/cpaw_strict1d_mpi_hdf5_regression.py`, which:
- launches strict-1D CPAW with MPI (`mpirun -n 2`) and HDF5 outputs enabled
- enforces strict-1D mesh overrides (`nx2=1`, `nx3=1`, with meshblock matches)
- checks completion and output presence:
  - zero exit status
  - `"Driver completed."` in logs
  - final cycle floor (`>= 1`)
  - final `parthenon.out0*.phdf` exists
- checks quantitative field metrics from final `cons` in `.phdf`:
  - all values finite
  - strict-1D geometry reflected in output layout
  - transverse magnetic activity floors:
    - `mean(|B2|) >= 1e-4`
    - `mean(|B3|) >= 1e-4`
  - transverse parity:
  - `|mean(|B2|)-mean(|B3|)| / max(mean(|B2|),mean(|B3|)) <= 1e-2`

### 25) Run MPI+HDF5 regression bundle

```bash
cmake --build /projects/fluid-engine/athenapk/build-mpi-hdf5 \
  --target modes4d_mpi_hdf5_regression
```

This target is available only in MPI+HDF5-enabled builds and runs, in order:
- `modes4d_cpaw_strict1d_mpi_hdf5_regression`
- `scripts/run_harris_full_mpi_hdf5.sh` with:
  - `--ranks 2`
  - `--tlim 0.02`
  - `--nlim 200`
  - `--output-dt 0.01`

It provides a one-command MPI smoke bundle covering:
- strict-1D CPAW quantitative parity in `.phdf`
- full-channel Harris runtime smoke with MPI/HDF5 outputs

CI integration:
- `.github/workflows/ci.yml` runs this bundle in the `mpi` matrix leg on
  pull requests before merge.

Current pressure-transport status for tuned Harris decks:
- default `plasma_pressure_transport_gain = 0.0` is regression-stable
- nonzero pressure transport now uses a bounded blended Rusanov path:
  - `plasma_pressure_transport_gain` is a blend factor in `[0,1]`
    between the baseline upwind transport (`0`) and pressure-coupled Rusanov
    transport (`1`)
  - `plasma_pressure_rusanov_gain` scales the Rusanov signal speed
  - `plasma_pressure_signal_speed_cap` clamps pressure-branch
    velocity/signal-speed contributions (default `10.0`)
  - `plasma_pressure_flux_relative_cap` limits pressure-branch fluxes
    relative to baseline advective fluxes (default `50.0`)
  - `plasma_pressure_transport_max_mode` limits pressure-coupled transport to
    low modes (default `0`, i.e. mode-0 only)
- tested full-case behavior:
  - stable through `1e-3` (`1e-8`, `1e-7`, `1e-6`, `1e-3` verified)
  - pressure-on/off Lundquist comparison (`gain=0` vs `gain=1e-3`) under
    the current scan gates shows no observable trend shift at printed precision
- calibration sweeps at `gain=1e-3` and stress checks at `gain=1e-1`:
    - `plasma_pressure_transport_max_mode = 0..3`: PASS
    - `plasma_pressure_rusanov_gain = 0.25..4.0`: PASS
    - `plasma_pressure_flux_relative_cap = 10..200`: PASS
    - `plasma_pressure_signal_speed_cap = 5..20`: PASS
- recommended pressure-stabilizer defaults (current package defaults):
  - `plasma_pressure_transport_max_mode = 0`
  - `plasma_pressure_rusanov_gain = 1.0`
  - `plasma_pressure_signal_speed_cap = 10.0`
  - `plasma_pressure_flux_relative_cap = 50.0`
- keep `plasma_pressure_transport_gain = 0.0` for baseline/CI unless explicitly
  testing the pressure-transport path

Current nonlinear plasma-transport status for tuned Harris decks:
- `modes4d/plasma_pseudospectral_transport` enables the docs-specified
  `reconstruct -> multiply -> project` path for conservative plasma fluxes.
- fallback mode-diagonal face fluxes remain available with
  `modes4d/plasma_pseudospectral_transport=false`.
- tuned full-case gate check with pseudospectral transport enabled:
  - `python3 scripts/harris_scan_matrix.py ... --fail-on-check --check-transport-closure --full-arg modes4d/plasma_pseudospectral_transport=true`: PASS
  - closure/transport/correlation/activity statuses remained `PASS`.
- dedicated CMake gate:
  - `cmake --build build-baseline --target modes4d_harris_pseudospectral_regression`: PASS
- Phase-10 bundle now includes that explicit pseudospectral gate run.

Current transverse two-fluid source status for tuned Harris decks:
- `modes4d/plasma_momw_pressure_source_gain` controls the ID-3-projected
  transverse pressure-gradient source in species `momw` updates:
  - source term: `-<phi_n, d_w p_s>_Z`
  - evaluated from node-space pressure reconstruction and projected through
    `ApplyID3Raising` (no nodal finite-difference in `w`).
- tuned Harris decks currently set:
  - `plasma_momw_pressure_source_gain = 1.0`
- `modes4d/plasma_w_flux_source_gain` controls optional projected `w`-flux
  source coupling for species momentum/energy modes:
  - applies leakage-style projected couplings
    `-<phi_n, d_w F_w(U)>_Z` to
    `momx/momy/momz/momw/energy`
  - `F_w` is reconstructed in node space from modal state:
    - `F_w(momx)=momx*v_w`
    - `F_w(momy)=momy*v_w`
    - `F_w(momz)=momz*v_w`
    - `F_w(momw)=momw*v_w + p`
    - `F_w(energy)=(energy + p)*v_w`
  - default is `0.0` (baseline unchanged), with focused smoke verified at
    `plasma_w_flux_source_gain=0.1`.
- dedicated CMake gate:
  - `cmake --build build-baseline --target modes4d_harris_wflux_regression`: PASS
- dedicated CMake gate:
  - `cmake --build build-baseline --target modes4d_em_conservative_transport_ledger_regression`: PASS
- dedicated CMake gate:
  - `cmake --build build-baseline --target modes4d_em_conservative_source_channels_regression`: PASS
- tuned Harris regression (`modes4d_harris_regression`): PASS

Current `N_w`-convergence status for tuned Harris decks:
- dedicated CMake gate:
  - `cmake --build build-baseline --target modes4d_harris_nw_convergence_regression`: PASS
- closure and transport checks pass for `N_w = 2,3,4`, and highest-mode
  activity floors plus `N_w=3 -> 4` relative-delta stabilization gates pass.
- Phase-10 bundle now includes this `N_w` convergence gate run.

## Dependency setup

From the AthenaPK repository root, install system + Python dependencies with:

```bash
./scripts/install_deps.sh --yes
```

Python-only dependencies are tracked in `requirements.txt` and can be installed
with:

```bash
./scripts/install_deps.sh --python-only
```

If MPI detection fails even with OpenMPI installed, configure AthenaPK with
wrapper compilers explicitly:

```bash
cmake -S . -B build-mpi -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
```

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
- `init_charge_rel_tol`
- `init_charge_abs_tol`

## Notes

- The `harris_4d` initializer currently sets:
  - hydro/MHD conservative state (`cons`) for a Harris-like zero-mode profile
  - mode-0 `A_y` perturbation in `em4d_a`
  - optional mode-1 `A_w` perturbation in `em4d_a` for mixed-sector triggering
  - species mode-0 plasma state in `plasma4d_cons` using q/m-aware splits that
    preserve mode-0 charge neutrality and target drift-current amplitude
- The provided pulse/Harris decks now carry explicit defaults for frequently
  overridden regression knobs:
  - `modes4d/em_conservative_transport = false`
  - `problem/em4d_pulse/component`
  - `problem/em4d_pulse/pi_pulse_scale`
- This is a bring-up scaffold and not yet the final physics model from the docs.
