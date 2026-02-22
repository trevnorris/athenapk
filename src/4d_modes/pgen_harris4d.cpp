#include "pgen_harris4d.hpp"

#include <algorithm>
#include <cmath>

#include "../main.hpp"
#include "utils/error_checking.hpp"

namespace Modes4D {
using namespace parthenon::package::prelude;

void InitializeHarrisModes(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {
  const Real b0 = pin->GetOrAddReal("problem/harris_4d", "b0", 1.0);
  const Real guide_bz = pin->GetOrAddReal("problem/harris_4d", "guide_bz", 0.0);
  const Real n_bg = pin->GetOrAddReal("problem/harris_4d", "n_bg", 0.2);
  const Real n_sheet = pin->GetOrAddReal("problem/harris_4d", "n_sheet", 1.0);
  const Real p_bg = pin->GetOrAddReal("problem/harris_4d", "p_bg", 0.2);
  const Real sheet_half_width =
      pin->GetOrAddReal("problem/harris_4d", "sheet_half_width", 0.1);
  const Real sheet_w_offset =
      pin->GetOrAddReal("problem/harris_4d", "sheet_w_offset", 0.0);
  const Real perturbation_amp =
      pin->GetOrAddReal("problem/harris_4d", "perturbation_amp", 1.0e-3);
  const Real perturb_phase_x =
      pin->GetOrAddReal("problem/harris_4d", "perturb_phase_x", 0.0);
  const Real perturb_phase_z =
      pin->GetOrAddReal("problem/harris_4d", "perturb_phase_z", 0.0);
  const Real drift_current_scale =
      pin->GetOrAddReal("problem/harris_4d", "drift_current_scale", 0.0);
  const Real drift_momw_scale =
      pin->GetOrAddReal("problem/harris_4d", "drift_momw_scale", 0.0);
  const Real aw_mode1_amp = pin->GetOrAddReal("problem/harris_4d", "aw_mode1_amp", 0.0);
  const Real piw_mode1_amp = pin->GetOrAddReal("problem/harris_4d", "piw_mode1_amp", 0.0);
  const Real ay_mode1_amp = pin->GetOrAddReal("problem/harris_4d", "ay_mode1_amp", 0.0);
  const Real piy_mode1_amp = pin->GetOrAddReal("problem/harris_4d", "piy_mode1_amp", 0.0);
  const Real aw_mode1_phase_x =
      pin->GetOrAddReal("problem/harris_4d", "aw_mode1_phase_x", 0.0);
  const Real aw_mode1_phase_z =
      pin->GetOrAddReal("problem/harris_4d", "aw_mode1_phase_z", 0.0);
  const Real piw_mode1_phase_x =
      pin->GetOrAddReal("problem/harris_4d", "piw_mode1_phase_x", 0.0);
  const Real piw_mode1_phase_z =
      pin->GetOrAddReal("problem/harris_4d", "piw_mode1_phase_z", 0.0);
  const Real ay_mode1_phase_x =
      pin->GetOrAddReal("problem/harris_4d", "ay_mode1_phase_x", 0.0);
  const Real ay_mode1_phase_z =
      pin->GetOrAddReal("problem/harris_4d", "ay_mode1_phase_z", 0.0);
  const Real piy_mode1_phase_x =
      pin->GetOrAddReal("problem/harris_4d", "piy_mode1_phase_x", 0.0);
  const Real piy_mode1_phase_z =
      pin->GetOrAddReal("problem/harris_4d", "piy_mode1_phase_z", 0.0);
  const Real ay_mode2_amp = pin->GetOrAddReal("problem/harris_4d", "ay_mode2_amp", 0.0);
  const Real a0_mode2_amp = pin->GetOrAddReal("problem/harris_4d", "a0_mode2_amp", 0.0);
  const Real aw_mode2_amp = pin->GetOrAddReal("problem/harris_4d", "aw_mode2_amp", 0.0);
  const Real piw_mode2_amp = pin->GetOrAddReal("problem/harris_4d", "piw_mode2_amp", 0.0);
  const Real init_charge_rel_tol =
      pin->GetOrAddReal("problem/harris_4d", "init_charge_rel_tol", 1.0e-12);
  const Real init_charge_abs_tol =
      pin->GetOrAddReal("problem/harris_4d", "init_charge_abs_tol", 1.0e-12);
  const Real gamma = pin->GetOrAddReal("hydro", "gamma", 5.0 / 3.0);

  PARTHENON_REQUIRE(sheet_half_width > 0.0,
                    "problem/harris_4d/sheet_half_width must be > 0");
  PARTHENON_REQUIRE(gamma > 1.0, "hydro/gamma must be > 1");

  auto modes_pkg = pmb->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real lambda = modes_pkg->Param<double>("lambda");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");

  PARTHENON_REQUIRE(lambda > 0.0, "modes4d/lambda must be > 0");
  PARTHENON_REQUIRE(
      (std::abs(sheet_w_offset) == 0.0) || (n_modes > 1),
      "problem/harris_4d/sheet_w_offset requires modes4d/n_modes > 1");

  const Real x1min = pin->GetReal("parthenon/mesh", "x1min");
  const Real x1max = pin->GetReal("parthenon/mesh", "x1max");
  const Real x3min = pin->GetReal("parthenon/mesh", "x3min");
  const Real x3max = pin->GetReal("parthenon/mesh", "x3max");
  const Real lx = x1max - x1min;
  const Real lz = x3max - x3min;
  const Real gm1 = gamma - 1.0;
  const Real qom_denom = qom_ion - qom_electron;
  const Real sheet_mode1_coeff =
      (std::abs(sheet_w_offset) > 0.0) ? (sheet_w_offset / (std::sqrt(2.0) * lambda)) : 0.0;
  PARTHENON_REQUIRE(std::abs(qom_denom) > 0.0,
                    "modes4d plasma_qom_ion and plasma_qom_electron must differ");
  PARTHENON_REQUIRE((qom_ion * qom_electron) < 0.0,
                    "Harris modes initializer expects opposite-sign ion/electron q/m");
  constexpr Real kTwoPi = 6.2831853071795864769;

  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto &rc = pmb->meshblock_data.Get();
  auto &cons_dev = rc->Get("cons").data;
  auto &em_a_dev = rc->Get("em4d_a").data;
  auto &em_pi_dev = rc->Get("em4d_pi").data;
  auto &plasma_dev = rc->Get("plasma4d_cons").data;

  auto cons = cons_dev.GetHostMirrorAndCopy();
  auto em_a = em_a_dev.GetHostMirrorAndCopy();
  auto em_pi = em_pi_dev.GetHostMirrorAndCopy();
  auto plasma = plasma_dev.GetHostMirrorAndCopy();
  auto &coords = pmb->coords;

  const int n_em_vars = em_a.GetDim(4);
  const int n_plasma_vars = plasma.GetDim(4);
  Real total_charge = 0.0;
  Real charge_scale = 0.0;
  Real max_abs_charge = 0.0;

  for (int k = kb.s; k <= kb.e; ++k) {
    for (int j = jb.s; j <= jb.e; ++j) {
      for (int i = ib.s; i <= ib.e; ++i) {
        const Real x = coords.Xc<1>(i);
        const Real y = coords.Xc<2>(j);
        const Real z = coords.Xc<3>(k);
        const Real yhat = y / sheet_half_width;
        const Real tanhy = std::tanh(yhat);
        const Real sech2y = 1.0 / (std::cosh(yhat) * std::cosh(yhat));

        const Real rho = n_bg + n_sheet * sech2y;
        const Real pressure = p_bg + 0.5 * b0 * b0 * sech2y;
        const Real bx = b0 * tanhy;
        const Real by = 0.0;
        const Real bz = guide_bz;
        const Real magnetic_energy = 0.5 * ((bx * bx) + (by * by) + (bz * bz));

        cons(IDN, k, j, i) = rho;
        cons(IM1, k, j, i) = 0.0;
        cons(IM2, k, j, i) = 0.0;
        cons(IM3, k, j, i) = 0.0;
        cons(IEN, k, j, i) = (pressure / gm1) + magnetic_energy;
        cons(IB1, k, j, i) = bx;
        cons(IB2, k, j, i) = by;
        cons(IB3, k, j, i) = bz;
        cons(IPS, k, j, i) = 0.0;

        for (int n = 0; n < n_em_vars; ++n) {
          em_a(n, k, j, i) = 0.0;
          em_pi(n, k, j, i) = 0.0;
        }
        for (int n = 0; n < n_plasma_vars; ++n) {
          plasma(n, k, j, i) = 0.0;
        }

        const Real phase_x = (lx > 0.0) ? (kTwoPi * (x - x1min) / lx) + perturb_phase_x : 0.0;
        const Real phase_z = (lz > 0.0) ? (kTwoPi * (z - x3min) / lz) + perturb_phase_z : 0.0;
        const Real ay_perturb = perturbation_amp * std::cos(phase_x) * std::cos(phase_z);
        em_a(2, k, j, i) = ay_perturb; // mode 0, A_y component

        const Real phase = std::cos(phase_x) * std::cos(phase_z);
        if (n_modes > 1) {
          const int aw_mode1_idx = 5 + 4;
          const int ay_mode1_idx = 5 + 2;
          const Real phase_aw =
              std::cos(phase_x + aw_mode1_phase_x) * std::cos(phase_z + aw_mode1_phase_z);
          const Real phase_piw =
              std::cos(phase_x + piw_mode1_phase_x) * std::cos(phase_z + piw_mode1_phase_z);
          const Real phase_ay =
              std::cos(phase_x + ay_mode1_phase_x) * std::cos(phase_z + ay_mode1_phase_z);
          const Real phase_piy =
              std::cos(phase_x + piy_mode1_phase_x) * std::cos(phase_z + piy_mode1_phase_z);
          em_a(aw_mode1_idx, k, j, i) = aw_mode1_amp * phase_aw;
          em_pi(aw_mode1_idx, k, j, i) = piw_mode1_amp * phase_piw;
          em_a(ay_mode1_idx, k, j, i) = ay_mode1_amp * phase_ay;
          em_pi(ay_mode1_idx, k, j, i) = piy_mode1_amp * phase_piy;
          if (sheet_mode1_coeff != 0.0) {
            // Small physical w-offset: phi_0(w-dw) ~= phi_0(w) + dw/(sqrt(2)*lambda) * phi_1(w).
            em_a(ay_mode1_idx, k, j, i) += sheet_mode1_coeff * em_a(2, k, j, i);
          }
        }
        if (n_modes > 2) {
          const int mode2_off = 10;
          em_a(mode2_off + 0, k, j, i) = a0_mode2_amp * phase;
          em_a(mode2_off + 2, k, j, i) = ay_mode2_amp * phase;
          em_a(mode2_off + 4, k, j, i) = aw_mode2_amp * phase;
          em_pi(mode2_off + 4, k, j, i) = piw_mode2_amp * phase;
        }

        // Two-species (ions/electrons) mode state. Initialize only the zero mode.
        const Real jz_sheet = drift_current_scale * (b0 / sheet_half_width) * sech2y;
        const Real jw_sheet = drift_momw_scale * (b0 / sheet_half_width) * sech2y;
        const Real rho_ion = (-qom_electron / qom_denom) * rho;
        const Real rho_electron = (qom_ion / qom_denom) * rho;
        const Real momz_ion = jz_sheet / qom_denom;
        const Real momz_electron = -momz_ion;
        const Real momw_ion = jw_sheet / qom_denom;
        const Real momw_electron = -momw_ion;
        const Real energy_density = pressure / gm1;
        const Real energy_ion_internal = (-qom_electron / qom_denom) * energy_density;
        const Real energy_electron_internal = (qom_ion / qom_denom) * energy_density;
        const Real energy_ion = energy_ion_internal +
                                (0.5 * ((momz_ion * momz_ion) + (momw_ion * momw_ion)) / rho_ion);
        const Real energy_electron = energy_electron_internal +
                                     (0.5 * ((momz_electron * momz_electron) +
                                             (momw_electron * momw_electron)) /
                                      rho_electron);
        const Real charge_density = (qom_ion * rho_ion) + (qom_electron * rho_electron);
        const Real cell_volume = coords.CellVolume(k, j, i);
        total_charge += charge_density * cell_volume;
        charge_scale +=
            (std::abs(qom_ion * rho_ion) + std::abs(qom_electron * rho_electron)) * cell_volume;
        max_abs_charge = std::max(max_abs_charge, std::abs(charge_density));

        const int ion_offset = 0;
        plasma(ion_offset + 0, k, j, i) = rho_ion;
        plasma(ion_offset + 1, k, j, i) = 0.0;
        plasma(ion_offset + 2, k, j, i) = 0.0;
        plasma(ion_offset + 3, k, j, i) = momz_ion;
        plasma(ion_offset + 4, k, j, i) = momw_ion;
        plasma(ion_offset + 5, k, j, i) = energy_ion;

        const int ele_offset = 6 * n_modes;
        plasma(ele_offset + 0, k, j, i) = rho_electron;
        plasma(ele_offset + 1, k, j, i) = 0.0;
        plasma(ele_offset + 2, k, j, i) = 0.0;
        plasma(ele_offset + 3, k, j, i) = momz_electron;
        plasma(ele_offset + 4, k, j, i) = momw_electron;
        plasma(ele_offset + 5, k, j, i) = energy_electron;

        if (n_modes > 1 && sheet_mode1_coeff != 0.0) {
          const int ion_mode1_offset = 6;
          plasma(ion_mode1_offset + 0, k, j, i) = sheet_mode1_coeff * rho_ion;
          plasma(ion_mode1_offset + 1, k, j, i) = 0.0;
          plasma(ion_mode1_offset + 2, k, j, i) = 0.0;
          plasma(ion_mode1_offset + 3, k, j, i) = sheet_mode1_coeff * momz_ion;
          plasma(ion_mode1_offset + 4, k, j, i) = sheet_mode1_coeff * momw_ion;
          plasma(ion_mode1_offset + 5, k, j, i) = sheet_mode1_coeff * energy_ion;

          const int ele_mode1_offset = (6 * n_modes) + 6;
          plasma(ele_mode1_offset + 0, k, j, i) = sheet_mode1_coeff * rho_electron;
          plasma(ele_mode1_offset + 1, k, j, i) = 0.0;
          plasma(ele_mode1_offset + 2, k, j, i) = 0.0;
          plasma(ele_mode1_offset + 3, k, j, i) = sheet_mode1_coeff * momz_electron;
          plasma(ele_mode1_offset + 4, k, j, i) = sheet_mode1_coeff * momw_electron;
          plasma(ele_mode1_offset + 5, k, j, i) = sheet_mode1_coeff * energy_electron;
        }
      }
    }
  }

  const Real charge_norm = std::abs(total_charge) / std::max(charge_scale, 1.0e-30);
  PARTHENON_REQUIRE(
      charge_norm <= init_charge_rel_tol,
      "harris_4d initialization failed charge-neutrality relative check; "
      "tune problem/harris_4d/init_charge_rel_tol if needed");
  PARTHENON_REQUIRE(
      max_abs_charge <= init_charge_abs_tol,
      "harris_4d initialization failed charge-neutrality local absolute check; "
      "tune problem/harris_4d/init_charge_abs_tol if needed");

  cons_dev.DeepCopy(cons);
  em_a_dev.DeepCopy(em_a);
  em_pi_dev.DeepCopy(em_pi);
  plasma_dev.DeepCopy(plasma);
}

} // namespace Modes4D
