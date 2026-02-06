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
  const Real perturbation_amp =
      pin->GetOrAddReal("problem/harris_4d", "perturbation_amp", 1.0e-3);
  const Real drift_current_scale =
      pin->GetOrAddReal("problem/harris_4d", "drift_current_scale", 0.0);
  const Real aw_mode1_amp = pin->GetOrAddReal("problem/harris_4d", "aw_mode1_amp", 0.0);
  const Real piw_mode1_amp = pin->GetOrAddReal("problem/harris_4d", "piw_mode1_amp", 0.0);
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
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");

  const Real x1min = pin->GetReal("parthenon/mesh", "x1min");
  const Real x1max = pin->GetReal("parthenon/mesh", "x1max");
  const Real x3min = pin->GetReal("parthenon/mesh", "x3min");
  const Real x3max = pin->GetReal("parthenon/mesh", "x3max");
  const Real lx = x1max - x1min;
  const Real lz = x3max - x3min;
  const Real gm1 = gamma - 1.0;
  const Real qom_denom = qom_ion - qom_electron;
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

        const Real phase_x = (lx > 0.0) ? (kTwoPi * (x - x1min) / lx) : 0.0;
        const Real phase_z = (lz > 0.0) ? (kTwoPi * (z - x3min) / lz) : 0.0;
        const Real ay_perturb = perturbation_amp * std::cos(phase_x) * std::cos(phase_z);
        em_a(2, k, j, i) = ay_perturb; // mode 0, A_y component

        if (n_modes > 1) {
          const int aw_mode1_idx = 5 + 4;
          const Real phase = std::cos(phase_x) * std::cos(phase_z);
          em_a(aw_mode1_idx, k, j, i) = aw_mode1_amp * phase;
          em_pi(aw_mode1_idx, k, j, i) = piw_mode1_amp * phase;
        }

        // Two-species (ions/electrons) mode state. Initialize only the zero mode.
        const Real jz_sheet = drift_current_scale * (b0 / sheet_half_width) * sech2y;
        const Real rho_ion = (-qom_electron / qom_denom) * rho;
        const Real rho_electron = (qom_ion / qom_denom) * rho;
        const Real momz_ion = jz_sheet / qom_denom;
        const Real momz_electron = -momz_ion;
        const Real energy_density = pressure / gm1;
        const Real energy_ion = (-qom_electron / qom_denom) * energy_density;
        const Real energy_electron = (qom_ion / qom_denom) * energy_density;
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
        plasma(ion_offset + 4, k, j, i) = 0.0;
        plasma(ion_offset + 5, k, j, i) = energy_ion;

        const int ele_offset = 6 * n_modes;
        plasma(ele_offset + 0, k, j, i) = rho_electron;
        plasma(ele_offset + 1, k, j, i) = 0.0;
        plasma(ele_offset + 2, k, j, i) = 0.0;
        plasma(ele_offset + 3, k, j, i) = momz_electron;
        plasma(ele_offset + 4, k, j, i) = 0.0;
        plasma(ele_offset + 5, k, j, i) = energy_electron;
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
