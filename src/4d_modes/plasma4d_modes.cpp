#include "plasma4d_modes.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "interface/update.hpp"

namespace Modes4D {
using namespace parthenon::package::prelude;

namespace {

constexpr int kSpeciesCount = 2;
constexpr int kVarsPerModePerSpecies = 6;
constexpr int kPlasmaRho = 0;
constexpr int kPlasmaMomX = 1;
constexpr int kPlasmaMomY = 2;
constexpr int kPlasmaMomZ = 3;
constexpr int kPlasmaMomW = 4;
constexpr int kPlasmaEnergy = 5;

KOKKOS_INLINE_FUNCTION int PlasmaIndex(const int species, const int mode, const int var,
                                       const int n_modes) {
  return (species * kVarsPerModePerSpecies * n_modes) + (mode * kVarsPerModePerSpecies) + var;
}

} // namespace

void RegisterPlasmaVariables(parthenon::StateDescriptor *pkg, const int n_modes) {
  constexpr int kSpecies = 2; // ions + electrons in the initial design.

  std::vector<std::string> labels(kSpecies * kVarsPerModePerSpecies * n_modes);
  int idx = 0;
  for (int s = 0; s < kSpecies; ++s) {
    const std::string prefix = (s == 0) ? "ion" : "electron";
    for (int n = 0; n < n_modes; ++n) {
      labels[idx++] = prefix + "_rho_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momx_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momy_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momz_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momw_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_energy_mode_" + std::to_string(n);
    }
  }

  Metadata m({Metadata::Cell, Metadata::Independent, Metadata::FillGhost,
              Metadata::WithFluxes},
             std::vector<int>({static_cast<int>(labels.size())}), labels);
  pkg->AddField("plasma4d_cons", m);
}

TaskStatus AddPlasmaTransportFluxes(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto modes_pkg = pmb->packages.Get("modes4d");
  if (!modes_pkg->Param<bool>("enabled")) {
    return TaskStatus::complete;
  }

  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real transport_gain = modes_pkg->Param<double>("plasma4d/rho_divj_gain");
  const Real rho_floor = modes_pkg->Param<double>("plasma4d/rho_floor");
  if (n_modes < 1 || transport_gain == 0.0) {
    return TaskStatus::complete;
  }

  const auto &plasma_pack =
      md->PackVariablesAndFluxes(std::vector<std::string>{"plasma4d_cons"},
                                 std::vector<std::string>{"plasma4d_cons"});
  if (plasma_pack.GetDim(4) == 0) {
    return TaskStatus::complete;
  }

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);
  const int ndim = pmb->pmy_mesh->ndim;

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Modes4DPlasmaRhoFluxX1", parthenon::DevExecSpace(), 0,
      plasma_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        auto &plasma = plasma_pack(b);
        for (int s = 0; s < kSpeciesCount; ++s) {
          for (int n = 0; n < n_modes; ++n) {
            const int rho_idx = PlasmaIndex(s, n, kPlasmaRho, n_modes);
            const int momx_idx = PlasmaIndex(s, n, kPlasmaMomX, n_modes);
            const int momy_idx = PlasmaIndex(s, n, kPlasmaMomY, n_modes);
            const int momz_idx = PlasmaIndex(s, n, kPlasmaMomZ, n_modes);
            const int momw_idx = PlasmaIndex(s, n, kPlasmaMomW, n_modes);
            const int energy_idx = PlasmaIndex(s, n, kPlasmaEnergy, n_modes);

            const Real rho_l = plasma(rho_idx, k, j, i - 1);
            const Real rho_r = plasma(rho_idx, k, j, i);
            const Real momx_l = plasma(momx_idx, k, j, i - 1);
            const Real momx_r = plasma(momx_idx, k, j, i);
            const Real vx_l = momx_l / fmax(rho_l, rho_floor);
            const Real vx_r = momx_r / fmax(rho_r, rho_floor);
            const Real v_face = 0.5 * (vx_l + vx_r);
            const bool take_left = (v_face >= 0.0);
            const int il = take_left ? (i - 1) : i;

            plasma.flux(X1DIR, rho_idx, k, j, i) = transport_gain * v_face * plasma(rho_idx, k, j, il);
            plasma.flux(X1DIR, momx_idx, k, j, i) =
                transport_gain * v_face * plasma(momx_idx, k, j, il);
            plasma.flux(X1DIR, momy_idx, k, j, i) =
                transport_gain * v_face * plasma(momy_idx, k, j, il);
            plasma.flux(X1DIR, momz_idx, k, j, i) =
                transport_gain * v_face * plasma(momz_idx, k, j, il);
            plasma.flux(X1DIR, momw_idx, k, j, i) =
                transport_gain * v_face * plasma(momw_idx, k, j, il);
            plasma.flux(X1DIR, energy_idx, k, j, i) =
                transport_gain * v_face * plasma(energy_idx, k, j, il);
          }
        }
      });

  if (ndim >= 2) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "Modes4DPlasmaRhoFluxX2", parthenon::DevExecSpace(), 0,
        plasma_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e + 1, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          auto &plasma = plasma_pack(b);
          for (int s = 0; s < kSpeciesCount; ++s) {
            for (int n = 0; n < n_modes; ++n) {
              const int rho_idx = PlasmaIndex(s, n, kPlasmaRho, n_modes);
              const int momx_idx = PlasmaIndex(s, n, kPlasmaMomX, n_modes);
              const int momy_idx = PlasmaIndex(s, n, kPlasmaMomY, n_modes);
              const int momz_idx = PlasmaIndex(s, n, kPlasmaMomZ, n_modes);
              const int momw_idx = PlasmaIndex(s, n, kPlasmaMomW, n_modes);
              const int energy_idx = PlasmaIndex(s, n, kPlasmaEnergy, n_modes);

              const Real rho_l = plasma(rho_idx, k, j - 1, i);
              const Real rho_r = plasma(rho_idx, k, j, i);
              const Real momy_l = plasma(momy_idx, k, j - 1, i);
              const Real momy_r = plasma(momy_idx, k, j, i);
              const Real vy_l = momy_l / fmax(rho_l, rho_floor);
              const Real vy_r = momy_r / fmax(rho_r, rho_floor);
              const Real v_face = 0.5 * (vy_l + vy_r);
              const bool take_left = (v_face >= 0.0);
              const int jl = take_left ? (j - 1) : j;

              plasma.flux(X2DIR, rho_idx, k, j, i) = transport_gain * v_face * plasma(rho_idx, k, jl, i);
              plasma.flux(X2DIR, momx_idx, k, j, i) =
                  transport_gain * v_face * plasma(momx_idx, k, jl, i);
              plasma.flux(X2DIR, momy_idx, k, j, i) =
                  transport_gain * v_face * plasma(momy_idx, k, jl, i);
              plasma.flux(X2DIR, momz_idx, k, j, i) =
                  transport_gain * v_face * plasma(momz_idx, k, jl, i);
              plasma.flux(X2DIR, momw_idx, k, j, i) =
                  transport_gain * v_face * plasma(momw_idx, k, jl, i);
              plasma.flux(X2DIR, energy_idx, k, j, i) =
                  transport_gain * v_face * plasma(energy_idx, k, jl, i);
            }
          }
        });
  }

  if (ndim >= 3) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "Modes4DPlasmaRhoFluxX3", parthenon::DevExecSpace(), 0,
        plasma_pack.GetDim(5) - 1, kb.s, kb.e + 1, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          auto &plasma = plasma_pack(b);
          for (int s = 0; s < kSpeciesCount; ++s) {
            for (int n = 0; n < n_modes; ++n) {
              const int rho_idx = PlasmaIndex(s, n, kPlasmaRho, n_modes);
              const int momx_idx = PlasmaIndex(s, n, kPlasmaMomX, n_modes);
              const int momy_idx = PlasmaIndex(s, n, kPlasmaMomY, n_modes);
              const int momz_idx = PlasmaIndex(s, n, kPlasmaMomZ, n_modes);
              const int momw_idx = PlasmaIndex(s, n, kPlasmaMomW, n_modes);
              const int energy_idx = PlasmaIndex(s, n, kPlasmaEnergy, n_modes);

              const Real rho_l = plasma(rho_idx, k - 1, j, i);
              const Real rho_r = plasma(rho_idx, k, j, i);
              const Real momz_l = plasma(momz_idx, k - 1, j, i);
              const Real momz_r = plasma(momz_idx, k, j, i);
              const Real vz_l = momz_l / fmax(rho_l, rho_floor);
              const Real vz_r = momz_r / fmax(rho_r, rho_floor);
              const Real v_face = 0.5 * (vz_l + vz_r);
              const bool take_left = (v_face >= 0.0);
              const int kl = take_left ? (k - 1) : k;

              plasma.flux(X3DIR, rho_idx, k, j, i) = transport_gain * v_face * plasma(rho_idx, kl, j, i);
              plasma.flux(X3DIR, momx_idx, k, j, i) =
                  transport_gain * v_face * plasma(momx_idx, kl, j, i);
              plasma.flux(X3DIR, momy_idx, k, j, i) =
                  transport_gain * v_face * plasma(momy_idx, kl, j, i);
              plasma.flux(X3DIR, momz_idx, k, j, i) =
                  transport_gain * v_face * plasma(momz_idx, kl, j, i);
              plasma.flux(X3DIR, momw_idx, k, j, i) =
                  transport_gain * v_face * plasma(momw_idx, kl, j, i);
              plasma.flux(X3DIR, energy_idx, k, j, i) =
                  transport_gain * v_face * plasma(energy_idx, kl, j, i);
            }
          }
        });
  }

  return TaskStatus::complete;
}

TaskStatus AccumulateTransportMode0Diagnostics(MeshData<Real> *md, const Real dt) {
  if (dt <= 0.0) {
    return TaskStatus::complete;
  }

  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto modes_pkg = pmb->packages.Get("modes4d");
  if (!modes_pkg->Param<bool>("enabled")) {
    return TaskStatus::complete;
  }

  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");
  const Real rho_transport_gain = modes_pkg->Param<double>("plasma4d/rho_divj_gain");
  if (n_modes < 1 || rho_transport_gain == 0.0) {
    return TaskStatus::complete;
  }

  const auto &plasma_pack =
      md->PackVariablesAndFluxes(std::vector<std::string>{"plasma4d_cons"},
                                 std::vector<std::string>{"plasma4d_cons"});
  if (plasma_pack.GetDim(4) == 0) {
    return TaskStatus::complete;
  }

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);
  const int ndim = plasma_pack.GetNdim();
  const int ion_rho_mode0 = PlasmaIndex(0, 0, kPlasmaRho, n_modes);
  const int ele_rho_mode0 = PlasmaIndex(1, 0, kPlasmaRho, n_modes);
  const int ion_momx_mode0 = PlasmaIndex(0, 0, kPlasmaMomX, n_modes);
  const int ele_momx_mode0 = PlasmaIndex(1, 0, kPlasmaMomX, n_modes);
  const int ion_momy_mode0 = PlasmaIndex(0, 0, kPlasmaMomY, n_modes);
  const int ele_momy_mode0 = PlasmaIndex(1, 0, kPlasmaMomY, n_modes);
  const int ion_momz_mode0 = PlasmaIndex(0, 0, kPlasmaMomZ, n_modes);
  const int ele_momz_mode0 = PlasmaIndex(1, 0, kPlasmaMomZ, n_modes);
  const int ion_momw_mode0 = PlasmaIndex(0, 0, kPlasmaMomW, n_modes);
  const int ele_momw_mode0 = PlasmaIndex(1, 0, kPlasmaMomW, n_modes);
  const int ion_energy_mode0 = PlasmaIndex(0, 0, kPlasmaEnergy, n_modes);
  const int ele_energy_mode0 = PlasmaIndex(1, 0, kPlasmaEnergy, n_modes);

  auto accumulate_transport = [&](const int ion_idx, const int ele_idx, const Real ion_weight,
                                  const Real ele_weight, const char *kernel_name) {
    Real step = 0.0;
    parthenon::par_reduce(
        DEFAULT_LOOP_PATTERN, kernel_name, parthenon::DevExecSpace(), 0,
        plasma_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &local_sum) {
          const auto &plasma = plasma_pack(b);
          const auto &coords = plasma_pack.GetCoords(b);
          const Real ion_flux_div =
              parthenon::Update::FluxDivHelper(ion_idx, k, j, i, ndim, coords, plasma);
          const Real ele_flux_div =
              parthenon::Update::FluxDivHelper(ele_idx, k, j, i, ndim, coords, plasma);
          const Real div_term = -((ion_weight * ion_flux_div) + (ele_weight * ele_flux_div));
          local_sum += dt * coords.CellVolume(k, j, i) * div_term;
        },
        step);
    return step;
  };

  const Real divj_step = accumulate_transport(ion_rho_mode0, ele_rho_mode0, qom_ion,
                                              qom_electron, "Modes4DAccumulateDivJMode0");
  const Real divmomx_step =
      accumulate_transport(ion_momx_mode0, ele_momx_mode0, 1.0, 1.0,
                           "Modes4DAccumulateDivMomXMode0");
  const Real divmomy_step =
      accumulate_transport(ion_momy_mode0, ele_momy_mode0, 1.0, 1.0,
                           "Modes4DAccumulateDivMomYMode0");
  const Real divmomz_step =
      accumulate_transport(ion_momz_mode0, ele_momz_mode0, 1.0, 1.0,
                           "Modes4DAccumulateDivMomZMode0");
  const Real divmomw_step =
      accumulate_transport(ion_momw_mode0, ele_momw_mode0, 1.0, 1.0,
                           "Modes4DAccumulateDivMomWMode0");
  const Real divenergy_step =
      accumulate_transport(ion_energy_mode0, ele_energy_mode0, 1.0, 1.0,
                           "Modes4DAccumulateDivEnergyMode0");

  auto *diag_divj_mode0 = modes_pkg->MutableParam<double>("diag/int_divj_mode0");
  auto *diag_divmomx_mode0 = modes_pkg->MutableParam<double>("diag/int_divmomx_mode0");
  auto *diag_divmomy_mode0 = modes_pkg->MutableParam<double>("diag/int_divmomy_mode0");
  auto *diag_divmomz_mode0 = modes_pkg->MutableParam<double>("diag/int_divmomz_mode0");
  auto *diag_divmomw_mode0 = modes_pkg->MutableParam<double>("diag/int_divmomw_mode0");
  auto *diag_divenergy_mode0 = modes_pkg->MutableParam<double>("diag/int_divenergy_mode0");
  *diag_divj_mode0 += divj_step;
  *diag_divmomx_mode0 += divmomx_step;
  *diag_divmomy_mode0 += divmomy_step;
  *diag_divmomz_mode0 += divmomz_step;
  *diag_divmomw_mode0 += divmomw_step;
  *diag_divenergy_mode0 += divenergy_step;
  return TaskStatus::complete;
}

} // namespace Modes4D
