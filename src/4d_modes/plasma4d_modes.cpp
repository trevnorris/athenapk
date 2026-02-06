#include "plasma4d_modes.hpp"

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

TaskStatus AddRhoTransportFluxes(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto modes_pkg = pmb->packages.Get("modes4d");
  if (!modes_pkg->Param<bool>("enabled")) {
    return TaskStatus::complete;
  }

  const int n_modes = modes_pkg->Param<int>("n_modes");
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
            plasma.flux(X1DIR, rho_idx, k, j, i) =
                0.5 * rho_transport_gain *
                (plasma(momx_idx, k, j, i - 1) + plasma(momx_idx, k, j, i));
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
              const int momy_idx = PlasmaIndex(s, n, kPlasmaMomY, n_modes);
              plasma.flux(X2DIR, rho_idx, k, j, i) =
                  0.5 * rho_transport_gain *
                  (plasma(momy_idx, k, j - 1, i) + plasma(momy_idx, k, j, i));
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
              const int momz_idx = PlasmaIndex(s, n, kPlasmaMomZ, n_modes);
              plasma.flux(X3DIR, rho_idx, k, j, i) =
                  0.5 * rho_transport_gain *
                  (plasma(momz_idx, k - 1, j, i) + plasma(momz_idx, k, j, i));
            }
          }
        });
  }

  return TaskStatus::complete;
}

TaskStatus AccumulateRhoTransportDivJ(MeshData<Real> *md, const Real dt) {
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

  Real divj_step = 0.0;
  parthenon::par_reduce(
      DEFAULT_LOOP_PATTERN, "Modes4DAccumulateRhoTransportDivJ", parthenon::DevExecSpace(),
      0, plasma_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &local_sum) {
        const auto &plasma = plasma_pack(b);
        const auto &coords = plasma_pack.GetCoords(b);
        const Real ion_flux_div =
            parthenon::Update::FluxDivHelper(ion_rho_mode0, k, j, i, ndim, coords, plasma);
        const Real ele_flux_div =
            parthenon::Update::FluxDivHelper(ele_rho_mode0, k, j, i, ndim, coords, plasma);
        const Real div_charge = -((qom_ion * ion_flux_div) + (qom_electron * ele_flux_div));
        local_sum += dt * coords.CellVolume(k, j, i) * div_charge;
      },
      divj_step);

  auto *diag_divj_mode0 = modes_pkg->MutableParam<double>("diag/int_divj_mode0");
  *diag_divj_mode0 += divj_step;
  return TaskStatus::complete;
}

} // namespace Modes4D
