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
constexpr Real kTransportSignalRhoFloor = 1.0e-3;
constexpr Real kTransportSignalSpeedCap = 10.0;
constexpr Real kPressureFluxRelativeCap = 50.0;

struct PlasmaState {
  Real rho;
  Real momx;
  Real momy;
  Real momz;
  Real momw;
  Real energy;
  Real pressure;
};

KOKKOS_INLINE_FUNCTION int PlasmaIndex(const int species, const int mode, const int var,
                                       const int n_modes) {
  return (species * kVarsPerModePerSpecies * n_modes) + (mode * kVarsPerModePerSpecies) + var;
}

KOKKOS_INLINE_FUNCTION Real SpeciesPressure(const Real rho, const Real momx, const Real momy,
                                            const Real momz, const Real momw, const Real energy,
                                            const Real rho_floor, const Real gm1,
                                            const Real pressure_floor) {
  const Real rho_safe = fmax(rho, rho_floor);
  const Real kinetic =
      0.5 * ((momx * momx) + (momy * momy) + (momz * momz) + (momw * momw)) / rho_safe;
  const Real internal = fmax(energy - kinetic, 0.0);
  return fmax(pressure_floor, gm1 * internal);
}

KOKKOS_INLINE_FUNCTION Real NormalVelocity(const PlasmaState &state, const int dir,
                                           const Real rho_floor) {
  const Real rho_safe = fmax(state.rho, rho_floor);
  if (dir == X1DIR) {
    return state.momx / rho_safe;
  }
  if (dir == X2DIR) {
    return state.momy / rho_safe;
  }
  return state.momz / rho_safe;
}

KOKKOS_INLINE_FUNCTION Real ClampMagnitude(const Real value, const Real cap) {
  if (value > cap) return cap;
  if (value < -cap) return -cap;
  return value;
}

KOKKOS_INLINE_FUNCTION Real LimitRelativeToAdvective(const Real candidate, const Real advective) {
  constexpr Real kSmall = 1.0e-12;
  const Real cap = kPressureFluxRelativeCap * (fabs(advective) + kSmall);
  return ClampMagnitude(candidate, cap);
}

KOKKOS_INLINE_FUNCTION void ComputeDirectionalFlux(
    const int dir, const PlasmaState &left, const PlasmaState &right, const Real rho_floor,
    const Real gamma, const Real pressure_transport_gain, const Real pressure_rusanov_gain,
    Real &flux_rho, Real &flux_momx, Real &flux_momy, Real &flux_momz, Real &flux_momw,
    Real &flux_energy) {
  const Real vn_l = NormalVelocity(left, dir, rho_floor);
  const Real vn_r = NormalVelocity(right, dir, rho_floor);
  const Real v_face = 0.5 * (vn_l + vn_r);
  const PlasmaState upwind = (v_face >= 0.0) ? left : right;

  // Baseline upwind transport used by existing tuned regressions.
  const Real adv_rho = v_face * upwind.rho;
  const Real adv_momx = v_face * upwind.momx;
  const Real adv_momy = v_face * upwind.momy;
  const Real adv_momz = v_face * upwind.momz;
  const Real adv_momw = v_face * upwind.momw;
  const Real adv_energy = v_face * upwind.energy;

  if (pressure_transport_gain <= 0.0) {
    flux_rho = adv_rho;
    flux_momx = adv_momx;
    flux_momy = adv_momy;
    flux_momz = adv_momz;
    flux_momw = adv_momw;
    flux_energy = adv_energy;
    return;
  }

  const Real blend = fmin(1.0, fmax(0.0, pressure_transport_gain));
  const Real vn_cap = kTransportSignalSpeedCap;
  const Real vn_l_rus = ClampMagnitude(vn_l, vn_cap);
  const Real vn_r_rus = ClampMagnitude(vn_r, vn_cap);
  const Real signal_rho_floor = fmax(rho_floor, kTransportSignalRhoFloor);
  const Real csl =
      sqrt(fmax(0.0, gamma * left.pressure / fmax(left.rho, signal_rho_floor)));
  const Real csr =
      sqrt(fmax(0.0, gamma * right.pressure / fmax(right.rho, signal_rho_floor)));
  const Real alpha_uncapped = fmax(0.0, pressure_rusanov_gain) *
                              fmax(fabs(vn_l_rus) + csl, fabs(vn_r_rus) + csr);
  const Real alpha = fmin(alpha_uncapped, kTransportSignalSpeedCap);

  Real fl_rho = 0.0;
  Real fl_momx = 0.0;
  Real fl_momy = 0.0;
  Real fl_momz = 0.0;
  Real fl_momw = 0.0;
  Real fl_energy = 0.0;
  Real fr_rho = 0.0;
  Real fr_momx = 0.0;
  Real fr_momy = 0.0;
  Real fr_momz = 0.0;
  Real fr_momw = 0.0;
  Real fr_energy = 0.0;
  if (dir == X1DIR) {
    fl_rho = left.momx;
    fl_momx = left.momx * vn_l_rus + left.pressure;
    fl_momy = left.momy * vn_l_rus;
    fl_momz = left.momz * vn_l_rus;
    fl_momw = left.momw * vn_l_rus;
    fl_energy = (left.energy + left.pressure) * vn_l_rus;
    fr_rho = right.momx;
    fr_momx = right.momx * vn_r_rus + right.pressure;
    fr_momy = right.momy * vn_r_rus;
    fr_momz = right.momz * vn_r_rus;
    fr_momw = right.momw * vn_r_rus;
    fr_energy = (right.energy + right.pressure) * vn_r_rus;
  } else if (dir == X2DIR) {
    fl_rho = left.momy;
    fl_momx = left.momx * vn_l_rus;
    fl_momy = left.momy * vn_l_rus + left.pressure;
    fl_momz = left.momz * vn_l_rus;
    fl_momw = left.momw * vn_l_rus;
    fl_energy = (left.energy + left.pressure) * vn_l_rus;
    fr_rho = right.momy;
    fr_momx = right.momx * vn_r_rus;
    fr_momy = right.momy * vn_r_rus + right.pressure;
    fr_momz = right.momz * vn_r_rus;
    fr_momw = right.momw * vn_r_rus;
    fr_energy = (right.energy + right.pressure) * vn_r_rus;
  } else {
    fl_rho = left.momz;
    fl_momx = left.momx * vn_l_rus;
    fl_momy = left.momy * vn_l_rus;
    fl_momz = left.momz * vn_l_rus + left.pressure;
    fl_momw = left.momw * vn_l_rus;
    fl_energy = (left.energy + left.pressure) * vn_l_rus;
    fr_rho = right.momz;
    fr_momx = right.momx * vn_r_rus;
    fr_momy = right.momy * vn_r_rus;
    fr_momz = right.momz * vn_r_rus + right.pressure;
    fr_momw = right.momw * vn_r_rus;
    fr_energy = (right.energy + right.pressure) * vn_r_rus;
  }

  const Real rus_rho = 0.5 * (fl_rho + fr_rho) - 0.5 * alpha * (right.rho - left.rho);
  const Real rus_momx = 0.5 * (fl_momx + fr_momx) - 0.5 * alpha * (right.momx - left.momx);
  const Real rus_momy = 0.5 * (fl_momy + fr_momy) - 0.5 * alpha * (right.momy - left.momy);
  const Real rus_momz = 0.5 * (fl_momz + fr_momz) - 0.5 * alpha * (right.momz - left.momz);
  const Real rus_momw = 0.5 * (fl_momw + fr_momw) - 0.5 * alpha * (right.momw - left.momw);
  const Real rus_energy =
      0.5 * (fl_energy + fr_energy) - 0.5 * alpha * (right.energy - left.energy);

  const Real rus_rho_limited = LimitRelativeToAdvective(rus_rho, adv_rho);
  const Real rus_momx_limited = LimitRelativeToAdvective(rus_momx, adv_momx);
  const Real rus_momy_limited = LimitRelativeToAdvective(rus_momy, adv_momy);
  const Real rus_momz_limited = LimitRelativeToAdvective(rus_momz, adv_momz);
  const Real rus_momw_limited = LimitRelativeToAdvective(rus_momw, adv_momw);
  const Real rus_energy_limited = LimitRelativeToAdvective(rus_energy, adv_energy);

  flux_rho = ((1.0 - blend) * adv_rho) + (blend * rus_rho_limited);
  flux_momx = ((1.0 - blend) * adv_momx) + (blend * rus_momx_limited);
  flux_momy = ((1.0 - blend) * adv_momy) + (blend * rus_momy_limited);
  flux_momz = ((1.0 - blend) * adv_momz) + (blend * rus_momz_limited);
  flux_momw = ((1.0 - blend) * adv_momw) + (blend * rus_momw_limited);
  flux_energy = ((1.0 - blend) * adv_energy) + (blend * rus_energy_limited);
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
  const Real gamma = modes_pkg->Param<double>("plasma4d/gamma");
  const Real gm1 = gamma - 1.0;
  const Real pressure_transport_gain =
      modes_pkg->Param<double>("plasma4d/pressure_transport_gain");
  const Real pressure_rusanov_gain =
      modes_pkg->Param<double>("plasma4d/pressure_rusanov_gain");
  const int pressure_transport_max_mode =
      modes_pkg->Param<int>("plasma4d/pressure_transport_max_mode");
  const Real pressure_floor = modes_pkg->Param<double>("plasma4d/pressure_floor");
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
            const PlasmaState left{
                rho_l, momx_l, plasma(momy_idx, k, j, i - 1), plasma(momz_idx, k, j, i - 1),
                plasma(momw_idx, k, j, i - 1), plasma(energy_idx, k, j, i - 1), 0.0};
            const PlasmaState right{
                rho_r, momx_r, plasma(momy_idx, k, j, i), plasma(momz_idx, k, j, i),
                plasma(momw_idx, k, j, i), plasma(energy_idx, k, j, i), 0.0};
            PlasmaState left_eos = left;
            PlasmaState right_eos = right;
            left_eos.pressure =
                SpeciesPressure(left.rho, left.momx, left.momy, left.momz, left.momw,
                                left.energy, rho_floor, gm1, pressure_floor);
            right_eos.pressure =
                SpeciesPressure(right.rho, right.momx, right.momy, right.momz, right.momw,
                                right.energy, rho_floor, gm1, pressure_floor);
            Real flux_rho = 0.0;
            Real flux_momx = 0.0;
            Real flux_momy = 0.0;
            Real flux_momz = 0.0;
            Real flux_momw = 0.0;
            Real flux_energy = 0.0;
            const Real mode_pressure_gain =
                (n <= pressure_transport_max_mode) ? pressure_transport_gain : 0.0;
            ComputeDirectionalFlux(X1DIR, left_eos, right_eos, rho_floor, gamma,
                                   mode_pressure_gain, pressure_rusanov_gain, flux_rho,
                                   flux_momx, flux_momy, flux_momz, flux_momw, flux_energy);
            plasma.flux(X1DIR, rho_idx, k, j, i) = transport_gain * flux_rho;
            plasma.flux(X1DIR, momx_idx, k, j, i) = transport_gain * flux_momx;
            plasma.flux(X1DIR, momy_idx, k, j, i) = transport_gain * flux_momy;
            plasma.flux(X1DIR, momz_idx, k, j, i) = transport_gain * flux_momz;
            plasma.flux(X1DIR, momw_idx, k, j, i) = transport_gain * flux_momw;
            plasma.flux(X1DIR, energy_idx, k, j, i) = transport_gain * flux_energy;
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
              const PlasmaState left{
                  rho_l, plasma(momx_idx, k, j - 1, i), momy_l, plasma(momz_idx, k, j - 1, i),
                  plasma(momw_idx, k, j - 1, i), plasma(energy_idx, k, j - 1, i), 0.0};
              const PlasmaState right{
                  rho_r, plasma(momx_idx, k, j, i), momy_r, plasma(momz_idx, k, j, i),
                  plasma(momw_idx, k, j, i), plasma(energy_idx, k, j, i), 0.0};
              PlasmaState left_eos = left;
              PlasmaState right_eos = right;
              left_eos.pressure =
                  SpeciesPressure(left.rho, left.momx, left.momy, left.momz, left.momw,
                                  left.energy, rho_floor, gm1, pressure_floor);
              right_eos.pressure =
                  SpeciesPressure(right.rho, right.momx, right.momy, right.momz, right.momw,
                                  right.energy, rho_floor, gm1, pressure_floor);
              Real flux_rho = 0.0;
              Real flux_momx = 0.0;
              Real flux_momy = 0.0;
              Real flux_momz = 0.0;
              Real flux_momw = 0.0;
              Real flux_energy = 0.0;
              const Real mode_pressure_gain =
                  (n <= pressure_transport_max_mode) ? pressure_transport_gain : 0.0;
              ComputeDirectionalFlux(X2DIR, left_eos, right_eos, rho_floor, gamma,
                                     mode_pressure_gain, pressure_rusanov_gain, flux_rho,
                                     flux_momx, flux_momy, flux_momz, flux_momw, flux_energy);
              plasma.flux(X2DIR, rho_idx, k, j, i) = transport_gain * flux_rho;
              plasma.flux(X2DIR, momx_idx, k, j, i) = transport_gain * flux_momx;
              plasma.flux(X2DIR, momy_idx, k, j, i) = transport_gain * flux_momy;
              plasma.flux(X2DIR, momz_idx, k, j, i) = transport_gain * flux_momz;
              plasma.flux(X2DIR, momw_idx, k, j, i) = transport_gain * flux_momw;
              plasma.flux(X2DIR, energy_idx, k, j, i) = transport_gain * flux_energy;
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
              const PlasmaState left{
                  rho_l, plasma(momx_idx, k - 1, j, i), plasma(momy_idx, k - 1, j, i), momz_l,
                  plasma(momw_idx, k - 1, j, i), plasma(energy_idx, k - 1, j, i), 0.0};
              const PlasmaState right{
                  rho_r, plasma(momx_idx, k, j, i), plasma(momy_idx, k, j, i), momz_r,
                  plasma(momw_idx, k, j, i), plasma(energy_idx, k, j, i), 0.0};
              PlasmaState left_eos = left;
              PlasmaState right_eos = right;
              left_eos.pressure =
                  SpeciesPressure(left.rho, left.momx, left.momy, left.momz, left.momw,
                                  left.energy, rho_floor, gm1, pressure_floor);
              right_eos.pressure =
                  SpeciesPressure(right.rho, right.momx, right.momy, right.momz, right.momw,
                                  right.energy, rho_floor, gm1, pressure_floor);
              Real flux_rho = 0.0;
              Real flux_momx = 0.0;
              Real flux_momy = 0.0;
              Real flux_momz = 0.0;
              Real flux_momw = 0.0;
              Real flux_energy = 0.0;
              const Real mode_pressure_gain =
                  (n <= pressure_transport_max_mode) ? pressure_transport_gain : 0.0;
              ComputeDirectionalFlux(X3DIR, left_eos, right_eos, rho_floor, gamma,
                                     mode_pressure_gain, pressure_rusanov_gain, flux_rho,
                                     flux_momx, flux_momy, flux_momz, flux_momw, flux_energy);
              plasma.flux(X3DIR, rho_idx, k, j, i) = transport_gain * flux_rho;
              plasma.flux(X3DIR, momx_idx, k, j, i) = transport_gain * flux_momx;
              plasma.flux(X3DIR, momy_idx, k, j, i) = transport_gain * flux_momy;
              plasma.flux(X3DIR, momz_idx, k, j, i) = transport_gain * flux_momz;
              plasma.flux(X3DIR, momw_idx, k, j, i) = transport_gain * flux_momw;
              plasma.flux(X3DIR, energy_idx, k, j, i) = transport_gain * flux_energy;
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
