#include "em4d_modes.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "mode_tables.hpp"

namespace Modes4D {
using namespace parthenon::package::prelude;

namespace {

constexpr int kCompA0 = 0;
constexpr int kCompAX = 1;
constexpr int kCompAY = 2;
constexpr int kCompAZ = 3;
constexpr int kCompAW = 4;
constexpr int kNumEMComponents = 5;

constexpr int kSpeciesCount = 2;
constexpr int kVarsPerModePerSpecies = 6;
constexpr int kPlasmaRho = 0;
constexpr int kPlasmaMomX = 1;
constexpr int kPlasmaMomY = 2;
constexpr int kPlasmaMomZ = 3;
constexpr int kPlasmaMomW = 4;

int EMIndex(const int mode, const int component) {
  return kNumEMComponents * mode + component;
}

int PlasmaIndex(const int species, const int mode, const int var, const int n_modes) {
  return (species * kVarsPerModePerSpecies * n_modes) + (mode * kVarsPerModePerSpecies) + var;
}

} // namespace

void RegisterEMVariables(parthenon::StateDescriptor *pkg, const int n_modes) {
  std::vector<std::string> labels(5 * n_modes);
  for (int n = 0; n < n_modes; ++n) {
    labels[5 * n + 0] = "a0_mode_" + std::to_string(n);
    labels[5 * n + 1] = "ax_mode_" + std::to_string(n);
    labels[5 * n + 2] = "ay_mode_" + std::to_string(n);
    labels[5 * n + 3] = "az_mode_" + std::to_string(n);
    labels[5 * n + 4] = "aw_mode_" + std::to_string(n);
  }

  Metadata m({Metadata::Cell, Metadata::Independent, Metadata::FillGhost},
             std::vector<int>({5 * n_modes}), labels);
  pkg->AddField("em4d_a", m);

  std::vector<std::string> pi_labels(5 * n_modes);
  for (int n = 0; n < n_modes; ++n) {
    pi_labels[5 * n + 0] = "pi0_mode_" + std::to_string(n);
    pi_labels[5 * n + 1] = "pix_mode_" + std::to_string(n);
    pi_labels[5 * n + 2] = "piy_mode_" + std::to_string(n);
    pi_labels[5 * n + 3] = "piz_mode_" + std::to_string(n);
    pi_labels[5 * n + 4] = "piw_mode_" + std::to_string(n);
  }

  m = Metadata({Metadata::Cell, Metadata::Independent, Metadata::FillGhost},
               std::vector<int>({5 * n_modes}), pi_labels);
  pkg->AddField("em4d_pi", m);
}

void SourceUnsplit(MeshData<Real> *md, const parthenon::SimTime &, const Real dt) {
  if (dt <= 0.0) {
    return;
  }

  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto modes_pkg = pmb->packages.Get("modes4d");
  if (!modes_pkg->Param<bool>("enabled")) {
    return;
  }

  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real lambda = modes_pkg->Param<double>("lambda");
  const Real c_wave = modes_pkg->Param<double>("em4d/c_wave");
  const Real damping = modes_pkg->Param<double>("em4d/damping");
  const Real mu0 = modes_pkg->Param<double>("em4d/mu0");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");
  const Real c2 = c_wave * c_wave;
  const Real inv_lambda_root2 = std::sqrt(2.0) / lambda;
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const auto &weights = tables.Weights();
  const auto &mass_squared = tables.MassSquared();

  Real diag_jw_ew_step = 0.0;
  Real diag_s_leak_step = 0.0;

  const int num_blocks = md->NumBlocks();
  for (int b = 0; b < num_blocks; ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pblock = bd->GetBlockPointer();

    auto &a_dev = bd->Get("em4d_a").data;
    auto &pi_dev = bd->Get("em4d_pi").data;
    auto &plasma_dev = bd->Get("plasma4d_cons").data;

    auto a_old = a_dev.GetHostMirrorAndCopy();
    auto pi_old = pi_dev.GetHostMirrorAndCopy();
    auto plasma = plasma_dev.GetHostMirrorAndCopy();

    auto a_new = a_old;
    auto pi_new = pi_old;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    const bool has_x = (ib.e > ib.s);
    const bool has_y = (jb.e > jb.s);
    const bool has_z = (kb.e > kb.s);

    auto &coords = pblock->coords;
    auto gradient_x = [&](const auto &field, const int idx, const int k, const int j,
                          const int i) -> Real {
      if (!has_x) return 0.0;
      const Real dx = coords.Dxc<1>(i);
      return (field(idx, k, j, i + 1) - field(idx, k, j, i - 1)) / (2.0 * dx);
    };
    auto gradient_y = [&](const auto &field, const int idx, const int k, const int j,
                          const int i) -> Real {
      if (!has_y) return 0.0;
      const Real dy = coords.Dxc<2>(j);
      return (field(idx, k, j + 1, i) - field(idx, k, j - 1, i)) / (2.0 * dy);
    };
    auto gradient_z = [&](const auto &field, const int idx, const int k, const int j,
                          const int i) -> Real {
      if (!has_z) return 0.0;
      const Real dz = coords.Dxc<3>(k);
      return (field(idx, k + 1, j, i) - field(idx, k - 1, j, i)) / (2.0 * dz);
    };
    auto laplacian = [&](const auto &field, const int idx, const int k, const int j,
                         const int i) -> Real {
      Real lap = 0.0;
      if (has_x) {
        const Real dx = coords.Dxc<1>(i);
        lap += (field(idx, k, j, i + 1) - (2.0 * field(idx, k, j, i)) +
                field(idx, k, j, i - 1)) /
               (dx * dx);
      }
      if (has_y) {
        const Real dy = coords.Dxc<2>(j);
        lap += (field(idx, k, j + 1, i) - (2.0 * field(idx, k, j, i)) +
                field(idx, k, j - 1, i)) /
               (dy * dy);
      }
      if (has_z) {
        const Real dz = coords.Dxc<3>(k);
        lap += (field(idx, k + 1, j, i) - (2.0 * field(idx, k, j, i)) +
                field(idx, k - 1, j, i)) /
               (dz * dz);
      }
      return lap;
    };

    std::vector<Real> j0_modes(n_modes, 0.0);
    std::vector<Real> jx_modes(n_modes, 0.0);
    std::vector<Real> jy_modes(n_modes, 0.0);
    std::vector<Real> jz_modes(n_modes, 0.0);
    std::vector<Real> jw_modes(n_modes, 0.0);
    std::vector<Real> a0_modes(n_modes, 0.0);
    std::vector<Real> ew_modes(n_modes, 0.0);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          for (int n = 0; n < n_modes; ++n) {
            const int ion_base = PlasmaIndex(0, n, 0, n_modes);
            const int ele_base = PlasmaIndex(1, n, 0, n_modes);

            const Real rho_ion = plasma(ion_base + kPlasmaRho, k, j, i);
            const Real rho_electron = plasma(ele_base + kPlasmaRho, k, j, i);
            const Real momx_ion = plasma(ion_base + kPlasmaMomX, k, j, i);
            const Real momx_electron = plasma(ele_base + kPlasmaMomX, k, j, i);
            const Real momy_ion = plasma(ion_base + kPlasmaMomY, k, j, i);
            const Real momy_electron = plasma(ele_base + kPlasmaMomY, k, j, i);
            const Real momz_ion = plasma(ion_base + kPlasmaMomZ, k, j, i);
            const Real momz_electron = plasma(ele_base + kPlasmaMomZ, k, j, i);
            const Real momw_ion = plasma(ion_base + kPlasmaMomW, k, j, i);
            const Real momw_electron = plasma(ele_base + kPlasmaMomW, k, j, i);

            j0_modes[n] = (qom_ion * rho_ion) + (qom_electron * rho_electron);
            jx_modes[n] = (qom_ion * momx_ion) + (qom_electron * momx_electron);
            jy_modes[n] = (qom_ion * momy_ion) + (qom_electron * momy_electron);
            jz_modes[n] = (qom_ion * momz_ion) + (qom_electron * momz_electron);
            jw_modes[n] = (qom_ion * momw_ion) + (qom_electron * momw_electron);
          }

          for (int n = 0; n < n_modes; ++n) {
            const int idx_a0 = EMIndex(n, kCompA0);
            const int idx_ax = EMIndex(n, kCompAX);
            const int idx_ay = EMIndex(n, kCompAY);
            const int idx_az = EMIndex(n, kCompAZ);
            const int idx_aw = EMIndex(n, kCompAW);

            const Real coupling_coeff =
                (n > 0) ? (std::sqrt(2.0 * static_cast<Real>(n)) / lambda) : 0.0;
            Real grad_aw_x = 0.0;
            Real grad_aw_y = 0.0;
            Real grad_aw_z = 0.0;
            if (n > 0) {
              const int idx_aw_prev = EMIndex(n - 1, kCompAW);
              grad_aw_x = gradient_x(a_old, idx_aw_prev, k, j, i);
              grad_aw_y = gradient_y(a_old, idx_aw_prev, k, j, i);
              grad_aw_z = gradient_z(a_old, idx_aw_prev, k, j, i);
            }

            const Real rhs_a0 = (c2 * laplacian(a_old, idx_a0, k, j, i)) -
                                (c2 * mass_squared[n] * a_old(idx_a0, k, j, i)) -
                                (mu0 * j0_modes[n]) - (damping * pi_old(idx_a0, k, j, i));
            const Real rhs_ax = (c2 * laplacian(a_old, idx_ax, k, j, i)) -
                                (c2 * mass_squared[n] * a_old(idx_ax, k, j, i)) +
                                (c2 * coupling_coeff * grad_aw_x) - (mu0 * jx_modes[n]) -
                                (damping * pi_old(idx_ax, k, j, i));
            const Real rhs_ay = (c2 * laplacian(a_old, idx_ay, k, j, i)) -
                                (c2 * mass_squared[n] * a_old(idx_ay, k, j, i)) +
                                (c2 * coupling_coeff * grad_aw_y) - (mu0 * jy_modes[n]) -
                                (damping * pi_old(idx_ay, k, j, i));
            const Real rhs_az = (c2 * laplacian(a_old, idx_az, k, j, i)) -
                                (c2 * mass_squared[n] * a_old(idx_az, k, j, i)) +
                                (c2 * coupling_coeff * grad_aw_z) - (mu0 * jz_modes[n]) -
                                (damping * pi_old(idx_az, k, j, i));
            const Real rhs_aw = (c2 * laplacian(a_old, idx_aw, k, j, i)) -
                                (mu0 * jw_modes[n]) - (damping * pi_old(idx_aw, k, j, i));

            pi_new(idx_a0, k, j, i) = pi_old(idx_a0, k, j, i) + (dt * rhs_a0);
            pi_new(idx_ax, k, j, i) = pi_old(idx_ax, k, j, i) + (dt * rhs_ax);
            pi_new(idx_ay, k, j, i) = pi_old(idx_ay, k, j, i) + (dt * rhs_ay);
            pi_new(idx_az, k, j, i) = pi_old(idx_az, k, j, i) + (dt * rhs_az);
            pi_new(idx_aw, k, j, i) = pi_old(idx_aw, k, j, i) + (dt * rhs_aw);

            a_new(idx_a0, k, j, i) = a_old(idx_a0, k, j, i) + (dt * pi_new(idx_a0, k, j, i));
            a_new(idx_ax, k, j, i) = a_old(idx_ax, k, j, i) + (dt * pi_new(idx_ax, k, j, i));
            a_new(idx_ay, k, j, i) = a_old(idx_ay, k, j, i) + (dt * pi_new(idx_ay, k, j, i));
            a_new(idx_az, k, j, i) = a_old(idx_az, k, j, i) + (dt * pi_new(idx_az, k, j, i));
            a_new(idx_aw, k, j, i) = a_old(idx_aw, k, j, i) + (dt * pi_new(idx_aw, k, j, i));
          }

          for (int n = 0; n < n_modes; ++n) {
            a0_modes[n] = a_new(EMIndex(n, kCompA0), k, j, i);
          }
          for (int n = 0; n < n_modes; ++n) {
            ew_modes[n] = -pi_new(EMIndex(n, kCompAW), k, j, i) -
                          tables.ApplyID3Raising(a0_modes, n);
          }

          Real jw_ew_density = 0.0;
          for (int q = 0; q < tables.NumQuadrature(); ++q) {
            Real jw_at_node = 0.0;
            Real ew_at_node = 0.0;
            for (int n = 0; n < n_modes; ++n) {
              const Real phi_nq = tables.Phi(n, q);
              jw_at_node += jw_modes[n] * phi_nq;
              ew_at_node += ew_modes[n] * phi_nq;
            }
            jw_ew_density += weights[q] * jw_at_node * ew_at_node;
          }

          const Real s_leak_density =
              (n_modes > 1) ? (-(inv_lambda_root2 * jw_modes[1])) : 0.0;
          const Real cell_volume = coords.CellVolume(k, j, i);
          diag_jw_ew_step += dt * cell_volume * jw_ew_density;
          diag_s_leak_step += dt * cell_volume * s_leak_density;
        }
      }
    }

    a_dev.DeepCopy(a_new);
    pi_dev.DeepCopy(pi_new);
  }

  auto *diag_jw_ew = modes_pkg->MutableParam<double>("diag/int_jw_ew");
  auto *diag_s_leak = modes_pkg->MutableParam<double>("diag/int_s_leak");
  *diag_jw_ew += diag_jw_ew_step;
  *diag_s_leak += diag_s_leak_step;
}

} // namespace Modes4D
