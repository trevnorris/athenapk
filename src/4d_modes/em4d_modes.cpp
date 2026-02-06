#include "em4d_modes.hpp"

#include <algorithm>
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
constexpr int kPlasmaEnergy = 5;

Real SpeciesQOM(const int species, const Real qom_ion, const Real qom_electron) {
  return (species == 0) ? qom_ion : qom_electron;
}

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
  const Real force_source_gain = modes_pkg->Param<double>("plasma4d/force_source_gain");
  const Real momw_source_gain = modes_pkg->Param<double>("plasma4d/momw_source_gain");
  const Real momw_damping = modes_pkg->Param<double>("plasma4d/momw_damping");
  const Real rho_floor = modes_pkg->Param<double>("plasma4d/rho_floor");
  const Real energy_source_gain = modes_pkg->Param<double>("plasma4d/energy_source_gain");
  const Real energy_floor = modes_pkg->Param<double>("plasma4d/energy_floor");
  const Real c2 = c_wave * c_wave;
  const Real inv_lambda_root2 = std::sqrt(2.0) / lambda;
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const auto &weights = tables.Weights();
  const auto &mass_squared = tables.MassSquared();

  Real diag_jw_ew_step = 0.0;
  Real diag_s_leak_step = 0.0;
  Real diag_s_leak_abs_step = 0.0;
  Real diag_cont_local_l1_step = 0.0;
  Real diag_cont_local_l2_step = 0.0;
  Real diag_cont_local_max_abs_step = 0.0;
  Real diag_cont_mode0_l1_step = 0.0;
  Real diag_cont_mode0_l2_step = 0.0;
  Real diag_cont_mode0_max_abs_step = 0.0;

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
    auto plasma_new = plasma;

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
    std::vector<Real> ax_modes(n_modes, 0.0);
    std::vector<Real> ay_modes(n_modes, 0.0);
    std::vector<Real> az_modes(n_modes, 0.0);
    std::vector<Real> piw_modes(n_modes, 0.0);
    std::vector<Real> pix_modes(n_modes, 0.0);
    std::vector<Real> piy_modes(n_modes, 0.0);
    std::vector<Real> piz_modes(n_modes, 0.0);
    std::vector<Real> da0_dx_modes(n_modes, 0.0);
    std::vector<Real> da0_dy_modes(n_modes, 0.0);
    std::vector<Real> da0_dz_modes(n_modes, 0.0);
    std::vector<Real> dax_dy_modes(n_modes, 0.0);
    std::vector<Real> dax_dz_modes(n_modes, 0.0);
    std::vector<Real> day_dx_modes(n_modes, 0.0);
    std::vector<Real> day_dz_modes(n_modes, 0.0);
    std::vector<Real> daz_dx_modes(n_modes, 0.0);
    std::vector<Real> daz_dy_modes(n_modes, 0.0);
    std::vector<Real> daw_dx_modes(n_modes, 0.0);
    std::vector<Real> daw_dy_modes(n_modes, 0.0);
    std::vector<Real> daw_dz_modes(n_modes, 0.0);
    std::vector<Real> ex_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> ey_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> ez_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> bx_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> by_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> bz_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> ew_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> cx_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> cy_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> cz_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> rho_modes(n_modes, 0.0);
    std::vector<Real> momx_modes(n_modes, 0.0);
    std::vector<Real> momy_modes(n_modes, 0.0);
    std::vector<Real> momz_modes(n_modes, 0.0);
    std::vector<Real> momw_modes(n_modes, 0.0);
    std::vector<Real> energy_modes(n_modes, 0.0);
    std::vector<Real> momw_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> energy_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> momx_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> momy_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> momz_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> momw_modes_new(n_modes, 0.0);
    std::vector<Real> momx_modes_new(n_modes, 0.0);
    std::vector<Real> momy_modes_new(n_modes, 0.0);
    std::vector<Real> momz_modes_new(n_modes, 0.0);
    std::vector<Real> energy_modes_new(n_modes, 0.0);
    std::vector<Real> energy_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> charge_modes_old(n_modes, 0.0);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          for (int n = 0; n < n_modes; ++n) {
            const int ion_base = PlasmaIndex(0, n, 0, n_modes);
            const int ele_base = PlasmaIndex(1, n, 0, n_modes);
            const Real ion_rho = plasma(ion_base + kPlasmaRho, k, j, i);
            const Real ele_rho = plasma(ele_base + kPlasmaRho, k, j, i);
            charge_modes_old[n] = (qom_ion * ion_rho) + (qom_electron * ele_rho);
          }

          // Update transverse species momentum from mixed-sector 4D forcing:
          // m_s dv_w/dt = q_s (E_w - v^a C_a), represented as mode coefficients.
          for (int n = 0; n < n_modes; ++n) {
            const int off = EMIndex(n, 0);
            a0_modes[n] = a_old(off + kCompA0, k, j, i);
            ax_modes[n] = a_old(off + kCompAX, k, j, i);
            ay_modes[n] = a_old(off + kCompAY, k, j, i);
            az_modes[n] = a_old(off + kCompAZ, k, j, i);
            piw_modes[n] = pi_old(off + kCompAW, k, j, i);
            pix_modes[n] = pi_old(off + kCompAX, k, j, i);
            piy_modes[n] = pi_old(off + kCompAY, k, j, i);
            piz_modes[n] = pi_old(off + kCompAZ, k, j, i);
            da0_dx_modes[n] = gradient_x(a_old, off + kCompA0, k, j, i);
            da0_dy_modes[n] = gradient_y(a_old, off + kCompA0, k, j, i);
            da0_dz_modes[n] = gradient_z(a_old, off + kCompA0, k, j, i);
            dax_dy_modes[n] = gradient_y(a_old, off + kCompAX, k, j, i);
            dax_dz_modes[n] = gradient_z(a_old, off + kCompAX, k, j, i);
            day_dx_modes[n] = gradient_x(a_old, off + kCompAY, k, j, i);
            day_dz_modes[n] = gradient_z(a_old, off + kCompAY, k, j, i);
            daz_dx_modes[n] = gradient_x(a_old, off + kCompAZ, k, j, i);
            daz_dy_modes[n] = gradient_y(a_old, off + kCompAZ, k, j, i);
            daw_dx_modes[n] = gradient_x(a_old, off + kCompAW, k, j, i);
            daw_dy_modes[n] = gradient_y(a_old, off + kCompAW, k, j, i);
            daw_dz_modes[n] = gradient_z(a_old, off + kCompAW, k, j, i);
          }

          for (int q = 0; q < tables.NumQuadrature(); ++q) {
            Real d_w_a0 = 0.0;
            Real d_w_ax = 0.0;
            Real d_w_ay = 0.0;
            Real d_w_az = 0.0;
            Real piw_node = 0.0;
            Real pix_node = 0.0;
            Real piy_node = 0.0;
            Real piz_node = 0.0;
            Real da0_dx_node = 0.0;
            Real da0_dy_node = 0.0;
            Real da0_dz_node = 0.0;
            Real dax_dy_node = 0.0;
            Real dax_dz_node = 0.0;
            Real day_dx_node = 0.0;
            Real day_dz_node = 0.0;
            Real daz_dx_node = 0.0;
            Real daz_dy_node = 0.0;
            Real daw_dx_node = 0.0;
            Real daw_dy_node = 0.0;
            Real daw_dz_node = 0.0;

            for (int n = 0; n < n_modes; ++n) {
              const Real phi_nq = tables.Phi(n, q);
              const Real dphi_nq = tables.DPhi(n, q);
              d_w_a0 += a0_modes[n] * dphi_nq;
              d_w_ax += ax_modes[n] * dphi_nq;
              d_w_ay += ay_modes[n] * dphi_nq;
              d_w_az += az_modes[n] * dphi_nq;
              piw_node += piw_modes[n] * phi_nq;
              pix_node += pix_modes[n] * phi_nq;
              piy_node += piy_modes[n] * phi_nq;
              piz_node += piz_modes[n] * phi_nq;
              da0_dx_node += da0_dx_modes[n] * phi_nq;
              da0_dy_node += da0_dy_modes[n] * phi_nq;
              da0_dz_node += da0_dz_modes[n] * phi_nq;
              dax_dy_node += dax_dy_modes[n] * phi_nq;
              dax_dz_node += dax_dz_modes[n] * phi_nq;
              day_dx_node += day_dx_modes[n] * phi_nq;
              day_dz_node += day_dz_modes[n] * phi_nq;
              daz_dx_node += daz_dx_modes[n] * phi_nq;
              daz_dy_node += daz_dy_modes[n] * phi_nq;
              daw_dx_node += daw_dx_modes[n] * phi_nq;
              daw_dy_node += daw_dy_modes[n] * phi_nq;
              daw_dz_node += daw_dz_modes[n] * phi_nq;
            }

            ex_nodes[q] = -pix_node - da0_dx_node;
            ey_nodes[q] = -piy_node - da0_dy_node;
            ez_nodes[q] = -piz_node - da0_dz_node;
            bx_nodes[q] = daz_dy_node - day_dz_node;
            by_nodes[q] = dax_dz_node - daz_dx_node;
            bz_nodes[q] = day_dx_node - dax_dy_node;
            ew_nodes[q] = -piw_node - d_w_a0;
            cx_nodes[q] = daw_dx_node - d_w_ax;
            cy_nodes[q] = daw_dy_node - d_w_ay;
            cz_nodes[q] = daw_dz_node - d_w_az;
          }

          for (int s = 0; s < kSpeciesCount; ++s) {
            const Real qom_s = SpeciesQOM(s, qom_ion, qom_electron);
            for (int n = 0; n < n_modes; ++n) {
              const int base = PlasmaIndex(s, n, 0, n_modes);
              rho_modes[n] = plasma(base + kPlasmaRho, k, j, i);
              momx_modes[n] = plasma(base + kPlasmaMomX, k, j, i);
              momy_modes[n] = plasma(base + kPlasmaMomY, k, j, i);
              momz_modes[n] = plasma(base + kPlasmaMomZ, k, j, i);
              momw_modes[n] = plasma(base + kPlasmaMomW, k, j, i);
              energy_modes[n] = plasma(base + kPlasmaEnergy, k, j, i);
            }

            for (int q = 0; q < tables.NumQuadrature(); ++q) {
              Real rho_node = 0.0;
              Real momx_node = 0.0;
              Real momy_node = 0.0;
              Real momz_node = 0.0;
              Real momw_node = 0.0;
              Real energy_node = 0.0;
              for (int n = 0; n < n_modes; ++n) {
                const Real phi_nq = tables.Phi(n, q);
                rho_node += rho_modes[n] * phi_nq;
                momx_node += momx_modes[n] * phi_nq;
                momy_node += momy_modes[n] * phi_nq;
                momz_node += momz_modes[n] * phi_nq;
                momw_node += momw_modes[n] * phi_nq;
                energy_node += energy_modes[n] * phi_nq;
              }

              const Real rho_safe = std::max(rho_node, rho_floor);
              const Real vx = momx_node / rho_safe;
              const Real vy = momy_node / rho_safe;
              const Real vz = momz_node / rho_safe;
              const Real vw = momw_node / rho_safe;
              const Real ex = ex_nodes[q];
              const Real ey = ey_nodes[q];
              const Real ez = ez_nodes[q];
              const Real bx = bx_nodes[q];
              const Real by = by_nodes[q];
              const Real bz = bz_nodes[q];
              const Real cx = cx_nodes[q];
              const Real cy = cy_nodes[q];
              const Real cz = cz_nodes[q];
              const Real ew = ew_nodes[q];
              const Real fx = ex + ((vy * bz) - (vz * by)) + (vw * cx);
              const Real fy = ey + ((vz * bx) - (vx * bz)) + (vw * cy);
              const Real fz = ez + ((vx * by) - (vy * bx)) + (vw * cz);
              const Real v_dot_c =
                  (vx * cx) + (vy * cy) + (vz * cz);
              const Real force_w = qom_s * rho_safe * (ew - v_dot_c);
              const Real rhs_momx = force_source_gain * qom_s * rho_safe * fx;
              const Real rhs_momy = force_source_gain * qom_s * rho_safe * fy;
              const Real rhs_momz = force_source_gain * qom_s * rho_safe * fz;
              const Real rhs_momw =
                  (momw_source_gain * force_w) - (momw_damping * momw_node);
              const Real rhs_energy = energy_source_gain * qom_s * rho_safe *
                                      ((vx * ex) + (vy * ey) + (vz * ez) + (vw * ew));

              momw_nodes[q] = momw_node + (dt * rhs_momw);
              energy_nodes[q] = energy_node;
              momx_nodes_new[q] = momx_node + (dt * rhs_momx);
              momy_nodes_new[q] = momy_node + (dt * rhs_momy);
              momz_nodes_new[q] = momz_node + (dt * rhs_momz);
              energy_nodes_new[q] = std::max(energy_floor, energy_node + (dt * rhs_energy));
            }

            for (int n = 0; n < n_modes; ++n) {
              Real projected_momx = 0.0;
              Real projected_momy = 0.0;
              Real projected_momz = 0.0;
              Real projected_momw = 0.0;
              Real projected_energy = 0.0;
              for (int q = 0; q < tables.NumQuadrature(); ++q) {
                projected_momx += weights[q] * tables.Phi(n, q) * momx_nodes_new[q];
                projected_momy += weights[q] * tables.Phi(n, q) * momy_nodes_new[q];
                projected_momz += weights[q] * tables.Phi(n, q) * momz_nodes_new[q];
                projected_momw += weights[q] * tables.Phi(n, q) * momw_nodes[q];
                projected_energy += weights[q] * tables.Phi(n, q) * energy_nodes_new[q];
              }
              momx_modes_new[n] = projected_momx;
              momy_modes_new[n] = projected_momy;
              momz_modes_new[n] = projected_momz;
              momw_modes_new[n] = projected_momw;
              energy_modes_new[n] = projected_energy;
            }
            for (int n = 0; n < n_modes; ++n) {
              const int base = PlasmaIndex(s, n, 0, n_modes);
              plasma_new(base + kPlasmaMomX, k, j, i) = momx_modes_new[n];
              plasma_new(base + kPlasmaMomY, k, j, i) = momy_modes_new[n];
              plasma_new(base + kPlasmaMomZ, k, j, i) = momz_modes_new[n];
              plasma_new(base + kPlasmaMomW, k, j, i) = momw_modes_new[n];
              plasma_new(base + kPlasmaEnergy, k, j, i) = energy_modes_new[n];
            }
          }

          // Source-step mode continuity coupling update for each species:
          // d_t rho^(n) = -(sqrt(2(n+1))/lambda) * j_w^(n+1)
          // Brane transport divergence is handled through the conservative flux path.
          for (int s = 0; s < kSpeciesCount; ++s) {
            for (int n = 0; n < n_modes; ++n) {
              const int base = PlasmaIndex(s, n, 0, n_modes);
              Real leak_rhs = 0.0;
              if (n + 1 < n_modes) {
                const int base_np1 = PlasmaIndex(s, n + 1, 0, n_modes);
                const Real coupling = std::sqrt(2.0 * static_cast<Real>(n + 1)) / lambda;
                leak_rhs = -coupling * plasma_new(base_np1 + kPlasmaMomW, k, j, i);
              }

              const Real rhs_rho = leak_rhs;
              const Real rho_old = plasma(base + kPlasmaRho, k, j, i);
              plasma_new(base + kPlasmaRho, k, j, i) =
                  std::max(rho_floor, rho_old + (dt * rhs_rho));
            }
          }

          for (int n = 0; n < n_modes; ++n) {
            const int ion_base = PlasmaIndex(0, n, 0, n_modes);
            const int ele_base = PlasmaIndex(1, n, 0, n_modes);

            const Real rho_ion = plasma_new(ion_base + kPlasmaRho, k, j, i);
            const Real rho_electron = plasma_new(ele_base + kPlasmaRho, k, j, i);
            const Real momx_ion = plasma_new(ion_base + kPlasmaMomX, k, j, i);
            const Real momx_electron = plasma_new(ele_base + kPlasmaMomX, k, j, i);
            const Real momy_ion = plasma_new(ion_base + kPlasmaMomY, k, j, i);
            const Real momy_electron = plasma_new(ele_base + kPlasmaMomY, k, j, i);
            const Real momz_ion = plasma_new(ion_base + kPlasmaMomZ, k, j, i);
            const Real momz_electron = plasma_new(ele_base + kPlasmaMomZ, k, j, i);
            const Real momw_ion = plasma_new(ion_base + kPlasmaMomW, k, j, i);
            const Real momw_electron = plasma_new(ele_base + kPlasmaMomW, k, j, i);

            j0_modes[n] = (qom_ion * rho_ion) + (qom_electron * rho_electron);
            jx_modes[n] = (qom_ion * momx_ion) + (qom_electron * momx_electron);
            jy_modes[n] = (qom_ion * momy_ion) + (qom_electron * momy_electron);
            jz_modes[n] = (qom_ion * momz_ion) + (qom_electron * momz_electron);
            jw_modes[n] = (qom_ion * momw_ion) + (qom_electron * momw_electron);
          }

          const Real cell_volume = coords.CellVolume(k, j, i);
          for (int n = 0; n < n_modes; ++n) {
            const Real coupling =
                (n + 1 < n_modes) ? (std::sqrt(2.0 * static_cast<Real>(n + 1)) / lambda)
                                  : 0.0;
            const Real jw_np1 = (n + 1 < n_modes) ? jw_modes[n + 1] : 0.0;
            const Real continuity_residual =
                ((j0_modes[n] - charge_modes_old[n]) / dt) + (coupling * jw_np1);
            const Real abs_residual = std::abs(continuity_residual);

            diag_cont_local_l1_step += cell_volume * abs_residual;
            diag_cont_local_l2_step +=
                cell_volume * continuity_residual * continuity_residual;
            diag_cont_local_max_abs_step =
                std::max(diag_cont_local_max_abs_step, abs_residual);

            if (n == 0) {
              diag_cont_mode0_l1_step += cell_volume * abs_residual;
              diag_cont_mode0_l2_step +=
                  cell_volume * continuity_residual * continuity_residual;
              diag_cont_mode0_max_abs_step =
                  std::max(diag_cont_mode0_max_abs_step, abs_residual);
            }
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
          diag_jw_ew_step += dt * cell_volume * jw_ew_density;
          diag_s_leak_step += dt * cell_volume * s_leak_density;
          diag_s_leak_abs_step += dt * cell_volume * std::abs(s_leak_density);
        }
      }
    }

    a_dev.DeepCopy(a_new);
    pi_dev.DeepCopy(pi_new);
    plasma_dev.DeepCopy(plasma_new);
  }

  auto *diag_jw_ew = modes_pkg->MutableParam<double>("diag/int_jw_ew");
  auto *diag_s_leak = modes_pkg->MutableParam<double>("diag/int_s_leak");
  auto *diag_s_leak_abs = modes_pkg->MutableParam<double>("diag/int_s_leak_abs");
  auto *diag_cont_local_l1 = modes_pkg->MutableParam<double>("diag/continuity_local_l1");
  auto *diag_cont_local_l2 = modes_pkg->MutableParam<double>("diag/continuity_local_l2");
  auto *diag_cont_local_max_abs =
      modes_pkg->MutableParam<double>("diag/continuity_local_max_abs");
  auto *diag_cont_mode0_l1 = modes_pkg->MutableParam<double>("diag/continuity_mode0_l1");
  auto *diag_cont_mode0_l2 = modes_pkg->MutableParam<double>("diag/continuity_mode0_l2");
  auto *diag_cont_mode0_max_abs =
      modes_pkg->MutableParam<double>("diag/continuity_mode0_max_abs");
  *diag_jw_ew += diag_jw_ew_step;
  *diag_s_leak += diag_s_leak_step;
  *diag_s_leak_abs += diag_s_leak_abs_step;
  *diag_cont_local_l1 = diag_cont_local_l1_step;
  *diag_cont_local_l2 = diag_cont_local_l2_step;
  *diag_cont_local_max_abs = diag_cont_local_max_abs_step;
  *diag_cont_mode0_l1 = diag_cont_mode0_l1_step;
  *diag_cont_mode0_l2 = diag_cont_mode0_l2_step;
  *diag_cont_mode0_max_abs = diag_cont_mode0_max_abs_step;
}

} // namespace Modes4D
