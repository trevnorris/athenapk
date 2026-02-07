#include "diagnostics4d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "mode_tables.hpp"
#include "outputs/outputs.hpp"

namespace Modes4D {
using namespace parthenon::package::prelude;

namespace {

constexpr int kCompA0 = 0;
constexpr int kCompAX = 1;
constexpr int kCompAY = 2;
constexpr int kCompAZ = 3;
constexpr int kCompAW = 4;
constexpr int kVarsPerModePerSpecies = 6;
constexpr int kPlasmaRho = 0;
constexpr int kPlasmaMomX = 1;
constexpr int kPlasmaMomY = 2;
constexpr int kPlasmaMomZ = 3;
constexpr int kPlasmaMomW = 4;
constexpr int kPlasmaEnergy = 5;

int PlasmaIndex(const int species, const int mode, const int var, const int n_modes) {
  return (species * kVarsPerModePerSpecies * n_modes) + (mode * kVarsPerModePerSpecies) +
         var;
}

enum class BraneMixedQuantity {
  BraneE2,
  BraneB2,
  BraneEParallel2,
  MixedEw2,
  MixedC2
};

enum class HelicitySubscaleQuantity {
  HSub,
  EdotBSub
};

enum class TransverseEnergyQuantity {
  PoyntingW,
  LeakageW
};

Real FieldL2Integral(MeshData<Real> *md, const std::string &field_name) {
  auto *pmb = md->GetBlockData(0)->GetBlockPointer();
  const auto &field_pack = md->PackVariables(std::vector<std::string>{field_name});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real l2 = 0.0;
  pmb->par_reduce(
      "Modes4DFieldL2Integral", 0, field_pack.GetDim(5) - 1, 0, field_pack.GetDim(4) - 1,
      kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i,
                    Real &local_sum) {
        const auto &field = field_pack(b);
        const auto &coords = field_pack.GetCoords(b);
        const Real val = field(v, k, j, i);
        local_sum += val * val * coords.CellVolume(k, j, i);
      },
      l2);
  return 0.5 * l2;
}

Real EMA2Hst(MeshData<Real> *md) { return FieldL2Integral(md, "em4d_a"); }

Real EMPi2Hst(MeshData<Real> *md) { return FieldL2Integral(md, "em4d_pi"); }

Real PlasmaCons2Hst(MeshData<Real> *md) { return FieldL2Integral(md, "plasma4d_cons"); }

Real JwEwAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_jw_ew");
}

Real JaEaAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_ja_ea");
}

Real LeakAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_s_leak");
}

Real LeakAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_s_leak_abs");
}

Real DivJMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divj_mode0");
}

Real DivMomXMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divmomx_mode0");
}

Real DivMomYMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divmomy_mode0");
}

Real DivMomZMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divmomz_mode0");
}

Real DivMomWMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divmomw_mode0");
}

Real DivEnergyMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divenergy_mode0");
}

Real DivPi0Mode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divpi0_mode0");
}

Real DivPixMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divpix_mode0");
}

Real DivPiyMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divpiy_mode0");
}

Real DivPizMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divpiz_mode0");
}

Real DivPiwMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_divpiw_mode0");
}

Real SrcMomXMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcmomx_mode0");
}

Real SrcMomYMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcmomy_mode0");
}

Real SrcMomZMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcmomz_mode0");
}

Real SrcMomWMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcmomw_mode0");
}

Real SrcEnergyMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcenergy_mode0");
}

Real SrcPi0Mode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcpi0_mode0");
}

Real SrcPixMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcpix_mode0");
}

Real SrcPiyMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcpiy_mode0");
}

Real SrcPizMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcpiz_mode0");
}

Real SrcPiwMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_srcpiw_mode0");
}

Real SrcTimelikeA0FromPiwAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_timelike_a0_from_piw");
}

Real SrcTimelikeA0FromPiwAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_timelike_a0_from_piw_abs");
}

Real SrcTimelikeAwFromPi0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_timelike_aw_from_pi0");
}

Real SrcTimelikeAwFromPi0AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_timelike_aw_from_pi0_abs");
}

Real SrcEMLaplacianAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_laplacian_abs");
}

Real SrcEMMassAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_mass_abs");
}

Real SrcEMCurrentAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_current_abs");
}

Real SrcEMDampingAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_damping_abs");
}

Real SrcEMSpatialMixedAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_spatial_mixed_abs");
}

Real SrcEMTimelikeAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_timelike_abs");
}

Real ContinuityLocalL1Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/continuity_local_l1");
}

Real ContinuityLocalL2Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/continuity_local_l2");
}

Real ContinuityLocalMaxAbsHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/continuity_local_max_abs");
}

Real ContinuityMode0L1Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/continuity_mode0_l1");
}

Real ContinuityMode0L2Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/continuity_mode0_l2");
}

Real ContinuityMode0MaxAbsHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/continuity_mode0_max_abs");
}

Real PulseCentroidXHst(MeshData<Real> *md) {
  auto *pmb = md->GetBlockData(0)->GetBlockPointer();
  const auto &field_pack = md->PackVariables(std::vector<std::string>{"em4d_a"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real weighted_x = 0.0;
  pmb->par_reduce(
      "Modes4DPulseCentroidX", 0, field_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e,
      ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &wx) {
        const auto &field = field_pack(b);
        const auto &coords = field_pack.GetCoords(b);
        const Real ay0 = field(2, k, j, i);
        const Real local_weight = ay0 * ay0 * coords.CellVolume(k, j, i);
        wx += coords.Xc<1>(k, j, i) * local_weight;
      },
      weighted_x);

  Real weight = 0.0;
  pmb->par_reduce(
      "Modes4DPulseWeight", 0, field_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s,
      ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &w) {
        const auto &field = field_pack(b);
        const auto &coords = field_pack.GetCoords(b);
        const Real ay0 = field(2, k, j, i);
        w += ay0 * ay0 * coords.CellVolume(k, j, i);
      },
      weight);

  return (weight > 0.0) ? (weighted_x / weight) : 0.0;
}

Real BraneMixedIntegral(MeshData<Real> *md, BraneMixedQuantity quantity) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();

  Real integral = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();
    auto &coords = pmb->coords;

    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);
    const bool has_x = (ib.e > ib.s);
    const bool has_y = (jb.e > jb.s);
    const bool has_z = (kb.e > kb.s);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real dx = has_x ? coords.Dxc<1>(i) : 1.0;
          const Real dy = has_y ? coords.Dxc<2>(j) : 1.0;
          const Real dz = has_z ? coords.Dxc<3>(k) : 1.0;

          const Real dA0_dx =
              has_x ? (a(kCompA0, k, j, i + 1) - a(kCompA0, k, j, i - 1)) / (2.0 * dx) : 0.0;
          const Real dA0_dy =
              has_y ? (a(kCompA0, k, j + 1, i) - a(kCompA0, k, j - 1, i)) / (2.0 * dy) : 0.0;
          const Real dA0_dz =
              has_z ? (a(kCompA0, k + 1, j, i) - a(kCompA0, k - 1, j, i)) / (2.0 * dz) : 0.0;

          const Real dAx_dy =
              has_y ? (a(kCompAX, k, j + 1, i) - a(kCompAX, k, j - 1, i)) / (2.0 * dy) : 0.0;
          const Real dAx_dz =
              has_z ? (a(kCompAX, k + 1, j, i) - a(kCompAX, k - 1, j, i)) / (2.0 * dz) : 0.0;
          const Real dAy_dx =
              has_x ? (a(kCompAY, k, j, i + 1) - a(kCompAY, k, j, i - 1)) / (2.0 * dx) : 0.0;
          const Real dAy_dz =
              has_z ? (a(kCompAY, k + 1, j, i) - a(kCompAY, k - 1, j, i)) / (2.0 * dz) : 0.0;
          const Real dAz_dx =
              has_x ? (a(kCompAZ, k, j, i + 1) - a(kCompAZ, k, j, i - 1)) / (2.0 * dx) : 0.0;
          const Real dAz_dy =
              has_y ? (a(kCompAZ, k, j + 1, i) - a(kCompAZ, k, j - 1, i)) / (2.0 * dy) : 0.0;

          const Real ex = -pi(kCompAX, k, j, i) - dA0_dx;
          const Real ey = -pi(kCompAY, k, j, i) - dA0_dy;
          const Real ez = -pi(kCompAZ, k, j, i) - dA0_dz;

          const Real bx = dAz_dy - dAy_dz;
          const Real by = dAx_dz - dAz_dx;
          const Real bz = dAy_dx - dAx_dy;

          const Real b2 = (bx * bx) + (by * by) + (bz * bz);
          const Real e2 = (ex * ex) + (ey * ey) + (ez * ez);
          const Real edotb = (ex * bx) + (ey * by) + (ez * bz);
          const Real epar2 = (edotb * edotb) / (b2 + 1.0e-30);

          std::vector<Real> a0_modes(n_modes, 0.0);
          std::vector<Real> ax_modes(n_modes, 0.0);
          std::vector<Real> ay_modes(n_modes, 0.0);
          std::vector<Real> az_modes(n_modes, 0.0);
          std::vector<Real> aw_modes(n_modes, 0.0);
          std::vector<Real> piw_modes(n_modes, 0.0);
          std::vector<Real> daw_dx_modes(n_modes, 0.0);
          std::vector<Real> daw_dy_modes(n_modes, 0.0);
          std::vector<Real> daw_dz_modes(n_modes, 0.0);

          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            a0_modes[n] = a(off + kCompA0, k, j, i);
            ax_modes[n] = a(off + kCompAX, k, j, i);
            ay_modes[n] = a(off + kCompAY, k, j, i);
            az_modes[n] = a(off + kCompAZ, k, j, i);
            aw_modes[n] = a(off + kCompAW, k, j, i);
            piw_modes[n] = pi(off + kCompAW, k, j, i);

            if (has_x) {
              daw_dx_modes[n] =
                  (a(off + kCompAW, k, j, i + 1) - a(off + kCompAW, k, j, i - 1)) /
                  (2.0 * dx);
            }
            if (has_y) {
              daw_dy_modes[n] =
                  (a(off + kCompAW, k, j + 1, i) - a(off + kCompAW, k, j - 1, i)) /
                  (2.0 * dy);
            }
            if (has_z) {
              daw_dz_modes[n] =
                  (a(off + kCompAW, k + 1, j, i) - a(off + kCompAW, k - 1, j, i)) /
                  (2.0 * dz);
            }
          }

          Real ew2_int = 0.0;
          Real c2_int = 0.0;
          for (int q = 0; q < n_quad; ++q) {
            Real d_w_a0 = 0.0;
            Real d_w_ax = 0.0;
            Real d_w_ay = 0.0;
            Real d_w_az = 0.0;
            Real piw_node = 0.0;
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
              daw_dx_node += daw_dx_modes[n] * phi_nq;
              daw_dy_node += daw_dy_modes[n] * phi_nq;
              daw_dz_node += daw_dz_modes[n] * phi_nq;
            }

            const Real ew = -piw_node - d_w_a0;
            const Real cx = daw_dx_node - d_w_ax;
            const Real cy = daw_dy_node - d_w_ay;
            const Real cz = daw_dz_node - d_w_az;

            ew2_int += weights[q] * (ew * ew);
            c2_int += weights[q] * ((cx * cx) + (cy * cy) + (cz * cz));
          }

          const Real vol = coords.CellVolume(k, j, i);
          Real val = 0.0;
          if (quantity == BraneMixedQuantity::BraneE2) {
            val = 0.5 * e2;
          } else if (quantity == BraneMixedQuantity::BraneB2) {
            val = 0.5 * b2;
          } else if (quantity == BraneMixedQuantity::BraneEParallel2) {
            val = epar2;
          } else if (quantity == BraneMixedQuantity::MixedEw2) {
            val = 0.5 * ew2_int;
          } else {
            val = 0.5 * c2_int;
          }
          integral += val * vol;
        }
      }
    }
  }
  return integral;
}

Real BulkEMEnergyIntegral(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const Real mu0 = modes_pkg->Param<double>("em4d/mu0");

  Real integral = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();
    auto &coords = pmb->coords;

    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);
    const bool has_x = (ib.e > ib.s);
    const bool has_y = (jb.e > jb.s);
    const bool has_z = (kb.e > kb.s);

    std::vector<Real> a0_modes(n_modes, 0.0);
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

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real dx = has_x ? coords.Dxc<1>(i) : 1.0;
          const Real dy = has_y ? coords.Dxc<2>(j) : 1.0;
          const Real dz = has_z ? coords.Dxc<3>(k) : 1.0;

          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            a0_modes[n] = a(off + kCompA0, k, j, i);
            ax_modes[n] = a(off + kCompAX, k, j, i);
            ay_modes[n] = a(off + kCompAY, k, j, i);
            az_modes[n] = a(off + kCompAZ, k, j, i);
            piw_modes[n] = pi(off + kCompAW, k, j, i);
            pix_modes[n] = pi(off + kCompAX, k, j, i);
            piy_modes[n] = pi(off + kCompAY, k, j, i);
            piz_modes[n] = pi(off + kCompAZ, k, j, i);

            if (has_x) {
              da0_dx_modes[n] =
                  (a(off + kCompA0, k, j, i + 1) - a(off + kCompA0, k, j, i - 1)) /
                  (2.0 * dx);
              day_dx_modes[n] =
                  (a(off + kCompAY, k, j, i + 1) - a(off + kCompAY, k, j, i - 1)) /
                  (2.0 * dx);
              daz_dx_modes[n] =
                  (a(off + kCompAZ, k, j, i + 1) - a(off + kCompAZ, k, j, i - 1)) /
                  (2.0 * dx);
              daw_dx_modes[n] =
                  (a(off + kCompAW, k, j, i + 1) - a(off + kCompAW, k, j, i - 1)) /
                  (2.0 * dx);
            }
            if (has_y) {
              da0_dy_modes[n] =
                  (a(off + kCompA0, k, j + 1, i) - a(off + kCompA0, k, j - 1, i)) /
                  (2.0 * dy);
              dax_dy_modes[n] =
                  (a(off + kCompAX, k, j + 1, i) - a(off + kCompAX, k, j - 1, i)) /
                  (2.0 * dy);
              daz_dy_modes[n] =
                  (a(off + kCompAZ, k, j + 1, i) - a(off + kCompAZ, k, j - 1, i)) /
                  (2.0 * dy);
              daw_dy_modes[n] =
                  (a(off + kCompAW, k, j + 1, i) - a(off + kCompAW, k, j - 1, i)) /
                  (2.0 * dy);
            }
            if (has_z) {
              da0_dz_modes[n] =
                  (a(off + kCompA0, k + 1, j, i) - a(off + kCompA0, k - 1, j, i)) /
                  (2.0 * dz);
              dax_dz_modes[n] =
                  (a(off + kCompAX, k + 1, j, i) - a(off + kCompAX, k - 1, j, i)) /
                  (2.0 * dz);
              day_dz_modes[n] =
                  (a(off + kCompAY, k + 1, j, i) - a(off + kCompAY, k - 1, j, i)) /
                  (2.0 * dz);
              daw_dz_modes[n] =
                  (a(off + kCompAW, k + 1, j, i) - a(off + kCompAW, k - 1, j, i)) /
                  (2.0 * dz);
            }
          }

          Real local_w_integral = 0.0;
          for (int q = 0; q < n_quad; ++q) {
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

            const Real ex = -pix_node - da0_dx_node;
            const Real ey = -piy_node - da0_dy_node;
            const Real ez = -piz_node - da0_dz_node;
            const Real bx = daz_dy_node - day_dz_node;
            const Real by = dax_dz_node - daz_dx_node;
            const Real bz = day_dx_node - dax_dy_node;
            const Real ew = -piw_node - d_w_a0;
            const Real cx = daw_dx_node - d_w_ax;
            const Real cy = daw_dy_node - d_w_ay;
            const Real cz = daw_dz_node - d_w_az;

            local_w_integral +=
                weights[q] *
                ((ex * ex) + (ey * ey) + (ez * ez) + (ew * ew) + (bx * bx) +
                 (by * by) + (bz * bz) + (cx * cx) + (cy * cy) + (cz * cz));
          }

          integral += (0.5 / mu0) * local_w_integral * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return integral;
}

Real BulkEMEnergyHst(MeshData<Real> *md) { return BulkEMEnergyIntegral(md); }

Real ResolvedEMEnergyIntegral(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const Real mu0 = modes_pkg->Param<double>("em4d/mu0");
  const Real z_int = tables.Lambda() * std::sqrt(std::acos(-1.0));
  const Real inv_z_int = 1.0 / z_int;

  Real integral = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();
    auto &coords = pmb->coords;

    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);
    const bool has_x = (ib.e > ib.s);
    const bool has_y = (jb.e > jb.s);
    const bool has_z = (kb.e > kb.s);

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

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real dx = has_x ? coords.Dxc<1>(i) : 1.0;
          const Real dy = has_y ? coords.Dxc<2>(j) : 1.0;
          const Real dz = has_z ? coords.Dxc<3>(k) : 1.0;

          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            pix_modes[n] = pi(off + kCompAX, k, j, i);
            piy_modes[n] = pi(off + kCompAY, k, j, i);
            piz_modes[n] = pi(off + kCompAZ, k, j, i);

            if (has_x) {
              da0_dx_modes[n] =
                  (a(off + kCompA0, k, j, i + 1) - a(off + kCompA0, k, j, i - 1)) /
                  (2.0 * dx);
              day_dx_modes[n] =
                  (a(off + kCompAY, k, j, i + 1) - a(off + kCompAY, k, j, i - 1)) /
                  (2.0 * dx);
              daz_dx_modes[n] =
                  (a(off + kCompAZ, k, j, i + 1) - a(off + kCompAZ, k, j, i - 1)) /
                  (2.0 * dx);
            }
            if (has_y) {
              da0_dy_modes[n] =
                  (a(off + kCompA0, k, j + 1, i) - a(off + kCompA0, k, j - 1, i)) /
                  (2.0 * dy);
              dax_dy_modes[n] =
                  (a(off + kCompAX, k, j + 1, i) - a(off + kCompAX, k, j - 1, i)) /
                  (2.0 * dy);
              daz_dy_modes[n] =
                  (a(off + kCompAZ, k, j + 1, i) - a(off + kCompAZ, k, j - 1, i)) /
                  (2.0 * dy);
            }
            if (has_z) {
              da0_dz_modes[n] =
                  (a(off + kCompA0, k + 1, j, i) - a(off + kCompA0, k - 1, j, i)) /
                  (2.0 * dz);
              dax_dz_modes[n] =
                  (a(off + kCompAX, k + 1, j, i) - a(off + kCompAX, k - 1, j, i)) /
                  (2.0 * dz);
              day_dz_modes[n] =
                  (a(off + kCompAY, k + 1, j, i) - a(off + kCompAY, k - 1, j, i)) /
                  (2.0 * dz);
            }
          }

          Real ex_bar = 0.0;
          Real ey_bar = 0.0;
          Real ez_bar = 0.0;
          Real bx_bar = 0.0;
          Real by_bar = 0.0;
          Real bz_bar = 0.0;

          for (int q = 0; q < n_quad; ++q) {
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

            for (int n = 0; n < n_modes; ++n) {
              const Real phi_nq = tables.Phi(n, q);
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
            }

            const Real ex = -pix_node - da0_dx_node;
            const Real ey = -piy_node - da0_dy_node;
            const Real ez = -piz_node - da0_dz_node;
            const Real bx = daz_dy_node - day_dz_node;
            const Real by = dax_dz_node - daz_dx_node;
            const Real bz = day_dx_node - dax_dy_node;

            ex_bar += weights[q] * ex;
            ey_bar += weights[q] * ey;
            ez_bar += weights[q] * ez;
            bx_bar += weights[q] * bx;
            by_bar += weights[q] * by;
            bz_bar += weights[q] * bz;
          }

          ex_bar *= inv_z_int;
          ey_bar *= inv_z_int;
          ez_bar *= inv_z_int;
          bx_bar *= inv_z_int;
          by_bar *= inv_z_int;
          bz_bar *= inv_z_int;

          const Real ebar2 = (ex_bar * ex_bar) + (ey_bar * ey_bar) + (ez_bar * ez_bar);
          const Real bbar2 = (bx_bar * bx_bar) + (by_bar * by_bar) + (bz_bar * bz_bar);
          const Real u_resolved = 0.5 * (z_int / mu0) * (ebar2 + bbar2);

          integral += u_resolved * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return integral;
}

Real ResolvedEMEnergyHst(MeshData<Real> *md) { return ResolvedEMEnergyIntegral(md); }

Real HelicitySubscaleIntegral(MeshData<Real> *md, HelicitySubscaleQuantity quantity) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const Real z_int = tables.Lambda() * std::sqrt(std::acos(-1.0));
  const Real inv_z_int = 1.0 / z_int;

  Real integral = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();
    auto &coords = pmb->coords;

    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);
    const bool has_x = (ib.e > ib.s);
    const bool has_y = (jb.e > jb.s);
    const bool has_z = (kb.e > kb.s);

    std::vector<Real> ax_modes(n_modes, 0.0);
    std::vector<Real> ay_modes(n_modes, 0.0);
    std::vector<Real> az_modes(n_modes, 0.0);
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

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real dx = has_x ? coords.Dxc<1>(i) : 1.0;
          const Real dy = has_y ? coords.Dxc<2>(j) : 1.0;
          const Real dz = has_z ? coords.Dxc<3>(k) : 1.0;

          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            ax_modes[n] = a(off + kCompAX, k, j, i);
            ay_modes[n] = a(off + kCompAY, k, j, i);
            az_modes[n] = a(off + kCompAZ, k, j, i);
            pix_modes[n] = pi(off + kCompAX, k, j, i);
            piy_modes[n] = pi(off + kCompAY, k, j, i);
            piz_modes[n] = pi(off + kCompAZ, k, j, i);

            if (has_x) {
              da0_dx_modes[n] =
                  (a(off + kCompA0, k, j, i + 1) - a(off + kCompA0, k, j, i - 1)) /
                  (2.0 * dx);
              day_dx_modes[n] =
                  (a(off + kCompAY, k, j, i + 1) - a(off + kCompAY, k, j, i - 1)) /
                  (2.0 * dx);
              daz_dx_modes[n] =
                  (a(off + kCompAZ, k, j, i + 1) - a(off + kCompAZ, k, j, i - 1)) /
                  (2.0 * dx);
            }
            if (has_y) {
              da0_dy_modes[n] =
                  (a(off + kCompA0, k, j + 1, i) - a(off + kCompA0, k, j - 1, i)) /
                  (2.0 * dy);
              dax_dy_modes[n] =
                  (a(off + kCompAX, k, j + 1, i) - a(off + kCompAX, k, j - 1, i)) /
                  (2.0 * dy);
              daz_dy_modes[n] =
                  (a(off + kCompAZ, k, j + 1, i) - a(off + kCompAZ, k, j - 1, i)) /
                  (2.0 * dy);
            }
            if (has_z) {
              da0_dz_modes[n] =
                  (a(off + kCompA0, k + 1, j, i) - a(off + kCompA0, k - 1, j, i)) /
                  (2.0 * dz);
              dax_dz_modes[n] =
                  (a(off + kCompAX, k + 1, j, i) - a(off + kCompAX, k - 1, j, i)) /
                  (2.0 * dz);
              day_dz_modes[n] =
                  (a(off + kCompAY, k + 1, j, i) - a(off + kCompAY, k - 1, j, i)) /
                  (2.0 * dz);
            }
          }

          Real ax_bar = 0.0;
          Real ay_bar = 0.0;
          Real az_bar = 0.0;
          Real ex_bar = 0.0;
          Real ey_bar = 0.0;
          Real ez_bar = 0.0;
          Real bx_bar = 0.0;
          Real by_bar = 0.0;
          Real bz_bar = 0.0;
          Real ab_proj = 0.0;
          Real edotb_proj = 0.0;

          for (int q = 0; q < n_quad; ++q) {
            Real ax_node = 0.0;
            Real ay_node = 0.0;
            Real az_node = 0.0;
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

            for (int n = 0; n < n_modes; ++n) {
              const Real phi_nq = tables.Phi(n, q);
              ax_node += ax_modes[n] * phi_nq;
              ay_node += ay_modes[n] * phi_nq;
              az_node += az_modes[n] * phi_nq;
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
            }

            const Real ex = -pix_node - da0_dx_node;
            const Real ey = -piy_node - da0_dy_node;
            const Real ez = -piz_node - da0_dz_node;
            const Real bx = daz_dy_node - day_dz_node;
            const Real by = dax_dz_node - daz_dx_node;
            const Real bz = day_dx_node - dax_dy_node;
            const Real ab = (ax_node * bx) + (ay_node * by) + (az_node * bz);
            const Real edotb = (ex * bx) + (ey * by) + (ez * bz);

            ax_bar += weights[q] * ax_node;
            ay_bar += weights[q] * ay_node;
            az_bar += weights[q] * az_node;
            ex_bar += weights[q] * ex;
            ey_bar += weights[q] * ey;
            ez_bar += weights[q] * ez;
            bx_bar += weights[q] * bx;
            by_bar += weights[q] * by;
            bz_bar += weights[q] * bz;
            ab_proj += weights[q] * ab;
            edotb_proj += weights[q] * edotb;
          }

          ax_bar *= inv_z_int;
          ay_bar *= inv_z_int;
          az_bar *= inv_z_int;
          ex_bar *= inv_z_int;
          ey_bar *= inv_z_int;
          ez_bar *= inv_z_int;
          bx_bar *= inv_z_int;
          by_bar *= inv_z_int;
          bz_bar *= inv_z_int;
          ab_proj *= inv_z_int;
          edotb_proj *= inv_z_int;

          const Real ab_res = (ax_bar * bx_bar) + (ay_bar * by_bar) + (az_bar * bz_bar);
          const Real edotb_res = (ex_bar * bx_bar) + (ey_bar * by_bar) + (ez_bar * bz_bar);

          if (quantity == HelicitySubscaleQuantity::HSub) {
            integral += (ab_proj - ab_res) * coords.CellVolume(k, j, i);
          } else {
            integral += (edotb_proj - edotb_res) * coords.CellVolume(k, j, i);
          }
        }
      }
    }
  }
  return integral;
}

Real HelicitySubHst(MeshData<Real> *md) {
  return HelicitySubscaleIntegral(md, HelicitySubscaleQuantity::HSub);
}

Real EdotBSubHst(MeshData<Real> *md) {
  return HelicitySubscaleIntegral(md, HelicitySubscaleQuantity::EdotBSub);
}

Real TransverseEnergyIntegral(MeshData<Real> *md, TransverseEnergyQuantity quantity) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const auto &nodes = tables.Nodes();
  const Real lambda = tables.Lambda();
  const Real z_int = lambda * std::sqrt(std::acos(-1.0));
  const Real mu0 = modes_pkg->Param<double>("em4d/mu0");

  Real integral = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();
    auto &coords = pmb->coords;

    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);
    const bool has_x = (ib.e > ib.s);
    const bool has_y = (jb.e > jb.s);
    const bool has_z = (kb.e > kb.s);

    std::vector<Real> a0_modes(n_modes, 0.0);
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
    std::vector<Real> daw_dx_modes(n_modes, 0.0);
    std::vector<Real> daw_dy_modes(n_modes, 0.0);
    std::vector<Real> daw_dz_modes(n_modes, 0.0);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real dx = has_x ? coords.Dxc<1>(i) : 1.0;
          const Real dy = has_y ? coords.Dxc<2>(j) : 1.0;
          const Real dz = has_z ? coords.Dxc<3>(k) : 1.0;

          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            a0_modes[n] = a(off + kCompA0, k, j, i);
            ax_modes[n] = a(off + kCompAX, k, j, i);
            ay_modes[n] = a(off + kCompAY, k, j, i);
            az_modes[n] = a(off + kCompAZ, k, j, i);
            piw_modes[n] = pi(off + kCompAW, k, j, i);
            pix_modes[n] = pi(off + kCompAX, k, j, i);
            piy_modes[n] = pi(off + kCompAY, k, j, i);
            piz_modes[n] = pi(off + kCompAZ, k, j, i);

            if (has_x) {
              da0_dx_modes[n] =
                  (a(off + kCompA0, k, j, i + 1) - a(off + kCompA0, k, j, i - 1)) /
                  (2.0 * dx);
              daw_dx_modes[n] =
                  (a(off + kCompAW, k, j, i + 1) - a(off + kCompAW, k, j, i - 1)) /
                  (2.0 * dx);
            }
            if (has_y) {
              da0_dy_modes[n] =
                  (a(off + kCompA0, k, j + 1, i) - a(off + kCompA0, k, j - 1, i)) /
                  (2.0 * dy);
              daw_dy_modes[n] =
                  (a(off + kCompAW, k, j + 1, i) - a(off + kCompAW, k, j - 1, i)) /
                  (2.0 * dy);
            }
            if (has_z) {
              da0_dz_modes[n] =
                  (a(off + kCompA0, k + 1, j, i) - a(off + kCompA0, k - 1, j, i)) /
                  (2.0 * dz);
              daw_dz_modes[n] =
                  (a(off + kCompAW, k + 1, j, i) - a(off + kCompAW, k - 1, j, i)) /
                  (2.0 * dz);
            }
          }

          Real local_w_integral = 0.0;
          for (int q = 0; q < n_quad; ++q) {
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
              daw_dx_node += daw_dx_modes[n] * phi_nq;
              daw_dy_node += daw_dy_modes[n] * phi_nq;
              daw_dz_node += daw_dz_modes[n] * phi_nq;
            }

            const Real ex = -pix_node - da0_dx_node;
            const Real ey = -piy_node - da0_dy_node;
            const Real ez = -piz_node - da0_dz_node;
            const Real cx = daw_dx_node - d_w_ax;
            const Real cy = daw_dy_node - d_w_ay;
            const Real cz = daw_dz_node - d_w_az;
            const Real edotc = (ex * cx) + (ey * cy) + (ez * cz);

            if (quantity == TransverseEnergyQuantity::PoyntingW) {
              local_w_integral += weights[q] * (edotc / mu0);
            } else {
              const Real w = nodes[q];
              const Real z = std::exp(-(w * w) / (lambda * lambda));
              const Real leakage_weight = (-2.0 * w) / (lambda * lambda * z_int);
              local_w_integral +=
                  weights[q] * leakage_weight * ((z / mu0) * edotc);
            }
          }

          integral += local_w_integral * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return integral;
}

Real EMPoyntingWHst(MeshData<Real> *md) {
  return TransverseEnergyIntegral(md, TransverseEnergyQuantity::PoyntingW);
}

Real EMLeakageWHst(MeshData<Real> *md) {
  return TransverseEnergyIntegral(md, TransverseEnergyQuantity::LeakageW);
}

Real BraneE2Hst(MeshData<Real> *md) {
  return BraneMixedIntegral(md, BraneMixedQuantity::BraneE2);
}

Real BraneB2Hst(MeshData<Real> *md) {
  return BraneMixedIntegral(md, BraneMixedQuantity::BraneB2);
}

Real BraneEParallel2Hst(MeshData<Real> *md) {
  return BraneMixedIntegral(md, BraneMixedQuantity::BraneEParallel2);
}

Real MixedEw2Hst(MeshData<Real> *md) {
  return BraneMixedIntegral(md, BraneMixedQuantity::MixedEw2);
}

Real MixedC2Hst(MeshData<Real> *md) {
  return BraneMixedIntegral(md, BraneMixedQuantity::MixedC2);
}

Real Ay0SpanHst(MeshData<Real> *md) {
  Real ay0_min = std::numeric_limits<Real>::max();
  Real ay0_max = std::numeric_limits<Real>::lowest();

  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real ay0 = a(kCompAY, k, j, i);
          ay0_min = std::min(ay0_min, ay0);
          ay0_max = std::max(ay0_max, ay0);
        }
      }
    }
  }

  return (ay0_max > ay0_min) ? (ay0_max - ay0_min) : 0.0;
}

Real AyProjectedSpanHst(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const Real z_int = tables.Lambda() * std::sqrt(std::acos(-1.0));

  std::vector<Real> proj_coeff(static_cast<size_t>(n_modes), 0.0);
  for (int n = 0; n < n_modes; ++n) {
    Real coeff = 0.0;
    for (int q = 0; q < n_quad; ++q) {
      coeff += weights[q] * tables.Phi(n, q);
    }
    proj_coeff[static_cast<size_t>(n)] = coeff / z_int;
  }

  Real ay_proj_min = std::numeric_limits<Real>::max();
  Real ay_proj_max = std::numeric_limits<Real>::lowest();

  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          Real ay_proj = 0.0;
          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            ay_proj += proj_coeff[static_cast<size_t>(n)] * a(off + kCompAY, k, j, i);
          }
          ay_proj_min = std::min(ay_proj_min, ay_proj);
          ay_proj_max = std::max(ay_proj_max, ay_proj);
        }
      }
    }
  }

  return (ay_proj_max > ay_proj_min) ? (ay_proj_max - ay_proj_min) : 0.0;
}

Real JwModeL2Integral(MeshData<Real> *md, const int mode_idx) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");
  if (mode_idx < 0 || mode_idx >= n_modes) {
    return 0.0;
  }

  Real sum = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto plasma = bd->Get("plasma4d_cons").data.GetHostMirrorAndCopy();
    auto &coords = bd->GetBlockPointer()->coords;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const int ion_base = PlasmaIndex(0, mode_idx, 0, n_modes);
          const int ele_base = PlasmaIndex(1, mode_idx, 0, n_modes);
          const Real ion_momw = plasma(ion_base + kPlasmaMomW, k, j, i);
          const Real ele_momw = plasma(ele_base + kPlasmaMomW, k, j, i);
          const Real jw = (qom_ion * ion_momw) + (qom_electron * ele_momw);
          sum += 0.5 * jw * jw * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real JwModeIntegral(MeshData<Real> *md, const int mode_idx) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");
  if (mode_idx < 0 || mode_idx >= n_modes) {
    return 0.0;
  }

  Real sum = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto plasma = bd->Get("plasma4d_cons").data.GetHostMirrorAndCopy();
    auto &coords = bd->GetBlockPointer()->coords;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const int ion_base = PlasmaIndex(0, mode_idx, 0, n_modes);
          const int ele_base = PlasmaIndex(1, mode_idx, 0, n_modes);
          const Real ion_momw = plasma(ion_base + kPlasmaMomW, k, j, i);
          const Real ele_momw = plasma(ele_base + kPlasmaMomW, k, j, i);
          const Real jw = (qom_ion * ion_momw) + (qom_electron * ele_momw);
          sum += jw * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real ChargeModeIntegral(MeshData<Real> *md, const int mode_idx) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");
  if (mode_idx < 0 || mode_idx >= n_modes) {
    return 0.0;
  }

  Real sum = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto plasma = bd->Get("plasma4d_cons").data.GetHostMirrorAndCopy();
    auto &coords = bd->GetBlockPointer()->coords;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const int ion_base = PlasmaIndex(0, mode_idx, 0, n_modes);
          const int ele_base = PlasmaIndex(1, mode_idx, 0, n_modes);
          const Real ion_rho = plasma(ion_base + kPlasmaRho, k, j, i);
          const Real ele_rho = plasma(ele_base + kPlasmaRho, k, j, i);
          const Real charge = (qom_ion * ion_rho) + (qom_electron * ele_rho);
          sum += charge * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real SpeciesSummedModeIntegral(MeshData<Real> *md, const int mode_idx, const int var_idx) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  if (mode_idx < 0 || mode_idx >= n_modes) {
    return 0.0;
  }

  Real sum = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto plasma = bd->Get("plasma4d_cons").data.GetHostMirrorAndCopy();
    auto &coords = bd->GetBlockPointer()->coords;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const int ion_base = PlasmaIndex(0, mode_idx, 0, n_modes);
          const int ele_base = PlasmaIndex(1, mode_idx, 0, n_modes);
          const Real ion_val = plasma(ion_base + var_idx, k, j, i);
          const Real ele_val = plasma(ele_base + var_idx, k, j, i);
          sum += (ion_val + ele_val) * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real EMModeComponentIntegral(MeshData<Real> *md, const int mode_idx, const int comp_idx) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  if (mode_idx < 0 || mode_idx >= n_modes || comp_idx < 0 || comp_idx >= 5) {
    return 0.0;
  }

  const int idx = (5 * mode_idx) + comp_idx;
  Real sum = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();
    auto &coords = bd->GetBlockPointer()->coords;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          sum += pi(idx, k, j, i) * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real FieldModeL2Integral(MeshData<Real> *md, const std::string &field_name,
                         const int mode_idx) {
  auto *pmb = md->GetBlockData(0)->GetBlockPointer();
  const auto &field_pack = md->PackVariables(std::vector<std::string>{field_name});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real sum = 0.0;
  pmb->par_reduce(
      "Modes4DFieldModeL2", 0, field_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s,
      ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i,
                    Real &local_sum) {
        const auto &field = field_pack(b);
        const auto &coords = field_pack.GetCoords(b);
        Real val2 = 0.0;
        for (int c = 0; c < 5; ++c) {
          const Real val = field(5 * mode_idx + c, k, j, i);
          val2 += val * val;
        }
        local_sum += 0.5 * val2 * coords.CellVolume(k, j, i);
      },
      sum);
  return sum;
}

} // namespace

void RegisterDiagnostics(parthenon::StateDescriptor *pkg) {
  // Placeholder accumulators for the required energy/leakage ledger channels.
  pkg->AddParam<double>("diag/int_jw_ew", 0.0, true);
  pkg->AddParam<double>("diag/int_ja_ea", 0.0, true);
  pkg->AddParam<double>("diag/int_s_leak", 0.0, true);
  pkg->AddParam<double>("diag/int_s_leak_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_divj_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divmomx_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divmomy_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divmomz_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divmomw_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divenergy_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divpi0_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divpix_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divpiy_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divpiz_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_divpiw_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcmomx_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcmomy_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcmomz_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcmomw_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcenergy_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcpi0_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcpix_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcpiy_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcpiz_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_srcpiw_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_src_timelike_a0_from_piw", 0.0, true);
  pkg->AddParam<double>("diag/int_src_timelike_a0_from_piw_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_timelike_aw_from_pi0", 0.0, true);
  pkg->AddParam<double>("diag/int_src_timelike_aw_from_pi0_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_laplacian_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_mass_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_current_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_damping_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_spatial_mixed_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_timelike_abs", 0.0, true);
  pkg->AddParam<double>("diag/continuity_local_l1", 0.0, true);
  pkg->AddParam<double>("diag/continuity_local_l2", 0.0, true);
  pkg->AddParam<double>("diag/continuity_local_max_abs", 0.0, true);
  pkg->AddParam<double>("diag/continuity_mode0_l1", 0.0, true);
  pkg->AddParam<double>("diag/continuity_mode0_l2", 0.0, true);
  pkg->AddParam<double>("diag/continuity_mode0_max_abs", 0.0, true);

  parthenon::HstVar_list hst_vars = {};
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMA2Hst, "m4d_em_a2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMPi2Hst, "m4d_em_pi2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    PlasmaCons2Hst,
                                                    "m4d_plasma_cons2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    JwEwAccumulatorHst,
                                                    "m4d_int_jw_ew"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    JaEaAccumulatorHst,
                                                    "m4d_int_ja_ea"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    LeakAccumulatorHst,
                                                    "m4d_int_s_leak"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    LeakAbsAccumulatorHst,
                                                    "m4d_int_s_leak_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivJMode0AccumulatorHst,
                                                    "m4d_int_divj_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivMomXMode0AccumulatorHst,
                                                    "m4d_int_divmomx_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivMomYMode0AccumulatorHst,
                                                    "m4d_int_divmomy_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivMomZMode0AccumulatorHst,
                                                    "m4d_int_divmomz_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivMomWMode0AccumulatorHst,
                                                    "m4d_int_divmomw_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivEnergyMode0AccumulatorHst,
                                                    "m4d_int_divenergy_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivPi0Mode0AccumulatorHst,
                                                    "m4d_int_divpi0_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivPixMode0AccumulatorHst,
                                                    "m4d_int_divpix_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivPiyMode0AccumulatorHst,
                                                    "m4d_int_divpiy_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivPizMode0AccumulatorHst,
                                                    "m4d_int_divpiz_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    DivPiwMode0AccumulatorHst,
                                                    "m4d_int_divpiw_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcMomXMode0AccumulatorHst,
                                                    "m4d_int_srcmomx_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcMomYMode0AccumulatorHst,
                                                    "m4d_int_srcmomy_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcMomZMode0AccumulatorHst,
                                                    "m4d_int_srcmomz_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcMomWMode0AccumulatorHst,
                                                    "m4d_int_srcmomw_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcEnergyMode0AccumulatorHst,
                                                    "m4d_int_srcenergy_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcPi0Mode0AccumulatorHst,
                                                    "m4d_int_srcpi0_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcPixMode0AccumulatorHst,
                                                    "m4d_int_srcpix_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcPiyMode0AccumulatorHst,
                                                    "m4d_int_srcpiy_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcPizMode0AccumulatorHst,
                                                    "m4d_int_srcpiz_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    SrcPiwMode0AccumulatorHst,
                                                    "m4d_int_srcpiw_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcTimelikeA0FromPiwAccumulatorHst,
      "m4d_int_src_timelike_a0_from_piw"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcTimelikeA0FromPiwAbsAccumulatorHst,
      "m4d_int_src_timelike_a0_from_piw_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcTimelikeAwFromPi0AccumulatorHst,
      "m4d_int_src_timelike_aw_from_pi0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcTimelikeAwFromPi0AbsAccumulatorHst,
      "m4d_int_src_timelike_aw_from_pi0_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMLaplacianAbsAccumulatorHst,
      "m4d_int_src_em_laplacian_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMMassAbsAccumulatorHst,
      "m4d_int_src_em_mass_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMCurrentAbsAccumulatorHst,
      "m4d_int_src_em_current_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMDampingAbsAccumulatorHst,
      "m4d_int_src_em_damping_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMSpatialMixedAbsAccumulatorHst,
      "m4d_int_src_em_spatial_mixed_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMTimelikeAbsAccumulatorHst,
      "m4d_int_src_em_timelike_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ContinuityLocalL1Hst, "m4d_cont_local_l1"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ContinuityLocalL2Hst, "m4d_cont_local_l2"));
  hst_vars.emplace_back(
      parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::max,
                                  ContinuityLocalMaxAbsHst, "m4d_cont_local_max_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ContinuityMode0L1Hst, "m4d_cont_mode0_l1"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ContinuityMode0L2Hst, "m4d_cont_mode0_l2"));
  hst_vars.emplace_back(
      parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::max,
                                  ContinuityMode0MaxAbsHst, "m4d_cont_mode0_max_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    BraneE2Hst, "m4d_brane_e2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    BraneB2Hst, "m4d_brane_b2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BraneEParallel2Hst, "m4d_brane_epar2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    MixedEw2Hst, "m4d_mixed_ew2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    MixedC2Hst, "m4d_mixed_c2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    BulkEMEnergyHst,
                                                    "m4d_em_u_bulk"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ResolvedEMEnergyHst,
                                                    "m4d_em_u_resolved"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    HelicitySubHst,
                                                    "m4d_helicity_sub"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EdotBSubHst,
                                                    "m4d_edotb_sub"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMPoyntingWHst,
                                                    "m4d_em_sw"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMLeakageWHst,
                                                    "m4d_em_leak_w"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    Ay0SpanHst, "m4d_psi0_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedSpanHst,
                                                    "m4d_psi_proj_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    PulseCentroidXHst,
                                                    "m4d_pulse_xc"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum,
      [](MeshData<Real> *md) { return EMModeComponentIntegral(md, 0, kCompA0); },
      "m4d_pi0_mode_0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum,
      [](MeshData<Real> *md) { return EMModeComponentIntegral(md, 0, kCompAX); },
      "m4d_pix_mode_0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum,
      [](MeshData<Real> *md) { return EMModeComponentIntegral(md, 0, kCompAY); },
      "m4d_piy_mode_0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum,
      [](MeshData<Real> *md) { return EMModeComponentIntegral(md, 0, kCompAZ); },
      "m4d_piz_mode_0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum,
      [](MeshData<Real> *md) { return EMModeComponentIntegral(md, 0, kCompAW); },
      "m4d_piw_mode_0"));

  const int n_modes = pkg->Param<int>("n_modes");
  for (int n = 0; n < n_modes; ++n) {
    const int mode_idx = n;
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return FieldModeL2Integral(md, "em4d_a", mode_idx);
        },
        "m4d_em_a2_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return FieldModeL2Integral(md, "em4d_pi", mode_idx);
        },
        "m4d_em_pi2_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) { return JwModeL2Integral(md, mode_idx); },
        "m4d_jw_mode_l2_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) { return JwModeIntegral(md, mode_idx); },
        "m4d_jw_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) { return ChargeModeIntegral(md, mode_idx); },
        "m4d_charge_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return SpeciesSummedModeIntegral(md, mode_idx, kPlasmaMomX);
        },
        "m4d_momx_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return SpeciesSummedModeIntegral(md, mode_idx, kPlasmaMomY);
        },
        "m4d_momy_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return SpeciesSummedModeIntegral(md, mode_idx, kPlasmaMomZ);
        },
        "m4d_momz_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return SpeciesSummedModeIntegral(md, mode_idx, kPlasmaMomW);
        },
        "m4d_momw_mode_" + std::to_string(mode_idx)));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::sum,
        [mode_idx](MeshData<Real> *md) {
          return SpeciesSummedModeIntegral(md, mode_idx, kPlasmaEnergy);
        },
        "m4d_energy_mode_" + std::to_string(mode_idx)));
  }

  pkg->AddParam<>(parthenon::hist_param_key, hst_vars, true);
}

} // namespace Modes4D
