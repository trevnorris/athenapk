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

enum class BraneMixedQuantity {
  BraneE2,
  BraneB2,
  BraneEParallel2,
  MixedEw2,
  MixedC2
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

Real LeakAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_s_leak");
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
  pkg->AddParam<double>("diag/int_s_leak", 0.0, true);

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
                                                    LeakAccumulatorHst,
                                                    "m4d_int_s_leak"));
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
                                                    Ay0SpanHst, "m4d_psi0_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    PulseCentroidXHst,
                                                    "m4d_pulse_xc"));

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
  }

  pkg->AddParam<>(parthenon::hist_param_key, hst_vars, true);
}

} // namespace Modes4D
