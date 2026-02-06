#include "diagnostics4d.hpp"

#include <string>
#include <vector>

#include "outputs/outputs.hpp"

namespace Modes4D {
using namespace parthenon::package::prelude;

namespace {

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
  pkg->AddParam<double>("diag/int_jw_ew", 0.0);
  pkg->AddParam<double>("diag/int_s_leak", 0.0);

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
