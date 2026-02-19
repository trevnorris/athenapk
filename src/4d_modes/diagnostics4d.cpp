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

enum class SpillbackEMFQuantity {
  VWCMagnitude,
  VWCParallelMagnitude,
  CovVxBMagnitude,
  CovVxBParallelMagnitude,
  JwL1,
  JwL2,
  EwL1,
  EwL2
};

enum class HelicitySubscaleQuantity {
  HSub,
  EdotBSub
};

enum class TransverseEnergyQuantity {
  PoyntingW,
  LeakageW
};

enum class ProjectionParity {
  All,
  EvenOnly,
  OddOnly
};

bool ProjectionParityAcceptsMode(const ProjectionParity parity, const int mode_idx) {
  if (parity == ProjectionParity::All) {
    return true;
  }
  if (parity == ProjectionParity::EvenOnly) {
    return (mode_idx % 2) == 0;
  }
  return (mode_idx % 2) != 0;
}

std::vector<Real> BranePointCoefficients(const ModeTables &tables, const int n_modes,
                                         const Real center_w);

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

Real SrcEMGeometryShiftAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_geometry_shift_abs");
}

Real SrcEMSpatialMixedAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_spatial_mixed_abs");
}

Real SrcEMTimelikeAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_timelike_abs");
}

Real SrcEMGaugeAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_gauge_abs");
}

Real SrcEMGaugeAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_gauge");
}

Real SrcEMGaugeMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_gauge_mode0");
}

Real SrcEMGaugeMode0AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_src_em_gauge_mode0_abs");
}

Real RHSA0Mode0LapAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_lap");
}

Real RHSA0Mode0MassAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_mass");
}

Real RHSA0Mode0CurrentAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_current");
}

Real RHSA0Mode0DampingAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_damping");
}

Real RHSA0Mode0TimelikeAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_timelike");
}

Real RHSA0Mode0GaugeAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_gauge");
}

Real RHSA0Mode0TotalAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_a0_mode0_total");
}

Real RHSAyMode0LapAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_lap");
}

Real RHSAyMode0MassAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_mass");
}

Real RHSAyMode0CurrentAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_current");
}

Real RHSAyMode0DampingAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_damping");
}

Real RHSAyMode0SpatialMixedAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_spatial_mixed");
}

Real RHSAyMode0EvenBridgeAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_even_bridge");
}

Real RHSAyMode0EvenBridgeAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_even_bridge_abs");
}

Real RHSAyMode0GeometryShiftAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_geometry_shift");
}

Real RHSAyMode0GeometryShiftAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_geometry_shift_abs");
}

Real BridgePowerMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_bridge_power_mode0");
}

Real BridgePowerMode0AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_bridge_power_mode0_abs");
}

Real GeometryShiftPowerMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_geometry_shift_power_mode0");
}

Real GeometryShiftPowerMode0AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_geometry_shift_power_mode0_abs");
}

Real RHSAwMode2EvenBridgeAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode2_even_bridge");
}

Real RHSAwMode2EvenBridgeAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode2_even_bridge_abs");
}

Real BridgePowerMode2AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_bridge_power_mode2");
}

Real BridgePowerMode2AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_bridge_power_mode2_abs");
}

Real BridgePowerSumAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_bridge_power_sum");
}

Real BridgePowerSumAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_bridge_power_sum_abs");
}

Real ResponseW0Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0");
}

Real ResponseW0DotHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_dot");
}

Real ResponseW0DriveHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_drive");
}

Real ResponseW0DriveReservoirHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_drive_reservoir");
}

Real ResponseW0DriveBridgeHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_drive_bridge");
}

Real ResponseW0DriveParityHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_drive_parity");
}

Real ResponseW0ForceHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_force");
}

Real ResponseW0EnergyHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_energy");
}

Real ResponseW0GeometryCenterHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_geometry_center");
}

Real ResponseW0LambdaFractionHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_lambda_fraction");
}

Real ResponseW0LambdaEffHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_w0_lambda_eff");
}

Real ResponseLambdaFractionHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_fraction");
}

Real ResponseLambdaDotHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_dot");
}

Real ResponseLambdaDriveHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_drive");
}

Real ResponseLambdaDriveReservoirHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_drive_reservoir");
}

Real ResponseLambdaDriveBridgeHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_drive_bridge");
}

Real ResponseLambdaForceHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_force");
}

Real ResponseLambdaEnergyHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/response_lambda_energy");
}

Real ResponseW0PowerDriveAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_w0_power_drive");
}

Real ResponseW0PowerStiffnessAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_w0_power_stiffness");
}

Real ResponseW0PowerDampingAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_w0_power_damping");
}

Real ResponseW0PowerNetAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_w0_power_net");
}

Real ResponseLambdaPowerDriveAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_lambda_power_drive");
}

Real ResponseLambdaPowerStiffnessAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_lambda_power_stiffness");
}

Real ResponseLambdaPowerDampingAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_lambda_power_damping");
}

Real ResponseLambdaPowerNetAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_lambda_power_net");
}

Real RHSAyMode0W0ResponseAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_w0_response");
}

Real RHSAyMode0W0ResponseAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_w0_response_abs");
}

Real RHSAwMode1W0ResponseAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode1_w0_response");
}

Real RHSAwMode1W0ResponseAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode1_w0_response_abs");
}

Real ResponseBridgePowerMode0AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_bridge_power_mode0");
}

Real ResponseBridgePowerMode0AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_bridge_power_mode0_abs");
}

Real ResponseBridgePowerMode1AccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_bridge_power_mode1");
}

Real ResponseBridgePowerMode1AbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_bridge_power_mode1_abs");
}

Real ResponseBridgePowerSumAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_bridge_power_sum");
}

Real ResponseBridgePowerSumAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_response_bridge_power_sum_abs");
}

Real RHSAyMode0TotalAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_ay_mode0_total");
}

Real RHSAwMode0LapAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode0_lap");
}

Real RHSAwMode0CurrentAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode0_current");
}

Real RHSAwMode0DampingAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode0_damping");
}

Real RHSAwMode0SpatialMixedAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode0_spatial_mixed");
}

Real RHSAwMode0TimelikeAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode0_timelike");
}

Real RHSAwMode0TotalAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/int_rhs_aw_mode0_total");
}

Real GaugeL1Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/gauge_l1");
}

Real GaugeL2Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/gauge_l2");
}

Real GaugeMaxAbsHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/gauge_max_abs");
}

Real GaugeMode0L1Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/gauge_mode0_l1");
}

Real GaugeMode0L2Hst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/gauge_mode0_l2");
}

Real GaugeMode0MaxAbsHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return pkg->Param<double>("diag/gauge_mode0_max_abs");
}

Real DivMode0PlasmaAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_divmomx_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomw_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divenergy_mode0"));
}

Real DivMode0EMAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_divpi0_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpix_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpiy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpiz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpiw_mode0"));
}

Real DivMode0TotalAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_divj_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomx_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divmomw_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divenergy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpi0_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpix_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpiy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpiz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_divpiw_mode0"));
}

Real SrcMode0PlasmaAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_srcmomx_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcmomy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcmomz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcmomw_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcenergy_mode0"));
}

Real SrcMode0EMAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_srcpi0_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpix_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpiy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpiz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpiw_mode0"));
}

Real SrcMode0TimelikeAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_src_timelike_a0_from_piw")) +
         std::abs(pkg->Param<double>("diag/int_src_timelike_aw_from_pi0"));
}

Real SrcMode0TotalAbsAccumulatorHst(MeshData<Real> *md) {
  auto pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return std::abs(pkg->Param<double>("diag/int_srcmomx_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcmomy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcmomz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcmomw_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcenergy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpi0_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpix_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpiy_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpiz_mode0")) +
         std::abs(pkg->Param<double>("diag/int_srcpiw_mode0")) +
         std::abs(pkg->Param<double>("diag/int_src_timelike_a0_from_piw")) +
         std::abs(pkg->Param<double>("diag/int_src_timelike_aw_from_pi0"));
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

Real SpillbackEMFIntegral(MeshData<Real> *md, SpillbackEMFQuantity quantity) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const Real z_int = tables.Lambda() * std::sqrt(std::acos(-1.0));
  const Real inv_z_int = 1.0 / z_int;
  const Real rho_floor = modes_pkg->Param<double>("plasma4d/rho_floor");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");

  Real integral = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();
    auto &coords = pmb->coords;

    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();
    auto pi = bd->Get("em4d_pi").data.GetHostMirrorAndCopy();
    auto plasma = bd->Get("plasma4d_cons").data.GetHostMirrorAndCopy();

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

    std::vector<Real> rho_e_modes(n_modes, 0.0);
    std::vector<Real> momx_e_modes(n_modes, 0.0);
    std::vector<Real> momy_e_modes(n_modes, 0.0);
    std::vector<Real> momz_e_modes(n_modes, 0.0);
    std::vector<Real> momw_e_modes(n_modes, 0.0);
    std::vector<Real> momw_i_modes(n_modes, 0.0);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real dx = has_x ? coords.Dxc<1>(i) : 1.0;
          const Real dy = has_y ? coords.Dxc<2>(j) : 1.0;
          const Real dz = has_z ? coords.Dxc<3>(k) : 1.0;

          for (int n = 0; n < n_modes; ++n) {
            const int off = 5 * n;
            const int ion_base = PlasmaIndex(0, n, 0, n_modes);
            const int ele_base = PlasmaIndex(1, n, 0, n_modes);

            a0_modes[n] = a(off + kCompA0, k, j, i);
            ax_modes[n] = a(off + kCompAX, k, j, i);
            ay_modes[n] = a(off + kCompAY, k, j, i);
            az_modes[n] = a(off + kCompAZ, k, j, i);
            piw_modes[n] = pi(off + kCompAW, k, j, i);
            pix_modes[n] = pi(off + kCompAX, k, j, i);
            piy_modes[n] = pi(off + kCompAY, k, j, i);
            piz_modes[n] = pi(off + kCompAZ, k, j, i);

            rho_e_modes[n] = plasma(ele_base + kPlasmaRho, k, j, i);
            momx_e_modes[n] = plasma(ele_base + kPlasmaMomX, k, j, i);
            momy_e_modes[n] = plasma(ele_base + kPlasmaMomY, k, j, i);
            momz_e_modes[n] = plasma(ele_base + kPlasmaMomZ, k, j, i);
            momw_e_modes[n] = plasma(ele_base + kPlasmaMomW, k, j, i);
            momw_i_modes[n] = plasma(ion_base + kPlasmaMomW, k, j, i);

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

          Real vx_bar = 0.0;
          Real vy_bar = 0.0;
          Real vz_bar = 0.0;
          Real bx_bar = 0.0;
          Real by_bar = 0.0;
          Real bz_bar = 0.0;
          Real vxb_x_bar = 0.0;
          Real vxb_y_bar = 0.0;
          Real vxb_z_bar = 0.0;
          Real vwc_x_bar = 0.0;
          Real vwc_y_bar = 0.0;
          Real vwc_z_bar = 0.0;
          Real jw_l1 = 0.0;
          Real jw_l2_sq = 0.0;
          Real ew_l1 = 0.0;
          Real ew_l2_sq = 0.0;

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
            Real rho_e_node = 0.0;
            Real momx_e_node = 0.0;
            Real momy_e_node = 0.0;
            Real momz_e_node = 0.0;
            Real momw_e_node = 0.0;
            Real momw_i_node = 0.0;

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
              rho_e_node += rho_e_modes[n] * phi_nq;
              momx_e_node += momx_e_modes[n] * phi_nq;
              momy_e_node += momy_e_modes[n] * phi_nq;
              momz_e_node += momz_e_modes[n] * phi_nq;
              momw_e_node += momw_e_modes[n] * phi_nq;
              momw_i_node += momw_i_modes[n] * phi_nq;
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

            const Real rho_e_safe = std::max(rho_e_node, rho_floor);
            const Real vx = momx_e_node / rho_e_safe;
            const Real vy = momy_e_node / rho_e_safe;
            const Real vz = momz_e_node / rho_e_safe;
            const Real vw = momw_e_node / rho_e_safe;
            const Real vxb_x = (vy * bz) - (vz * by);
            const Real vxb_y = (vz * bx) - (vx * bz);
            const Real vxb_z = (vx * by) - (vy * bx);
            const Real vwc_x = -vw * cx;
            const Real vwc_y = -vw * cy;
            const Real vwc_z = -vw * cz;
            const Real jw_node = (qom_ion * momw_i_node) + (qom_electron * momw_e_node);

            vx_bar += weights[q] * vx;
            vy_bar += weights[q] * vy;
            vz_bar += weights[q] * vz;
            bx_bar += weights[q] * bx;
            by_bar += weights[q] * by;
            bz_bar += weights[q] * bz;
            vxb_x_bar += weights[q] * vxb_x;
            vxb_y_bar += weights[q] * vxb_y;
            vxb_z_bar += weights[q] * vxb_z;
            vwc_x_bar += weights[q] * vwc_x;
            vwc_y_bar += weights[q] * vwc_y;
            vwc_z_bar += weights[q] * vwc_z;
            jw_l1 += weights[q] * std::abs(jw_node);
            jw_l2_sq += weights[q] * jw_node * jw_node;
            ew_l1 += weights[q] * std::abs(ew);
            ew_l2_sq += weights[q] * ew * ew;
          }

          vx_bar *= inv_z_int;
          vy_bar *= inv_z_int;
          vz_bar *= inv_z_int;
          bx_bar *= inv_z_int;
          by_bar *= inv_z_int;
          bz_bar *= inv_z_int;
          vxb_x_bar *= inv_z_int;
          vxb_y_bar *= inv_z_int;
          vxb_z_bar *= inv_z_int;
          vwc_x_bar *= inv_z_int;
          vwc_y_bar *= inv_z_int;
          vwc_z_bar *= inv_z_int;
          jw_l1 *= inv_z_int;
          jw_l2_sq *= inv_z_int;
          ew_l1 *= inv_z_int;
          ew_l2_sq *= inv_z_int;

          const Real bbar2 = (bx_bar * bx_bar) + (by_bar * by_bar) + (bz_bar * bz_bar);
          const Real bbar_mag = std::sqrt(bbar2) + 1.0e-30;
          const Real vbar_x_b_x = (vy_bar * bz_bar) - (vz_bar * by_bar);
          const Real vbar_x_b_y = (vz_bar * bx_bar) - (vx_bar * bz_bar);
          const Real vbar_x_b_z = (vx_bar * by_bar) - (vy_bar * bx_bar);
          const Real cov_x = vxb_x_bar - vbar_x_b_x;
          const Real cov_y = vxb_y_bar - vbar_x_b_y;
          const Real cov_z = vxb_z_bar - vbar_x_b_z;

          const Real vwc_mag =
              std::sqrt((vwc_x_bar * vwc_x_bar) + (vwc_y_bar * vwc_y_bar) +
                        (vwc_z_bar * vwc_z_bar));
          const Real vwc_par_mag =
              std::abs((vwc_x_bar * bx_bar) + (vwc_y_bar * by_bar) + (vwc_z_bar * bz_bar)) /
              bbar_mag;
          const Real cov_mag =
              std::sqrt((cov_x * cov_x) + (cov_y * cov_y) + (cov_z * cov_z));
          const Real cov_par_mag =
              std::abs((cov_x * bx_bar) + (cov_y * by_bar) + (cov_z * bz_bar)) / bbar_mag;
          const Real jw_l2 = std::sqrt(std::max(jw_l2_sq, 0.0));
          const Real ew_l2 = std::sqrt(std::max(ew_l2_sq, 0.0));

          Real val = 0.0;
          if (quantity == SpillbackEMFQuantity::VWCMagnitude) {
            val = vwc_mag;
          } else if (quantity == SpillbackEMFQuantity::VWCParallelMagnitude) {
            val = vwc_par_mag;
          } else if (quantity == SpillbackEMFQuantity::CovVxBMagnitude) {
            val = cov_mag;
          } else if (quantity == SpillbackEMFQuantity::CovVxBParallelMagnitude) {
            val = cov_par_mag;
          } else if (quantity == SpillbackEMFQuantity::JwL1) {
            val = jw_l1;
          } else if (quantity == SpillbackEMFQuantity::JwL2) {
            val = jw_l2;
          } else if (quantity == SpillbackEMFQuantity::EwL1) {
            val = ew_l1;
          } else {
            val = ew_l2;
          }

          integral += val * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return integral;
}

Real EMFVwCAbsHst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::VWCMagnitude);
}

Real EMFVwCParallelAbsHst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::VWCParallelMagnitude);
}

Real EMFCovVxBAbsHst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::CovVxBMagnitude);
}

Real EMFCovVxBParallelAbsHst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::CovVxBParallelMagnitude);
}

Real JwL1Hst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::JwL1);
}

Real JwL2Hst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::JwL2);
}

Real EwL1Hst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::EwL1);
}

Real EwL2Hst(MeshData<Real> *md) {
  return SpillbackEMFIntegral(md, SpillbackEMFQuantity::EwL2);
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

Real ProjectionCenterWHst(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  return modes_pkg->Param<double>("diag/projection_center_w");
}

std::vector<Real> ProjectionCoefficients(const ModeTables &tables, const int n_modes,
                                         const std::string &kernel,
                                         const Real sigma_factor,
                                         const Real center_w) {
  if (kernel == "point") {
    return BranePointCoefficients(tables, n_modes, center_w);
  }

  const int n_quad = tables.NumQuadrature();
  const auto &weights = tables.Weights();
  const Real lambda = tables.Lambda();
  const Real z_int = lambda * std::sqrt(std::acos(-1.0));
  std::vector<Real> coeff(static_cast<size_t>(n_modes), 0.0);

  if (kernel == "gaussian") {
    const auto &nodes = tables.Nodes();
    const Real lambda_sq = lambda * lambda;
    const Real sigma_sq = sigma_factor * sigma_factor;
    const Real w_norm = sigma_factor * z_int;
    for (int n = 0; n < n_modes; ++n) {
      Real mode_coeff = 0.0;
      for (int q = 0; q < n_quad; ++q) {
        const Real w = nodes[q];
        const Real z = std::exp(-(w * w) / lambda_sq);
        if (z <= 0.0) {
          continue;
        }
        const Real dw = w - center_w;
        const Real w_kernel =
            std::exp(-(dw * dw) / (sigma_sq * lambda_sq)) / w_norm;
        mode_coeff += weights[q] * (w_kernel / z) * tables.Phi(n, q);
      }
      coeff[static_cast<size_t>(n)] = mode_coeff;
    }
    return coeff;
  }

  // Default: matched kernel W = Z(w-center_w) / Z_int.
  const auto &nodes = tables.Nodes();
  const Real lambda_sq = lambda * lambda;
  for (int n = 0; n < n_modes; ++n) {
    Real mode_coeff = 0.0;
    for (int q = 0; q < n_quad; ++q) {
      if (std::abs(center_w) <= 1.0e-14) {
        mode_coeff += weights[q] * tables.Phi(n, q);
      } else {
        const Real w = nodes[q];
        const Real z = std::exp(-(w * w) / lambda_sq);
        if (z <= 0.0) {
          continue;
        }
        const Real dw = w - center_w;
        const Real z_shift = std::exp(-(dw * dw) / lambda_sq);
        mode_coeff += weights[q] * ((z_shift / z) / z_int) * tables.Phi(n, q);
      }
    }
    if (std::abs(center_w) <= 1.0e-14) {
      coeff[static_cast<size_t>(n)] = mode_coeff / z_int;
    } else {
      coeff[static_cast<size_t>(n)] = mode_coeff;
    }
  }
  return coeff;
}

Real ProjectionSpanWithKernelHst(MeshData<Real> *md, const std::string &kernel,
                                 const Real sigma_factor,
                                 const ProjectionParity parity) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const Real center_w = modes_pkg->Param<double>("diag/projection_center_w");
  const auto proj_coeff =
      ProjectionCoefficients(tables, n_modes, kernel, sigma_factor, center_w);

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
            if (!ProjectionParityAcceptsMode(parity, n)) {
              continue;
            }
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

Real AyProjectedSpanHst(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const auto &projection_kernel = modes_pkg->Param<std::string>("diag/projection_kernel");
  const Real projection_sigma_factor =
      modes_pkg->Param<double>("diag/projection_sigma_factor");
  return ProjectionSpanWithKernelHst(md, projection_kernel, projection_sigma_factor,
                                     ProjectionParity::All);
}

Real AyProjectedEvenSpanHst(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const auto &projection_kernel = modes_pkg->Param<std::string>("diag/projection_kernel");
  const Real projection_sigma_factor =
      modes_pkg->Param<double>("diag/projection_sigma_factor");
  return ProjectionSpanWithKernelHst(md, projection_kernel, projection_sigma_factor,
                                     ProjectionParity::EvenOnly);
}

Real AyProjectedOddSpanHst(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const auto &projection_kernel = modes_pkg->Param<std::string>("diag/projection_kernel");
  const Real projection_sigma_factor =
      modes_pkg->Param<double>("diag/projection_sigma_factor");
  return ProjectionSpanWithKernelHst(md, projection_kernel, projection_sigma_factor,
                                     ProjectionParity::OddOnly);
}

Real AyProjectedMatchedSpanHst(MeshData<Real> *md) {
  return ProjectionSpanWithKernelHst(md, "matched", 1.0, ProjectionParity::All);
}

Real AyProjectedPointSpanHst(MeshData<Real> *md) {
  return ProjectionSpanWithKernelHst(md, "point", 1.0, ProjectionParity::All);
}

Real AyProjectedGaussianSpanHst(MeshData<Real> *md) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const Real projection_sigma_factor =
      modes_pkg->Param<double>("diag/projection_sigma_factor");
  return ProjectionSpanWithKernelHst(md, "gaussian", projection_sigma_factor,
                                     ProjectionParity::All);
}

std::vector<Real> BranePointCoefficients(const ModeTables &tables, const int n_modes,
                                         const Real center_w) {
  std::vector<Real> coeff(static_cast<size_t>(n_modes), 0.0);
  if (n_modes <= 0) {
    return coeff;
  }

  // Evaluate phi_n(center_w) via Hermite recurrence in physicists' convention.
  std::vector<Real> hermite_zero(static_cast<size_t>(n_modes), 0.0);
  const Real x = center_w / tables.Lambda();
  hermite_zero[0] = 1.0;
  if (n_modes > 1) {
    hermite_zero[1] = 2.0 * x;
  }
  for (int n = 1; n + 1 < n_modes; ++n) {
    hermite_zero[static_cast<size_t>(n + 1)] =
        (2.0 * x * hermite_zero[static_cast<size_t>(n)]) -
        (2.0 * static_cast<Real>(n) * hermite_zero[static_cast<size_t>(n - 1)]);
  }

  const Real normalizer_factor = tables.Lambda() * std::sqrt(std::acos(-1.0));
  for (int n = 0; n < n_modes; ++n) {
    const Real two_to_n = std::ldexp(1.0, n);
    const Real n_factorial = std::tgamma(static_cast<Real>(n) + 1.0);
    const Real normalization = std::sqrt(normalizer_factor * two_to_n * n_factorial);
    coeff[static_cast<size_t>(n)] = hermite_zero[static_cast<size_t>(n)] / normalization;
  }
  return coeff;
}

Real AyBranePointSpanWithParityHst(MeshData<Real> *md, const ProjectionParity parity) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const Real center_w = modes_pkg->Param<double>("diag/projection_center_w");
  const auto brane_coeff = BranePointCoefficients(tables, n_modes, center_w);

  Real ay_w0_min = std::numeric_limits<Real>::max();
  Real ay_w0_max = std::numeric_limits<Real>::lowest();

  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto a = bd->Get("em4d_a").data.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          Real ay_w0 = 0.0;
          for (int n = 0; n < n_modes; ++n) {
            if (!ProjectionParityAcceptsMode(parity, n)) {
              continue;
            }
            const int off = 5 * n;
            ay_w0 += brane_coeff[static_cast<size_t>(n)] * a(off + kCompAY, k, j, i);
          }
          ay_w0_min = std::min(ay_w0_min, ay_w0);
          ay_w0_max = std::max(ay_w0_max, ay_w0);
        }
      }
    }
  }

  return (ay_w0_max > ay_w0_min) ? (ay_w0_max - ay_w0_min) : 0.0;
}

Real AyBranePointSpanHst(MeshData<Real> *md) {
  return AyBranePointSpanWithParityHst(md, ProjectionParity::All);
}

Real AyBranePointEvenSpanHst(MeshData<Real> *md) {
  return AyBranePointSpanWithParityHst(md, ProjectionParity::EvenOnly);
}

Real AyBranePointOddSpanHst(MeshData<Real> *md) {
  return AyBranePointSpanWithParityHst(md, ProjectionParity::OddOnly);
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

Real FieldParityL2Integral(MeshData<Real> *md, const std::string &field_name,
                           const bool even_modes) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");

  Real sum = 0.0;
  for (int b = 0; b < md->NumBlocks(); ++b) {
    auto &bd = md->GetBlockData(b);
    auto field = bd->Get(field_name).data.GetHostMirrorAndCopy();
    auto &coords = bd->GetBlockPointer()->coords;

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          Real val2 = 0.0;
          for (int n = even_modes ? 0 : 1; n < n_modes; n += 2) {
            const int off = 5 * n;
            for (int c = 0; c < 5; ++c) {
              const Real val = field(off + c, k, j, i);
              val2 += val * val;
            }
          }
          sum += 0.5 * val2 * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real JwParityL2Integral(MeshData<Real> *md, const bool even_modes) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");

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
          Real jw2 = 0.0;
          for (int n = even_modes ? 0 : 1; n < n_modes; n += 2) {
            const int ion_base = PlasmaIndex(0, n, 0, n_modes);
            const int ele_base = PlasmaIndex(1, n, 0, n_modes);
            const Real ion_momw = plasma(ion_base + kPlasmaMomW, k, j, i);
            const Real ele_momw = plasma(ele_base + kPlasmaMomW, k, j, i);
            const Real jw = (qom_ion * ion_momw) + (qom_electron * ele_momw);
            jw2 += jw * jw;
          }
          sum += 0.5 * jw2 * coords.CellVolume(k, j, i);
        }
      }
    }
  }
  return sum;
}

Real EMA2EvenHst(MeshData<Real> *md) {
  return FieldParityL2Integral(md, "em4d_a", true);
}

Real EMA2OddHst(MeshData<Real> *md) {
  return FieldParityL2Integral(md, "em4d_a", false);
}

Real EMPi2EvenHst(MeshData<Real> *md) {
  return FieldParityL2Integral(md, "em4d_pi", true);
}

Real EMPi2OddHst(MeshData<Real> *md) {
  return FieldParityL2Integral(md, "em4d_pi", false);
}

Real JwL2EvenHst(MeshData<Real> *md) {
  return JwParityL2Integral(md, true);
}

Real JwL2OddHst(MeshData<Real> *md) {
  return JwParityL2Integral(md, false);
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
  pkg->AddParam<double>("diag/int_src_em_geometry_shift_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_spatial_mixed_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_timelike_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_gauge", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_gauge_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_gauge_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_src_em_gauge_mode0_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_lap", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_mass", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_current", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_damping", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_timelike", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_gauge", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_a0_mode0_total", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_lap", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_mass", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_current", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_damping", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_spatial_mixed", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_even_bridge", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_even_bridge_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_geometry_shift", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_geometry_shift_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_bridge_power_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_bridge_power_mode0_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_geometry_shift_power_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_geometry_shift_power_mode0_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode2_even_bridge", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode2_even_bridge_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_bridge_power_mode2", 0.0, true);
  pkg->AddParam<double>("diag/int_bridge_power_mode2_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_bridge_power_sum", 0.0, true);
  pkg->AddParam<double>("diag/int_bridge_power_sum_abs", 0.0, true);
  pkg->AddParam<double>("diag/response_w0", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_dot", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_drive", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_drive_reservoir", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_drive_bridge", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_drive_parity", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_force", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_energy", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_geometry_center", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_lambda_fraction", 0.0, true);
  pkg->AddParam<double>("diag/response_w0_lambda_eff", 0.0, true);
  pkg->AddParam<bool>("diag/response_lambda_initialized", false, true);
  pkg->AddParam<double>("diag/response_lambda_fraction", 0.0, true);
  pkg->AddParam<double>("diag/response_lambda_dot", 0.0, true);
  pkg->AddParam<double>("diag/response_lambda_drive", 0.0, true);
  pkg->AddParam<double>("diag/response_lambda_drive_reservoir", 0.0, true);
  pkg->AddParam<double>("diag/response_lambda_drive_bridge", 0.0, true);
  pkg->AddParam<double>("diag/response_lambda_force", 0.0, true);
  pkg->AddParam<double>("diag/response_lambda_energy", 0.0, true);
  pkg->AddParam<double>("diag/projection_center_w", 0.0, true);
  pkg->AddParam<bool>("diag/response_w0_initialized", false, true);
  pkg->AddParam<double>("diag/int_response_w0_power_drive", 0.0, true);
  pkg->AddParam<double>("diag/int_response_w0_power_stiffness", 0.0, true);
  pkg->AddParam<double>("diag/int_response_w0_power_damping", 0.0, true);
  pkg->AddParam<double>("diag/int_response_w0_power_net", 0.0, true);
  pkg->AddParam<double>("diag/int_response_lambda_power_drive", 0.0, true);
  pkg->AddParam<double>("diag/int_response_lambda_power_stiffness", 0.0, true);
  pkg->AddParam<double>("diag/int_response_lambda_power_damping", 0.0, true);
  pkg->AddParam<double>("diag/int_response_lambda_power_net", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_w0_response", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_w0_response_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode1_w0_response", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode1_w0_response_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_response_bridge_power_mode0", 0.0, true);
  pkg->AddParam<double>("diag/int_response_bridge_power_mode0_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_response_bridge_power_mode1", 0.0, true);
  pkg->AddParam<double>("diag/int_response_bridge_power_mode1_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_response_bridge_power_sum", 0.0, true);
  pkg->AddParam<double>("diag/int_response_bridge_power_sum_abs", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_ay_mode0_total", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode0_lap", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode0_current", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode0_damping", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode0_spatial_mixed", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode0_timelike", 0.0, true);
  pkg->AddParam<double>("diag/int_rhs_aw_mode0_total", 0.0, true);
  pkg->AddParam<double>("diag/continuity_local_l1", 0.0, true);
  pkg->AddParam<double>("diag/continuity_local_l2", 0.0, true);
  pkg->AddParam<double>("diag/continuity_local_max_abs", 0.0, true);
  pkg->AddParam<double>("diag/continuity_mode0_l1", 0.0, true);
  pkg->AddParam<double>("diag/continuity_mode0_l2", 0.0, true);
  pkg->AddParam<double>("diag/continuity_mode0_max_abs", 0.0, true);
  pkg->AddParam<double>("diag/gauge_l1", 0.0, true);
  pkg->AddParam<double>("diag/gauge_l2", 0.0, true);
  pkg->AddParam<double>("diag/gauge_max_abs", 0.0, true);
  pkg->AddParam<double>("diag/gauge_mode0_l1", 0.0, true);
  pkg->AddParam<double>("diag/gauge_mode0_l2", 0.0, true);
  pkg->AddParam<double>("diag/gauge_mode0_max_abs", 0.0, true);

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
      parthenon::UserHistoryOperation::sum, SrcEMGeometryShiftAbsAccumulatorHst,
      "m4d_int_src_em_geometry_shift_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMSpatialMixedAbsAccumulatorHst,
      "m4d_int_src_em_spatial_mixed_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMTimelikeAbsAccumulatorHst,
      "m4d_int_src_em_timelike_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMGaugeAccumulatorHst,
      "m4d_int_src_em_gauge"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMGaugeAbsAccumulatorHst,
      "m4d_int_src_em_gauge_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMGaugeMode0AccumulatorHst,
      "m4d_int_src_em_gauge_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcEMGaugeMode0AbsAccumulatorHst,
      "m4d_int_src_em_gauge_mode0_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0LapAccumulatorHst,
      "m4d_int_rhs_a0_mode0_lap"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0MassAccumulatorHst,
      "m4d_int_rhs_a0_mode0_mass"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0CurrentAccumulatorHst,
      "m4d_int_rhs_a0_mode0_current"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0DampingAccumulatorHst,
      "m4d_int_rhs_a0_mode0_damping"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0TimelikeAccumulatorHst,
      "m4d_int_rhs_a0_mode0_timelike"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0GaugeAccumulatorHst,
      "m4d_int_rhs_a0_mode0_gauge"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSA0Mode0TotalAccumulatorHst,
      "m4d_int_rhs_a0_mode0_total"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0LapAccumulatorHst,
      "m4d_int_rhs_ay_mode0_lap"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0MassAccumulatorHst,
      "m4d_int_rhs_ay_mode0_mass"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0CurrentAccumulatorHst,
      "m4d_int_rhs_ay_mode0_current"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0DampingAccumulatorHst,
      "m4d_int_rhs_ay_mode0_damping"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0SpatialMixedAccumulatorHst,
      "m4d_int_rhs_ay_mode0_spatial_mixed"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0EvenBridgeAccumulatorHst,
      "m4d_int_rhs_ay_mode0_even_bridge"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0EvenBridgeAbsAccumulatorHst,
      "m4d_int_rhs_ay_mode0_even_bridge_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0GeometryShiftAccumulatorHst,
      "m4d_int_rhs_ay_mode0_geometry_shift"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0GeometryShiftAbsAccumulatorHst,
      "m4d_int_rhs_ay_mode0_geometry_shift_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BridgePowerMode0AccumulatorHst,
      "m4d_int_bridge_power_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BridgePowerMode0AbsAccumulatorHst,
      "m4d_int_bridge_power_mode0_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, GeometryShiftPowerMode0AccumulatorHst,
      "m4d_int_geometry_shift_power_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, GeometryShiftPowerMode0AbsAccumulatorHst,
      "m4d_int_geometry_shift_power_mode0_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode2EvenBridgeAccumulatorHst,
      "m4d_int_rhs_aw_mode2_even_bridge"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode2EvenBridgeAbsAccumulatorHst,
      "m4d_int_rhs_aw_mode2_even_bridge_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BridgePowerMode2AccumulatorHst,
      "m4d_int_bridge_power_mode2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BridgePowerMode2AbsAccumulatorHst,
      "m4d_int_bridge_power_mode2_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BridgePowerSumAccumulatorHst,
      "m4d_int_bridge_power_sum"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, BridgePowerSumAbsAccumulatorHst,
      "m4d_int_bridge_power_sum_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ResponseW0Hst, "m4d_response_w0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ResponseW0DotHst,
                                                    "m4d_response_w0_dot"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ResponseW0DriveHst,
                                                    "m4d_response_w0_drive"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0DriveReservoirHst,
      "m4d_response_w0_drive_reservoir"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0DriveBridgeHst,
      "m4d_response_w0_drive_bridge"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0DriveParityHst,
      "m4d_response_w0_drive_parity"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ResponseW0ForceHst,
                                                    "m4d_response_w0_force"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ResponseW0EnergyHst,
                                                    "m4d_response_w0_energy"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0GeometryCenterHst,
      "m4d_response_w0_geometry_center"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0LambdaFractionHst,
      "m4d_response_w0_lambda_fraction"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0LambdaEffHst,
      "m4d_response_w0_lambda_eff"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaFractionHst,
      "m4d_response_lambda_fraction"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaDotHst,
      "m4d_response_lambda_dot"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaDriveHst,
      "m4d_response_lambda_drive"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaDriveReservoirHst,
      "m4d_response_lambda_drive_reservoir"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaDriveBridgeHst,
      "m4d_response_lambda_drive_bridge"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaForceHst,
      "m4d_response_lambda_force"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaEnergyHst,
      "m4d_response_lambda_energy"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    ProjectionCenterWHst,
                                                    "m4d_projection_center_w"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0PowerDriveAccumulatorHst,
      "m4d_int_response_w0_power_drive"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0PowerStiffnessAccumulatorHst,
      "m4d_int_response_w0_power_stiffness"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0PowerDampingAccumulatorHst,
      "m4d_int_response_w0_power_damping"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseW0PowerNetAccumulatorHst,
      "m4d_int_response_w0_power_net"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaPowerDriveAccumulatorHst,
      "m4d_int_response_lambda_power_drive"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaPowerStiffnessAccumulatorHst,
      "m4d_int_response_lambda_power_stiffness"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaPowerDampingAccumulatorHst,
      "m4d_int_response_lambda_power_damping"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseLambdaPowerNetAccumulatorHst,
      "m4d_int_response_lambda_power_net"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0W0ResponseAccumulatorHst,
      "m4d_int_rhs_ay_mode0_w0_response"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0W0ResponseAbsAccumulatorHst,
      "m4d_int_rhs_ay_mode0_w0_response_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode1W0ResponseAccumulatorHst,
      "m4d_int_rhs_aw_mode1_w0_response"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode1W0ResponseAbsAccumulatorHst,
      "m4d_int_rhs_aw_mode1_w0_response_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseBridgePowerMode0AccumulatorHst,
      "m4d_int_response_bridge_power_mode0"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseBridgePowerMode0AbsAccumulatorHst,
      "m4d_int_response_bridge_power_mode0_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseBridgePowerMode1AccumulatorHst,
      "m4d_int_response_bridge_power_mode1"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseBridgePowerMode1AbsAccumulatorHst,
      "m4d_int_response_bridge_power_mode1_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseBridgePowerSumAccumulatorHst,
      "m4d_int_response_bridge_power_sum"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, ResponseBridgePowerSumAbsAccumulatorHst,
      "m4d_int_response_bridge_power_sum_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAyMode0TotalAccumulatorHst,
      "m4d_int_rhs_ay_mode0_total"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode0LapAccumulatorHst,
      "m4d_int_rhs_aw_mode0_lap"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode0CurrentAccumulatorHst,
      "m4d_int_rhs_aw_mode0_current"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode0DampingAccumulatorHst,
      "m4d_int_rhs_aw_mode0_damping"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode0SpatialMixedAccumulatorHst,
      "m4d_int_rhs_aw_mode0_spatial_mixed"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode0TimelikeAccumulatorHst,
      "m4d_int_rhs_aw_mode0_timelike"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, RHSAwMode0TotalAccumulatorHst,
      "m4d_int_rhs_aw_mode0_total"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, DivMode0PlasmaAbsAccumulatorHst,
      "m4d_int_div_mode0_plasma_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, DivMode0EMAbsAccumulatorHst,
      "m4d_int_div_mode0_em_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, DivMode0TotalAbsAccumulatorHst,
      "m4d_int_div_mode0_total_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcMode0PlasmaAbsAccumulatorHst,
      "m4d_int_src_mode0_plasma_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcMode0EMAbsAccumulatorHst,
      "m4d_int_src_mode0_em_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcMode0TimelikeAbsAccumulatorHst,
      "m4d_int_src_mode0_timelike_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, SrcMode0TotalAbsAccumulatorHst,
      "m4d_int_src_mode0_total_abs"));
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
  hst_vars.emplace_back(
      parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum, GaugeL1Hst,
                                  "m4d_gauge_l1"));
  hst_vars.emplace_back(
      parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum, GaugeL2Hst,
                                  "m4d_gauge_l2"));
  hst_vars.emplace_back(
      parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::max, GaugeMaxAbsHst,
                                  "m4d_gauge_max_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, GaugeMode0L1Hst, "m4d_gauge_mode0_l1"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::sum, GaugeMode0L2Hst, "m4d_gauge_mode0_l2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(
      parthenon::UserHistoryOperation::max, GaugeMode0MaxAbsHst,
      "m4d_gauge_mode0_max_abs"));
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
                                                    EMFVwCAbsHst,
                                                    "m4d_emf_vw_c_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMFVwCParallelAbsHst,
                                                    "m4d_emf_vw_c_par_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMFCovVxBAbsHst,
                                                    "m4d_emf_cov_vxb_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMFCovVxBParallelAbsHst,
                                                    "m4d_emf_cov_vxb_par_abs"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    JwL1Hst, "m4d_jw_l1"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    JwL2Hst, "m4d_jw_l2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EwL1Hst, "m4d_ew_l1"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EwL2Hst, "m4d_ew_l2"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMA2EvenHst, "m4d_em_a2_even"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMA2OddHst, "m4d_em_a2_odd"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMPi2EvenHst, "m4d_em_pi2_even"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    EMPi2OddHst, "m4d_em_pi2_odd"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    JwL2EvenHst, "m4d_jw_l2_even"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    JwL2OddHst, "m4d_jw_l2_odd"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    Ay0SpanHst, "m4d_psi0_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedSpanHst,
                                                    "m4d_psi_proj_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedEvenSpanHst,
                                                    "m4d_psi_proj_even_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedOddSpanHst,
                                                    "m4d_psi_proj_odd_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedMatchedSpanHst,
                                                    "m4d_psi_proj_matched_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedPointSpanHst,
                                                    "m4d_psi_proj_point_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyProjectedGaussianSpanHst,
                                                    "m4d_psi_proj_gaussian_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyBranePointSpanHst,
                                                    "m4d_psi_w0_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyBranePointEvenSpanHst,
                                                    "m4d_psi_w0_even_span"));
  hst_vars.emplace_back(parthenon::HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                                    AyBranePointOddSpanHst,
                                                    "m4d_psi_w0_odd_span"));
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
