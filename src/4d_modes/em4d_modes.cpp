#include "em4d_modes.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "mode_tables.hpp"
#include "interface/update.hpp"

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

void RegisterEMVariables(parthenon::StateDescriptor *pkg, const int n_modes,
                         const bool enable_conservative_transport) {
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

  std::vector<parthenon::MetadataFlag> pi_metadata{
      Metadata::Cell, Metadata::Independent, Metadata::FillGhost};
  if (enable_conservative_transport) {
    pi_metadata.push_back(Metadata::WithFluxes);
  }
  m = Metadata(pi_metadata, std::vector<int>({5 * n_modes}), pi_labels);
  pkg->AddField("em4d_pi", m);
}

TaskStatus AddEMTransportFluxes(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto modes_pkg = pmb->packages.Get("modes4d");
  if (!modes_pkg->Param<bool>("enabled")) {
    return TaskStatus::complete;
  }
  if (!modes_pkg->Param<bool>("em4d/use_conservative_transport")) {
    return TaskStatus::complete;
  }

  const Real c_wave = modes_pkg->Param<double>("em4d/c_wave");
  const Real c2 = c_wave * c_wave;
  const Real lambda_base = modes_pkg->Param<double>("lambda");
  Real lambda_eff = lambda_base;
  const Real lambda_diag = modes_pkg->Param<double>("diag/response_w0_lambda_eff");
  if (std::isfinite(lambda_diag) && (lambda_diag > 0.0)) {
    lambda_eff = lambda_diag;
  }
  const Real inv_lambda = 1.0 / lambda_eff;
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const auto &pi_pack =
      md->PackVariablesAndFluxes(std::vector<std::string>{"em4d_pi"},
                                 std::vector<std::string>{"em4d_pi"});
  const auto &a_pack = md->PackVariables(std::vector<std::string>{"em4d_a"});
  if (pi_pack.GetDim(4) == 0 || a_pack.GetDim(4) == 0) {
    return TaskStatus::complete;
  }

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);
  const int ndim = pmb->pmy_mesh->ndim;

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Modes4DEMFluxX1", parthenon::DevExecSpace(), 0,
      pi_pack.GetDim(5) - 1, 0, pi_pack.GetDim(4) - 1, kb.s, kb.e, jb.s, jb.e, ib.s,
      ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
        auto &pi = pi_pack(b);
        const auto &a = a_pack(b);
        const auto &coords = pi_pack.GetCoords(b);
        const Real dx = 0.5 * (coords.Dxc<1>(i - 1) + coords.Dxc<1>(i));
        const int mode = v / kNumEMComponents;
        const int comp = v - (mode * kNumEMComponents);
        Real flux = -c2 * (a(v, k, j, i) - a(v, k, j, i - 1)) / dx;

        if (comp == kCompAX && mode > 0) {
          const Real coupling = std::sqrt(2.0 * static_cast<Real>(mode)) * inv_lambda;
          const int idx_aw_prev = EMIndex(mode - 1, kCompAW);
          const Real aw_face = 0.5 * (a(idx_aw_prev, k, j, i) + a(idx_aw_prev, k, j, i - 1));
          flux -= c2 * coupling * aw_face;
        } else if (comp == kCompAW && (mode + 1) < n_modes) {
          const Real coupling =
              std::sqrt(2.0 * static_cast<Real>(mode + 1)) * inv_lambda;
          const int idx_ax_next = EMIndex(mode + 1, kCompAX);
          const Real ax_face = 0.5 * (a(idx_ax_next, k, j, i) + a(idx_ax_next, k, j, i - 1));
          flux += c2 * coupling * ax_face;
        }

        pi.flux(X1DIR, v, k, j, i) = flux;
      });

  if (ndim >= 2) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "Modes4DEMFluxX2", parthenon::DevExecSpace(), 0,
        pi_pack.GetDim(5) - 1, 0, pi_pack.GetDim(4) - 1, kb.s, kb.e, jb.s, jb.e + 1,
        ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
          auto &pi = pi_pack(b);
          const auto &a = a_pack(b);
          const auto &coords = pi_pack.GetCoords(b);
          const Real dy = 0.5 * (coords.Dxc<2>(j - 1) + coords.Dxc<2>(j));
          const int mode = v / kNumEMComponents;
          const int comp = v - (mode * kNumEMComponents);
          Real flux = -c2 * (a(v, k, j, i) - a(v, k, j - 1, i)) / dy;

          if (comp == kCompAY && mode > 0) {
            const Real coupling = std::sqrt(2.0 * static_cast<Real>(mode)) * inv_lambda;
            const int idx_aw_prev = EMIndex(mode - 1, kCompAW);
            const Real aw_face = 0.5 * (a(idx_aw_prev, k, j, i) + a(idx_aw_prev, k, j - 1, i));
            flux -= c2 * coupling * aw_face;
          } else if (comp == kCompAW && (mode + 1) < n_modes) {
            const Real coupling =
                std::sqrt(2.0 * static_cast<Real>(mode + 1)) * inv_lambda;
            const int idx_ay_next = EMIndex(mode + 1, kCompAY);
            const Real ay_face = 0.5 * (a(idx_ay_next, k, j, i) + a(idx_ay_next, k, j - 1, i));
            flux += c2 * coupling * ay_face;
          }

          pi.flux(X2DIR, v, k, j, i) = flux;
        });
  }

  if (ndim >= 3) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "Modes4DEMFluxX3", parthenon::DevExecSpace(), 0,
        pi_pack.GetDim(5) - 1, 0, pi_pack.GetDim(4) - 1, kb.s, kb.e + 1, jb.s, jb.e,
        ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
          auto &pi = pi_pack(b);
          const auto &a = a_pack(b);
          const auto &coords = pi_pack.GetCoords(b);
          const Real dz = 0.5 * (coords.Dxc<3>(k - 1) + coords.Dxc<3>(k));
          const int mode = v / kNumEMComponents;
          const int comp = v - (mode * kNumEMComponents);
          Real flux = -c2 * (a(v, k, j, i) - a(v, k - 1, j, i)) / dz;

          if (comp == kCompAZ && mode > 0) {
            const Real coupling = std::sqrt(2.0 * static_cast<Real>(mode)) * inv_lambda;
            const int idx_aw_prev = EMIndex(mode - 1, kCompAW);
            const Real aw_face = 0.5 * (a(idx_aw_prev, k, j, i) + a(idx_aw_prev, k - 1, j, i));
            flux -= c2 * coupling * aw_face;
          } else if (comp == kCompAW && (mode + 1) < n_modes) {
            const Real coupling =
                std::sqrt(2.0 * static_cast<Real>(mode + 1)) * inv_lambda;
            const int idx_az_next = EMIndex(mode + 1, kCompAZ);
            const Real az_face = 0.5 * (a(idx_az_next, k, j, i) + a(idx_az_next, k - 1, j, i));
            flux += c2 * coupling * az_face;
          }

          pi.flux(X3DIR, v, k, j, i) = flux;
        });
  }

  return TaskStatus::complete;
}

TaskStatus AccumulateEMTransportMode0Diagnostics(MeshData<Real> *md, const Real dt) {
  if (dt <= 0.0) {
    return TaskStatus::complete;
  }

  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto modes_pkg = pmb->packages.Get("modes4d");
  if (!modes_pkg->Param<bool>("enabled")) {
    return TaskStatus::complete;
  }
  if (!modes_pkg->Param<bool>("em4d/use_conservative_transport")) {
    return TaskStatus::complete;
  }

  const int n_modes = modes_pkg->Param<int>("n_modes");
  if (n_modes < 1) {
    return TaskStatus::complete;
  }

  const auto &pi_pack =
      md->PackVariablesAndFluxes(std::vector<std::string>{"em4d_pi"},
                                 std::vector<std::string>{"em4d_pi"});
  if (pi_pack.GetDim(4) == 0) {
    return TaskStatus::complete;
  }

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);
  const int ndim = pi_pack.GetNdim();

  auto accumulate_transport = [&](const int idx, const char *kernel_name) {
    Real step = 0.0;
    parthenon::par_reduce(
        DEFAULT_LOOP_PATTERN, kernel_name, parthenon::DevExecSpace(), 0,
        pi_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i,
                      Real &local_sum) {
          const auto &pi = pi_pack(b);
          const auto &coords = pi_pack.GetCoords(b);
          const Real flux_div =
              parthenon::Update::FluxDivHelper(idx, k, j, i, ndim, coords, pi);
          local_sum += -dt * coords.CellVolume(k, j, i) * flux_div;
        },
        step);
    return step;
  };

  const int pi0_mode0 = EMIndex(0, kCompA0);
  const int pix_mode0 = EMIndex(0, kCompAX);
  const int piy_mode0 = EMIndex(0, kCompAY);
  const int piz_mode0 = EMIndex(0, kCompAZ);
  const int piw_mode0 = EMIndex(0, kCompAW);

  const Real divpi0_step =
      accumulate_transport(pi0_mode0, "Modes4DAccumulateEMDivPi0Mode0");
  const Real divpix_step =
      accumulate_transport(pix_mode0, "Modes4DAccumulateEMDivPixMode0");
  const Real divpiy_step =
      accumulate_transport(piy_mode0, "Modes4DAccumulateEMDivPiyMode0");
  const Real divpiz_step =
      accumulate_transport(piz_mode0, "Modes4DAccumulateEMDivPizMode0");
  const Real divpiw_step =
      accumulate_transport(piw_mode0, "Modes4DAccumulateEMDivPiwMode0");

  auto *diag_divpi0_mode0 = modes_pkg->MutableParam<double>("diag/int_divpi0_mode0");
  auto *diag_divpix_mode0 = modes_pkg->MutableParam<double>("diag/int_divpix_mode0");
  auto *diag_divpiy_mode0 = modes_pkg->MutableParam<double>("diag/int_divpiy_mode0");
  auto *diag_divpiz_mode0 = modes_pkg->MutableParam<double>("diag/int_divpiz_mode0");
  auto *diag_divpiw_mode0 = modes_pkg->MutableParam<double>("diag/int_divpiw_mode0");
  *diag_divpi0_mode0 += divpi0_step;
  *diag_divpix_mode0 += divpix_step;
  *diag_divpiy_mode0 += divpiy_step;
  *diag_divpiz_mode0 += divpiz_step;
  *diag_divpiw_mode0 += divpiw_step;
  return TaskStatus::complete;
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
  const Real em_source_mass_gain = modes_pkg->Param<double>("em4d/source_mass_gain");
  const Real em_source_current_gain = modes_pkg->Param<double>("em4d/source_current_gain");
  const Real em_source_damping_gain = modes_pkg->Param<double>("em4d/source_damping_gain");
  const Real em_source_timelike_gain = modes_pkg->Param<double>("em4d/source_timelike_gain");
  const Real em_source_gauge_gain = modes_pkg->Param<double>("em4d/source_gauge_gain");
  const Real em_source_mode0_even_bridge_gain =
      modes_pkg->Param<double>("em4d/source_mode0_even_bridge_gain");
  const bool response_w0_enable = modes_pkg->Param<bool>("em4d/response_w0_enable");
  const Real response_w0_drive_gain = modes_pkg->Param<double>("em4d/response_w0_drive_gain");
  const Real response_w0_drive_from_bridge_gain =
      modes_pkg->Param<double>("em4d/response_w0_drive_from_bridge_gain");
  const std::string response_w0_drive_from_bridge_channel =
      modes_pkg->Param<std::string>("em4d/response_w0_drive_from_bridge_channel");
  const bool response_w0_drive_from_bridge_abs =
      modes_pkg->Param<bool>("em4d/response_w0_drive_from_bridge_abs");
  const Real response_w0_drive_from_parity_gain =
      modes_pkg->Param<double>("em4d/response_w0_drive_from_parity_gain");
  const Real response_w0_drive_from_work_gain =
      modes_pkg->Param<double>("em4d/response_w0_drive_from_work_gain");
  const bool response_w0_drive_from_work_abs =
      modes_pkg->Param<bool>("em4d/response_w0_drive_from_work_abs");
  const Real response_w0_drive_from_leak_gain =
      modes_pkg->Param<double>("em4d/response_w0_drive_from_leak_gain");
  const bool response_w0_drive_from_leak_abs =
      modes_pkg->Param<bool>("em4d/response_w0_drive_from_leak_abs");
  const Real response_w0_bias = modes_pkg->Param<double>("em4d/response_w0_bias");
  const Real response_w0_init = modes_pkg->Param<double>("em4d/response_w0_init");
  const Real response_w0_mass = modes_pkg->Param<double>("em4d/response_w0_mass");
  const Real response_w0_stiffness = modes_pkg->Param<double>("em4d/response_w0_stiffness");
  const Real response_w0_damping = modes_pkg->Param<double>("em4d/response_w0_damping");
  const Real response_w0_bridge_gain = modes_pkg->Param<double>("em4d/response_w0_bridge_gain");
  const Real response_w0_velocity_bridge_gain =
      modes_pkg->Param<double>("em4d/response_w0_velocity_bridge_gain");
  const Real response_w0_max_abs = modes_pkg->Param<double>("em4d/response_w0_max_abs");
  const bool response_w0_projection_enable =
      modes_pkg->Param<bool>("em4d/response_w0_projection_enable");
  const Real response_w0_projection_gain =
      modes_pkg->Param<double>("em4d/response_w0_projection_gain");
  const Real response_w0_projection_max_abs =
      modes_pkg->Param<double>("em4d/response_w0_projection_max_abs");
  const bool response_w0_geometry_shift_enable =
      modes_pkg->Param<bool>("em4d/response_w0_geometry_shift_enable");
  const Real response_w0_geometry_shift_gain =
      modes_pkg->Param<double>("em4d/response_w0_geometry_shift_gain");
  const bool response_w0_geometry_shift_include_constant =
      modes_pkg->Param<bool>("em4d/response_w0_geometry_shift_include_constant");
  const bool response_w0_lambda_shift_enable =
      modes_pkg->Param<bool>("em4d/response_w0_lambda_shift_enable");
  const Real response_w0_lambda_shift_gain =
      modes_pkg->Param<double>("em4d/response_w0_lambda_shift_gain");
  const Real response_w0_lambda_shift_max_frac =
      modes_pkg->Param<double>("em4d/response_w0_lambda_shift_max_frac");
  const bool response_w0_lambda_shift_apply_mass =
      modes_pkg->Param<bool>("em4d/response_w0_lambda_shift_apply_mass");
  const bool response_lambda_enable =
      modes_pkg->Param<bool>("em4d/response_lambda_enable");
  const Real response_lambda_drive_gain =
      modes_pkg->Param<double>("em4d/response_lambda_drive_gain");
  const Real response_lambda_drive_from_bridge_gain =
      modes_pkg->Param<double>("em4d/response_lambda_drive_from_bridge_gain");
  const std::string response_lambda_drive_from_bridge_channel =
      modes_pkg->Param<std::string>("em4d/response_lambda_drive_from_bridge_channel");
  const bool response_lambda_drive_from_bridge_abs =
      modes_pkg->Param<bool>("em4d/response_lambda_drive_from_bridge_abs");
  const Real response_lambda_drive_from_work_gain =
      modes_pkg->Param<double>("em4d/response_lambda_drive_from_work_gain");
  const bool response_lambda_drive_from_work_abs =
      modes_pkg->Param<bool>("em4d/response_lambda_drive_from_work_abs");
  const Real response_lambda_drive_from_leak_gain =
      modes_pkg->Param<double>("em4d/response_lambda_drive_from_leak_gain");
  const bool response_lambda_drive_from_leak_abs =
      modes_pkg->Param<bool>("em4d/response_lambda_drive_from_leak_abs");
  const Real response_lambda_bias =
      modes_pkg->Param<double>("em4d/response_lambda_bias");
  const Real response_lambda_init_fraction =
      modes_pkg->Param<double>("em4d/response_lambda_init_fraction");
  const Real response_lambda_mass =
      modes_pkg->Param<double>("em4d/response_lambda_mass");
  const Real response_lambda_stiffness =
      modes_pkg->Param<double>("em4d/response_lambda_stiffness");
  const Real response_lambda_damping =
      modes_pkg->Param<double>("em4d/response_lambda_damping");
  const Real response_lambda_max_abs_frac =
      modes_pkg->Param<double>("em4d/response_lambda_max_abs_frac");
  const bool response_lambda_apply_mass =
      modes_pkg->Param<bool>("em4d/response_lambda_apply_mass");
  const bool use_conservative_transport =
      modes_pkg->Param<bool>("em4d/use_conservative_transport");
  const Real qom_ion = modes_pkg->Param<double>("plasma4d/qom_ion");
  const Real qom_electron = modes_pkg->Param<double>("plasma4d/qom_electron");
  const Real force_source_gain = modes_pkg->Param<double>("plasma4d/force_source_gain");
  const Real momw_source_gain = modes_pkg->Param<double>("plasma4d/momw_source_gain");
  const Real momw_pressure_source_gain =
      modes_pkg->Param<double>("plasma4d/momw_pressure_source_gain");
  const Real momw_damping = modes_pkg->Param<double>("plasma4d/momw_damping");
  const Real w_flux_source_gain = modes_pkg->Param<double>("plasma4d/w_flux_source_gain");
  const Real rho_floor = modes_pkg->Param<double>("plasma4d/rho_floor");
  const Real energy_source_gain = modes_pkg->Param<double>("plasma4d/energy_source_gain");
  const Real energy_floor = modes_pkg->Param<double>("plasma4d/energy_floor");
  const Real gamma = modes_pkg->Param<double>("plasma4d/gamma");
  const Real gm1 = gamma - 1.0;
  const Real pressure_floor = modes_pkg->Param<double>("plasma4d/pressure_floor");
  const Real c2 = c_wave * c_wave;
  const Real inv_lambda2_base = 1.0 / (lambda * lambda);
  const auto &tables = modes_pkg->Param<ModeTables>("mode_tables");
  const auto &weights = tables.Weights();
  const auto &mass_squared = tables.MassSquared();

  Real response_w0 = modes_pkg->Param<double>("diag/response_w0");
  Real response_w0_dot = modes_pkg->Param<double>("diag/response_w0_dot");
  bool response_w0_initialized = modes_pkg->Param<bool>("diag/response_w0_initialized");
  Real response_w0_drive = 0.0;
  Real response_w0_drive_reservoir = 0.0;
  Real response_w0_drive_bridge = 0.0;
  Real response_w0_drive_parity = 0.0;
  Real response_w0_drive_work = 0.0;
  Real response_w0_drive_leak = 0.0;
  Real response_w0_force = 0.0;
  Real response_w0_energy = 0.0;
  Real response_w0_projection_center = 0.0;
  if (response_w0_enable && !response_w0_initialized) {
    response_w0 = response_w0_init;
    response_w0_dot = 0.0;
    response_w0_initialized = true;
  } else if (!response_w0_enable) {
    response_w0 = 0.0;
    response_w0_dot = 0.0;
    response_w0_initialized = false;
  }
  Real response_lambda_fraction =
      modes_pkg->Param<double>("diag/response_lambda_fraction");
  Real response_lambda_dot = modes_pkg->Param<double>("diag/response_lambda_dot");
  bool response_lambda_initialized =
      modes_pkg->Param<bool>("diag/response_lambda_initialized");
  Real response_lambda_drive = 0.0;
  Real response_lambda_drive_reservoir = 0.0;
  Real response_lambda_drive_bridge = 0.0;
  Real response_lambda_drive_work = 0.0;
  Real response_lambda_drive_leak = 0.0;
  Real response_lambda_force = 0.0;
  Real response_lambda_energy = 0.0;
  if (response_lambda_enable && !response_lambda_initialized) {
    response_lambda_fraction = response_lambda_init_fraction;
    response_lambda_dot = 0.0;
    response_lambda_initialized = true;
  } else if (!response_lambda_enable) {
    response_lambda_fraction = 0.0;
    response_lambda_dot = 0.0;
    response_lambda_initialized = false;
  }
  if (response_lambda_max_abs_frac > 0.0) {
    response_lambda_fraction = std::clamp(response_lambda_fraction, -response_lambda_max_abs_frac,
                                          response_lambda_max_abs_frac);
  }
  const Real response_w0_effective =
      response_w0 + (response_w0_velocity_bridge_gain * response_w0_dot);
  Real response_w0_lambda_fraction_from_w0 = 0.0;
  if (response_w0_enable && response_w0_lambda_shift_enable &&
      (response_w0_lambda_shift_gain != 0.0)) {
    response_w0_lambda_fraction_from_w0 =
        response_w0_lambda_shift_gain * response_w0_effective;
    if (response_w0_lambda_shift_max_frac > 0.0) {
      response_w0_lambda_fraction_from_w0 =
          std::clamp(response_w0_lambda_fraction_from_w0,
                     -response_w0_lambda_shift_max_frac,
                     response_w0_lambda_shift_max_frac);
    }
  }
  const Real response_lambda_fraction_total =
      response_w0_lambda_fraction_from_w0 + response_lambda_fraction;
  Real response_w0_lambda = lambda * (1.0 + response_lambda_fraction_total);
  response_w0_lambda = std::max(response_w0_lambda, 1.0e-6 * lambda);
  const Real inv_lambda2 = 1.0 / (response_w0_lambda * response_w0_lambda);
  const Real inv_lambda3 = inv_lambda2 / response_w0_lambda;
  const Real inv_lambda4 = inv_lambda2 * inv_lambda2;
  const Real inv_lambda_root2 = std::sqrt(2.0) / response_w0_lambda;
  const bool apply_lambda_mass_scale =
      (response_w0_lambda_shift_apply_mass &&
       (response_w0_lambda_fraction_from_w0 != 0.0)) ||
      (response_lambda_apply_mass && (response_lambda_fraction != 0.0));
  const Real response_w0_lambda_mass_scale =
      apply_lambda_mass_scale ? (inv_lambda2 / inv_lambda2_base) : 1.0;
  const bool response_w0_geometry_shift_active =
      response_w0_enable && response_w0_geometry_shift_enable &&
      (response_w0_geometry_shift_gain != 0.0);
  const Real response_w0_geometry_linear_coeff =
      response_w0_geometry_shift_active
          ? (-std::sqrt(2.0) * response_w0_geometry_shift_gain * response_w0 * inv_lambda3)
          : 0.0;
  const Real response_w0_geometry_diag_coeff =
      (response_w0_geometry_shift_active && response_w0_geometry_shift_include_constant)
          ? (response_w0_geometry_shift_gain * response_w0 * response_w0 * inv_lambda4)
          : 0.0;

  Real diag_jw_ew_step = 0.0;
  Real diag_ja_ea_step = 0.0;
  Real diag_s_leak_step = 0.0;
  Real diag_s_leak_abs_step = 0.0;
  Real diag_cont_local_l1_step = 0.0;
  Real diag_cont_local_l2_step = 0.0;
  Real diag_cont_local_max_abs_step = 0.0;
  Real diag_cont_mode0_l1_step = 0.0;
  Real diag_cont_mode0_l2_step = 0.0;
  Real diag_cont_mode0_max_abs_step = 0.0;
  Real diag_srcmomx_mode0_step = 0.0;
  Real diag_srcmomy_mode0_step = 0.0;
  Real diag_srcmomz_mode0_step = 0.0;
  Real diag_srcmomw_mode0_step = 0.0;
  Real diag_srcenergy_mode0_step = 0.0;
  Real diag_srcpi0_mode0_step = 0.0;
  Real diag_srcpix_mode0_step = 0.0;
  Real diag_srcpiy_mode0_step = 0.0;
  Real diag_srcpiz_mode0_step = 0.0;
  Real diag_srcpiw_mode0_step = 0.0;
  Real diag_src_timelike_a0_from_piw_step = 0.0;
  Real diag_src_timelike_a0_from_piw_abs_step = 0.0;
  Real diag_src_timelike_aw_from_pi0_step = 0.0;
  Real diag_src_timelike_aw_from_pi0_abs_step = 0.0;
  Real diag_src_em_laplacian_abs_step = 0.0;
  Real diag_src_em_mass_abs_step = 0.0;
  Real diag_src_em_current_abs_step = 0.0;
  Real diag_src_em_damping_abs_step = 0.0;
  Real diag_src_em_spatial_mixed_abs_step = 0.0;
  Real diag_src_em_timelike_abs_step = 0.0;
  Real diag_src_em_geometry_shift_abs_step = 0.0;
  Real diag_src_em_gauge_step = 0.0;
  Real diag_src_em_gauge_abs_step = 0.0;
  Real diag_src_em_gauge_mode0_step = 0.0;
  Real diag_src_em_gauge_mode0_abs_step = 0.0;
  Real diag_rhs_a0_mode0_lap_step = 0.0;
  Real diag_rhs_a0_mode0_mass_step = 0.0;
  Real diag_rhs_a0_mode0_current_step = 0.0;
  Real diag_rhs_a0_mode0_damping_step = 0.0;
  Real diag_rhs_a0_mode0_timelike_step = 0.0;
  Real diag_rhs_a0_mode0_gauge_step = 0.0;
  Real diag_rhs_a0_mode0_total_step = 0.0;
  Real diag_rhs_ay_mode0_lap_step = 0.0;
  Real diag_rhs_ay_mode0_mass_step = 0.0;
  Real diag_rhs_ay_mode0_current_step = 0.0;
  Real diag_rhs_ay_mode0_damping_step = 0.0;
  Real diag_rhs_ay_mode0_spatial_mixed_step = 0.0;
  Real diag_rhs_ay_mode0_even_bridge_step = 0.0;
  Real diag_rhs_ay_mode0_even_bridge_abs_step = 0.0;
  Real diag_rhs_ay_mode0_geometry_shift_step = 0.0;
  Real diag_rhs_ay_mode0_geometry_shift_abs_step = 0.0;
  Real diag_bridge_power_mode0_step = 0.0;
  Real diag_bridge_power_mode0_abs_step = 0.0;
  Real diag_geometry_shift_power_mode0_step = 0.0;
  Real diag_geometry_shift_power_mode0_abs_step = 0.0;
  Real diag_rhs_aw_mode2_even_bridge_step = 0.0;
  Real diag_rhs_aw_mode2_even_bridge_abs_step = 0.0;
  Real diag_bridge_power_mode2_step = 0.0;
  Real diag_bridge_power_mode2_abs_step = 0.0;
  Real diag_bridge_power_sum_step = 0.0;
  Real diag_bridge_power_sum_abs_step = 0.0;
  Real diag_rhs_ay_mode0_w0_response_step = 0.0;
  Real diag_rhs_ay_mode0_w0_response_abs_step = 0.0;
  Real diag_rhs_aw_mode1_w0_response_step = 0.0;
  Real diag_rhs_aw_mode1_w0_response_abs_step = 0.0;
  Real diag_response_bridge_power_mode0_step = 0.0;
  Real diag_response_bridge_power_mode0_abs_step = 0.0;
  Real diag_response_bridge_power_mode1_step = 0.0;
  Real diag_response_bridge_power_mode1_abs_step = 0.0;
  Real diag_response_bridge_power_sum_step = 0.0;
  Real diag_response_bridge_power_sum_abs_step = 0.0;
  Real diag_response_w0_power_drive_step = 0.0;
  Real diag_response_w0_power_stiffness_step = 0.0;
  Real diag_response_w0_power_damping_step = 0.0;
  Real diag_response_w0_power_net_step = 0.0;
  Real diag_response_lambda_power_drive_step = 0.0;
  Real diag_response_lambda_power_stiffness_step = 0.0;
  Real diag_response_lambda_power_damping_step = 0.0;
  Real diag_response_lambda_power_net_step = 0.0;
  Real diag_rhs_ay_mode0_total_step = 0.0;
  Real diag_rhs_aw_mode0_lap_step = 0.0;
  Real diag_rhs_aw_mode0_current_step = 0.0;
  Real diag_rhs_aw_mode0_damping_step = 0.0;
  Real diag_rhs_aw_mode0_spatial_mixed_step = 0.0;
  Real diag_rhs_aw_mode0_timelike_step = 0.0;
  Real diag_rhs_aw_mode0_total_step = 0.0;
  Real diag_gauge_l1_step = 0.0;
  Real diag_gauge_l2_step = 0.0;
  Real diag_gauge_max_abs_step = 0.0;
  Real diag_gauge_mode0_l1_step = 0.0;
  Real diag_gauge_mode0_l2_step = 0.0;
  Real diag_gauge_mode0_max_abs_step = 0.0;
  Real response_drive_energy_integral_step = 0.0;
  Real response_drive_odd_energy_integral_step = 0.0;
  Real response_drive_even_energy_integral_step = 0.0;
  Real response_drive_volume_step = 0.0;

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
    std::vector<Real> pressure_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> energy_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> momx_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> momy_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> momz_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> pressure_modes(n_modes, 0.0);
    std::vector<Real> momw_pressure_rhs_modes(n_modes, 0.0);
    std::vector<Real> momw_modes_new(n_modes, 0.0);
    std::vector<Real> momx_modes_new(n_modes, 0.0);
    std::vector<Real> momy_modes_new(n_modes, 0.0);
    std::vector<Real> momz_modes_new(n_modes, 0.0);
    std::vector<Real> energy_modes_new(n_modes, 0.0);
    std::vector<Real> energy_nodes_new(tables.NumQuadrature(), 0.0);
    std::vector<Real> wflux_momx_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> wflux_momy_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> wflux_momz_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> wflux_momw_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> wflux_energy_nodes(tables.NumQuadrature(), 0.0);
    std::vector<Real> wflux_momx_modes(n_modes, 0.0);
    std::vector<Real> wflux_momy_modes(n_modes, 0.0);
    std::vector<Real> wflux_momz_modes(n_modes, 0.0);
    std::vector<Real> wflux_momw_modes(n_modes, 0.0);
    std::vector<Real> wflux_energy_modes(n_modes, 0.0);
    std::vector<Real> charge_modes_old(n_modes, 0.0);
    const int pi0_mode0_idx = EMIndex(0, kCompA0);
    const int pix_mode0_idx = EMIndex(0, kCompAX);
    const int piy_mode0_idx = EMIndex(0, kCompAY);
    const int piz_mode0_idx = EMIndex(0, kCompAZ);
    const int piw_mode0_idx = EMIndex(0, kCompAW);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real cell_volume = coords.CellVolume(k, j, i);
          for (int n = 0; n < n_modes; ++n) {
            const int ion_base = PlasmaIndex(0, n, 0, n_modes);
            const int ele_base = PlasmaIndex(1, n, 0, n_modes);
            const Real ion_rho = plasma(ion_base + kPlasmaRho, k, j, i);
            const Real ele_rho = plasma(ele_base + kPlasmaRho, k, j, i);
            charge_modes_old[n] = (qom_ion * ion_rho) + (qom_electron * ele_rho);
          }

          const int ion_mode0 = PlasmaIndex(0, 0, 0, n_modes);
          const int ele_mode0 = PlasmaIndex(1, 0, 0, n_modes);
          const Real momx_mode0_old = plasma(ion_mode0 + kPlasmaMomX, k, j, i) +
                                      plasma(ele_mode0 + kPlasmaMomX, k, j, i);
          const Real momy_mode0_old = plasma(ion_mode0 + kPlasmaMomY, k, j, i) +
                                      plasma(ele_mode0 + kPlasmaMomY, k, j, i);
          const Real momz_mode0_old = plasma(ion_mode0 + kPlasmaMomZ, k, j, i) +
                                      plasma(ele_mode0 + kPlasmaMomZ, k, j, i);
          const Real momw_mode0_old = plasma(ion_mode0 + kPlasmaMomW, k, j, i) +
                                      plasma(ele_mode0 + kPlasmaMomW, k, j, i);
          const Real energy_mode0_old = plasma(ion_mode0 + kPlasmaEnergy, k, j, i) +
                                        plasma(ele_mode0 + kPlasmaEnergy, k, j, i);
          const Real pi0_mode0_old = pi_old(pi0_mode0_idx, k, j, i);
          const Real pix_mode0_old = pi_old(pix_mode0_idx, k, j, i);
          const Real piy_mode0_old = pi_old(piy_mode0_idx, k, j, i);
          const Real piz_mode0_old = pi_old(piz_mode0_idx, k, j, i);
          const Real piw_mode0_old = pi_old(piw_mode0_idx, k, j, i);

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
              const Real kinetic = 0.5 *
                                   ((momx_node * momx_node) + (momy_node * momy_node) +
                                    (momz_node * momz_node) + (momw_node * momw_node)) /
                                   rho_safe;
              const Real pressure_node =
                  std::max(pressure_floor, gm1 * std::max(energy_node - kinetic, 0.0));
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
              pressure_nodes[q] = pressure_node;
              energy_nodes[q] = energy_node;
              momx_nodes_new[q] = momx_node + (dt * rhs_momx);
              momy_nodes_new[q] = momy_node + (dt * rhs_momy);
              momz_nodes_new[q] = momz_node + (dt * rhs_momz);
              energy_nodes_new[q] = std::max(energy_floor, energy_node + (dt * rhs_energy));
            }

            if (momw_pressure_source_gain != 0.0) {
              for (int n = 0; n < n_modes; ++n) {
                Real projected_pressure = 0.0;
                for (int q = 0; q < tables.NumQuadrature(); ++q) {
                  projected_pressure += weights[q] * tables.Phi(n, q) * pressure_nodes[q];
                }
                pressure_modes[n] = projected_pressure;
              }
              for (int n = 0; n < n_modes; ++n) {
                // Two-fluid transverse momentum closure:
                // d_t (rho v_w) includes -<phi_n, d_w p>_Z; evaluate by ID-3.
                momw_pressure_rhs_modes[n] =
                    -momw_pressure_source_gain * tables.ApplyID3Raising(pressure_modes, n);
              }
            } else {
              for (int n = 0; n < n_modes; ++n) {
                momw_pressure_rhs_modes[n] = 0.0;
              }
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
              momw_modes_new[n] = projected_momw + (dt * momw_pressure_rhs_modes[n]);
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

          if (w_flux_source_gain != 0.0) {
            // Optional two-fluid closure extension: apply projected w-flux couplings
            // to momentum/energy modes via the same leakage-style mode raising.
            for (int s = 0; s < kSpeciesCount; ++s) {
              for (int n = 0; n < n_modes; ++n) {
                const int base = PlasmaIndex(s, n, 0, n_modes);
                rho_modes[n] = plasma_new(base + kPlasmaRho, k, j, i);
                momx_modes[n] = plasma_new(base + kPlasmaMomX, k, j, i);
                momy_modes[n] = plasma_new(base + kPlasmaMomY, k, j, i);
                momz_modes[n] = plasma_new(base + kPlasmaMomZ, k, j, i);
                momw_modes[n] = plasma_new(base + kPlasmaMomW, k, j, i);
                energy_modes[n] = plasma_new(base + kPlasmaEnergy, k, j, i);
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
                const Real vw = momw_node / rho_safe;
                const Real kinetic = 0.5 *
                                     ((momx_node * momx_node) + (momy_node * momy_node) +
                                      (momz_node * momz_node) + (momw_node * momw_node)) /
                                     rho_safe;
                const Real pressure_node =
                    std::max(pressure_floor, gm1 * std::max(energy_node - kinetic, 0.0));

                wflux_momx_nodes[q] = momx_node * vw;
                wflux_momy_nodes[q] = momy_node * vw;
                wflux_momz_nodes[q] = momz_node * vw;
                wflux_momw_nodes[q] = (momw_node * vw) + pressure_node;
                wflux_energy_nodes[q] = (energy_node + pressure_node) * vw;
              }

              for (int n = 0; n < n_modes; ++n) {
                Real projected_wflux_momx = 0.0;
                Real projected_wflux_momy = 0.0;
                Real projected_wflux_momz = 0.0;
                Real projected_wflux_momw = 0.0;
                Real projected_wflux_energy = 0.0;
                for (int q = 0; q < tables.NumQuadrature(); ++q) {
                  const Real proj = weights[q] * tables.Phi(n, q);
                  projected_wflux_momx += proj * wflux_momx_nodes[q];
                  projected_wflux_momy += proj * wflux_momy_nodes[q];
                  projected_wflux_momz += proj * wflux_momz_nodes[q];
                  projected_wflux_momw += proj * wflux_momw_nodes[q];
                  projected_wflux_energy += proj * wflux_energy_nodes[q];
                }
                wflux_momx_modes[n] = projected_wflux_momx;
                wflux_momy_modes[n] = projected_wflux_momy;
                wflux_momz_modes[n] = projected_wflux_momz;
                wflux_momw_modes[n] = projected_wflux_momw;
                wflux_energy_modes[n] = projected_wflux_energy;
              }

              for (int n = 0; n < n_modes; ++n) {
                const int base = PlasmaIndex(s, n, 0, n_modes);
                const Real coupling =
                    (n + 1 < n_modes)
                        ? (std::sqrt(2.0 * static_cast<Real>(n + 1)) / response_w0_lambda)
                        : 0.0;
                const Real leak_momx =
                    (n + 1 < n_modes) ? (-coupling * wflux_momx_modes[n + 1]) : 0.0;
                const Real leak_momy =
                    (n + 1 < n_modes) ? (-coupling * wflux_momy_modes[n + 1]) : 0.0;
                const Real leak_momz =
                    (n + 1 < n_modes) ? (-coupling * wflux_momz_modes[n + 1]) : 0.0;
                const Real leak_momw =
                    (n + 1 < n_modes) ? (-coupling * wflux_momw_modes[n + 1]) : 0.0;
                const Real leak_energy =
                    (n + 1 < n_modes) ? (-coupling * wflux_energy_modes[n + 1]) : 0.0;
                plasma_new(base + kPlasmaMomX, k, j, i) += dt * w_flux_source_gain * leak_momx;
                plasma_new(base + kPlasmaMomY, k, j, i) += dt * w_flux_source_gain * leak_momy;
                plasma_new(base + kPlasmaMomZ, k, j, i) += dt * w_flux_source_gain * leak_momz;
                plasma_new(base + kPlasmaMomW, k, j, i) += dt * w_flux_source_gain * leak_momw;
                plasma_new(base + kPlasmaEnergy, k, j, i) = std::max(
                    energy_floor,
                    plasma_new(base + kPlasmaEnergy, k, j, i) +
                        dt * w_flux_source_gain * leak_energy);
              }
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
                const Real coupling =
                    std::sqrt(2.0 * static_cast<Real>(n + 1)) / response_w0_lambda;
                leak_rhs = -coupling * plasma_new(base_np1 + kPlasmaMomW, k, j, i);
              }

              const Real rhs_rho = leak_rhs;
              const Real rho_old = plasma(base + kPlasmaRho, k, j, i);
              plasma_new(base + kPlasmaRho, k, j, i) =
                  std::max(rho_floor, rho_old + (dt * rhs_rho));
            }
          }

          const Real momx_mode0_new = plasma_new(ion_mode0 + kPlasmaMomX, k, j, i) +
                                      plasma_new(ele_mode0 + kPlasmaMomX, k, j, i);
          const Real momy_mode0_new = plasma_new(ion_mode0 + kPlasmaMomY, k, j, i) +
                                      plasma_new(ele_mode0 + kPlasmaMomY, k, j, i);
          const Real momz_mode0_new = plasma_new(ion_mode0 + kPlasmaMomZ, k, j, i) +
                                      plasma_new(ele_mode0 + kPlasmaMomZ, k, j, i);
          const Real momw_mode0_new = plasma_new(ion_mode0 + kPlasmaMomW, k, j, i) +
                                      plasma_new(ele_mode0 + kPlasmaMomW, k, j, i);
          const Real energy_mode0_new = plasma_new(ion_mode0 + kPlasmaEnergy, k, j, i) +
                                        plasma_new(ele_mode0 + kPlasmaEnergy, k, j, i);
          diag_srcmomx_mode0_step += cell_volume * (momx_mode0_new - momx_mode0_old);
          diag_srcmomy_mode0_step += cell_volume * (momy_mode0_new - momy_mode0_old);
          diag_srcmomz_mode0_step += cell_volume * (momz_mode0_new - momz_mode0_old);
          diag_srcmomw_mode0_step += cell_volume * (momw_mode0_new - momw_mode0_old);
          diag_srcenergy_mode0_step += cell_volume * (energy_mode0_new - energy_mode0_old);

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

          for (int n = 0; n < n_modes; ++n) {
            const Real coupling =
                (n + 1 < n_modes)
                    ? (std::sqrt(2.0 * static_cast<Real>(n + 1)) / response_w0_lambda)
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

          Real bridge_power_mode0_cell = 0.0;
          Real response_bridge_power_mode0_cell = 0.0;
          for (int n = 0; n < n_modes; ++n) {
            const int idx_a0 = EMIndex(n, kCompA0);
            const int idx_ax = EMIndex(n, kCompAX);
            const int idx_ay = EMIndex(n, kCompAY);
            const int idx_az = EMIndex(n, kCompAZ);
            const int idx_aw = EMIndex(n, kCompAW);

            const Real coupling_coeff =
                (n > 0) ? (std::sqrt(2.0 * static_cast<Real>(n)) / response_w0_lambda)
                        : 0.0;
            Real mixed_grad_aw_x = 0.0;
            Real mixed_grad_aw_y = 0.0;
            Real mixed_grad_aw_z = 0.0;
            Real mixed_dt_aw_prev = 0.0;
            if ((n > 0) && !use_conservative_transport) {
              const int idx_aw_prev = EMIndex(n - 1, kCompAW);
              mixed_grad_aw_x = gradient_x(a_old, idx_aw_prev, k, j, i);
              mixed_grad_aw_y = gradient_y(a_old, idx_aw_prev, k, j, i);
              mixed_grad_aw_z = gradient_z(a_old, idx_aw_prev, k, j, i);
            }
            if (n > 0) {
              const int idx_aw_prev = EMIndex(n - 1, kCompAW);
              // Time-like mixed coupling from +sqrt(2n)/lambda * d^0 a_w^(n-1).
              mixed_dt_aw_prev = pi_old(idx_aw_prev, k, j, i);
            }

            Real mixed_div_a_next = 0.0;
            const Real coupling_raise =
                (n + 1 < n_modes)
                    ? (std::sqrt(2.0 * static_cast<Real>(n + 1)) / response_w0_lambda)
                    : 0.0;
            const int idx_aw_next = (n + 1 < n_modes) ? EMIndex(n + 1, kCompAW) : -1;
            if ((n + 1 < n_modes) && !use_conservative_transport) {
              const int idx_ax_next = EMIndex(n + 1, kCompAX);
              const int idx_ay_next = EMIndex(n + 1, kCompAY);
              const int idx_az_next = EMIndex(n + 1, kCompAZ);
              mixed_div_a_next = coupling_raise *
                                 (gradient_x(a_old, idx_ax_next, k, j, i) +
                                  gradient_y(a_old, idx_ay_next, k, j, i) +
                                  gradient_z(a_old, idx_az_next, k, j, i));
            }
            Real mixed_dt_a0_next = 0.0;
            if (n + 1 < n_modes) {
              const int idx_a0_next = EMIndex(n + 1, kCompA0);
              // Time-like piece of -sqrt(2(n+1))/lambda * d_mu a^{mu,(n+1)}.
              mixed_dt_a0_next = pi_old(idx_a0_next, k, j, i);
            }

            // Mode-projected Lorenz-like gauge residual:
            // G^(n) = d_t a0^(n) + div a^(n) + sqrt(2(n+1))/lambda * a_w^(n+1).
            const Real gauge_div_spatial = gradient_x(a_old, idx_ax, k, j, i) +
                                           gradient_y(a_old, idx_ay, k, j, i) +
                                           gradient_z(a_old, idx_az, k, j, i);
            const Real aw_next =
                (idx_aw_next >= 0) ? a_old(idx_aw_next, k, j, i) : 0.0;
            const Real gauge_residual =
                pi_old(idx_a0, k, j, i) + gauge_div_spatial + (coupling_raise * aw_next);
            const Real abs_gauge_residual = std::abs(gauge_residual);
            diag_gauge_l1_step += cell_volume * abs_gauge_residual;
            diag_gauge_l2_step += cell_volume * gauge_residual * gauge_residual;
            diag_gauge_max_abs_step =
                std::max(diag_gauge_max_abs_step, abs_gauge_residual);
            if (n == 0) {
              diag_gauge_mode0_l1_step += cell_volume * abs_gauge_residual;
              diag_gauge_mode0_l2_step += cell_volume * gauge_residual * gauge_residual;
              diag_gauge_mode0_max_abs_step =
                  std::max(diag_gauge_mode0_max_abs_step, abs_gauge_residual);
            }

            const Real lap_a0 = use_conservative_transport
                                    ? 0.0
                                    : (c2 * laplacian(a_old, idx_a0, k, j, i));
            const Real lap_ax = use_conservative_transport
                                    ? 0.0
                                    : (c2 * laplacian(a_old, idx_ax, k, j, i));
            const Real lap_ay = use_conservative_transport
                                    ? 0.0
                                    : (c2 * laplacian(a_old, idx_ay, k, j, i));
            const Real lap_az = use_conservative_transport
                                    ? 0.0
                                    : (c2 * laplacian(a_old, idx_az, k, j, i));
            const Real lap_aw = use_conservative_transport
                                    ? 0.0
                                    : (c2 * laplacian(a_old, idx_aw, k, j, i));
            const Real src_a0_lap = lap_a0;
            const Real src_ax_lap = lap_ax;
            const Real src_ay_lap = lap_ay;
            const Real src_az_lap = lap_az;
            const Real src_aw_lap = lap_aw;

            const Real mass_squared_eff = mass_squared[n] * response_w0_lambda_mass_scale;
            const Real src_a0_mass =
                -em_source_mass_gain * (c2 * mass_squared_eff * a_old(idx_a0, k, j, i));
            const Real src_ax_mass =
                -em_source_mass_gain * (c2 * mass_squared_eff * a_old(idx_ax, k, j, i));
            const Real src_ay_mass =
                -em_source_mass_gain * (c2 * mass_squared_eff * a_old(idx_ay, k, j, i));
            const Real src_az_mass =
                -em_source_mass_gain * (c2 * mass_squared_eff * a_old(idx_az, k, j, i));
            Real src_a0_geometry_shift = 0.0;
            Real src_ax_geometry_shift = 0.0;
            Real src_ay_geometry_shift = 0.0;
            Real src_az_geometry_shift = 0.0;
            if (response_w0_geometry_shift_active) {
              const Real a0_prev =
                  (n > 0) ? a_old(EMIndex(n - 1, kCompA0), k, j, i) : 0.0;
              const Real a0_next =
                  (n + 1 < n_modes) ? a_old(EMIndex(n + 1, kCompA0), k, j, i) : 0.0;
              const Real ax_prev =
                  (n > 0) ? a_old(EMIndex(n - 1, kCompAX), k, j, i) : 0.0;
              const Real ax_next =
                  (n + 1 < n_modes) ? a_old(EMIndex(n + 1, kCompAX), k, j, i) : 0.0;
              const Real ay_prev =
                  (n > 0) ? a_old(EMIndex(n - 1, kCompAY), k, j, i) : 0.0;
              const Real ay_next =
                  (n + 1 < n_modes) ? a_old(EMIndex(n + 1, kCompAY), k, j, i) : 0.0;
              const Real az_prev =
                  (n > 0) ? a_old(EMIndex(n - 1, kCompAZ), k, j, i) : 0.0;
              const Real az_next =
                  (n + 1 < n_modes) ? a_old(EMIndex(n + 1, kCompAZ), k, j, i) : 0.0;

              const Real a0_shift_mix = (coupling_coeff * a0_prev) + (coupling_raise * a0_next);
              const Real ax_shift_mix = (coupling_coeff * ax_prev) + (coupling_raise * ax_next);
              const Real ay_shift_mix = (coupling_coeff * ay_prev) + (coupling_raise * ay_next);
              const Real az_shift_mix = (coupling_coeff * az_prev) + (coupling_raise * az_next);

              src_a0_geometry_shift = -em_source_mass_gain * c2 *
                                      ((response_w0_geometry_linear_coeff * a0_shift_mix) +
                                       (response_w0_geometry_diag_coeff *
                                        a_old(idx_a0, k, j, i)));
              src_ax_geometry_shift = -em_source_mass_gain * c2 *
                                      ((response_w0_geometry_linear_coeff * ax_shift_mix) +
                                       (response_w0_geometry_diag_coeff *
                                        a_old(idx_ax, k, j, i)));
              src_ay_geometry_shift = -em_source_mass_gain * c2 *
                                      ((response_w0_geometry_linear_coeff * ay_shift_mix) +
                                       (response_w0_geometry_diag_coeff *
                                        a_old(idx_ay, k, j, i)));
              src_az_geometry_shift = -em_source_mass_gain * c2 *
                                      ((response_w0_geometry_linear_coeff * az_shift_mix) +
                                       (response_w0_geometry_diag_coeff *
                                        a_old(idx_az, k, j, i)));
            }

            const Real src_a0_current = -em_source_current_gain * (mu0 * j0_modes[n]);
            const Real src_ax_current = -em_source_current_gain * (mu0 * jx_modes[n]);
            const Real src_ay_current = -em_source_current_gain * (mu0 * jy_modes[n]);
            const Real src_az_current = -em_source_current_gain * (mu0 * jz_modes[n]);
            const Real src_aw_current = -em_source_current_gain * (mu0 * jw_modes[n]);

            const Real src_a0_damping =
                -em_source_damping_gain * (damping * pi_old(idx_a0, k, j, i));
            const Real src_ax_damping =
                -em_source_damping_gain * (damping * pi_old(idx_ax, k, j, i));
            const Real src_ay_damping =
                -em_source_damping_gain * (damping * pi_old(idx_ay, k, j, i));
            const Real src_az_damping =
                -em_source_damping_gain * (damping * pi_old(idx_az, k, j, i));
            const Real src_aw_damping =
                -em_source_damping_gain * (damping * pi_old(idx_aw, k, j, i));
            const Real src_a0_gauge =
                -em_source_gauge_gain * (c2 * gauge_residual);

            const Real src_ax_spatial_mixed = c2 * coupling_coeff * mixed_grad_aw_x;
            const Real src_ay_spatial_mixed = c2 * coupling_coeff * mixed_grad_aw_y;
            const Real src_az_spatial_mixed = c2 * coupling_coeff * mixed_grad_aw_z;
            const Real src_aw_spatial_mixed = -(c2 * mixed_div_a_next);
            Real src_ay_w0_response = 0.0;
            if (response_w0_enable && (response_w0_bridge_gain != 0.0) &&
                (n + 1 < n_modes)) {
              // Dynamic center-shift response: w0 mixes adjacent Hermite modes
              // through an antisymmetric n <-> n+1 coupling.
              const int idx_aw_next = EMIndex(n + 1, kCompAW);
              src_ay_w0_response = response_w0_bridge_gain * response_w0_effective * c2 *
                                   coupling_raise *
                                   gradient_y(a_old, idx_aw_next, k, j, i);
            }
            Real src_ay_mode0_even_bridge = 0.0;
            if ((n == 0) && (n_modes > 2) && (em_source_mode0_even_bridge_gain != 0.0)) {
              // Optional explicit even-parity bridge: mode-2 A_w gradient into mode-0 A_y RHS.
              // Disabled by default (gain=0) so baseline behavior is unchanged.
              const int idx_aw_mode2 = EMIndex(2, kCompAW);
              const Real mode0_even_bridge_coeff =
                  std::sqrt(8.0) / (response_w0_lambda * response_w0_lambda);
              src_ay_mode0_even_bridge =
                  em_source_mode0_even_bridge_gain * c2 * mode0_even_bridge_coeff *
                  gradient_y(a_old, idx_aw_mode2, k, j, i);
            }
            Real src_aw_mode2_even_bridge = 0.0;
            if ((n == 2) && (em_source_mode0_even_bridge_gain != 0.0)) {
              // Paired antisymmetric companion term so bridge transfer is internal (0 <-> 2).
              const int idx_ay_mode0 = EMIndex(0, kCompAY);
              const Real mode0_even_bridge_coeff =
                  std::sqrt(8.0) / (response_w0_lambda * response_w0_lambda);
              src_aw_mode2_even_bridge =
                  -em_source_mode0_even_bridge_gain * c2 * mode0_even_bridge_coeff *
                  gradient_y(a_old, idx_ay_mode0, k, j, i);
            }
            Real src_aw_w0_response = 0.0;
            if (response_w0_enable && (response_w0_bridge_gain != 0.0) && (n > 0)) {
              const int idx_ay_prev = EMIndex(n - 1, kCompAY);
              src_aw_w0_response = -response_w0_bridge_gain * response_w0_effective * c2 *
                                   coupling_coeff *
                                   gradient_y(a_old, idx_ay_prev, k, j, i);
            }

            const Real rhs_a0_timelike_from_piw =
                -em_source_timelike_gain * (c2 * coupling_coeff * mixed_dt_aw_prev);
            const Real rhs_aw_timelike_from_pi0 =
                -em_source_timelike_gain * (c2 * coupling_raise * mixed_dt_a0_next);

            const Real rhs_a0 = src_a0_lap + src_a0_mass + src_a0_current +
                                src_a0_damping + src_a0_geometry_shift +
                                rhs_a0_timelike_from_piw +
                                src_a0_gauge;
            const Real rhs_ax = src_ax_lap + src_ax_mass + src_ax_geometry_shift +
                                src_ax_spatial_mixed +
                                src_ax_current + src_ax_damping;
            const Real rhs_ay = src_ay_lap + src_ay_mass + src_ay_geometry_shift +
                                src_ay_spatial_mixed +
                                src_ay_current + src_ay_damping +
                                src_ay_w0_response +
                                src_ay_mode0_even_bridge;
            const Real rhs_az = src_az_lap + src_az_mass + src_az_geometry_shift +
                                src_az_spatial_mixed +
                                src_az_current + src_az_damping;
            const Real rhs_aw = src_aw_lap + src_aw_spatial_mixed + src_aw_current +
                                src_aw_damping + rhs_aw_timelike_from_pi0 +
                                src_aw_w0_response +
                                src_aw_mode2_even_bridge;

            // Track the time-like mixed couplings explicitly so conservative-path runs
            // can be diagnosed before tightening behavioral gates.
            diag_src_timelike_a0_from_piw_step +=
                dt * cell_volume * rhs_a0_timelike_from_piw;
            diag_src_timelike_a0_from_piw_abs_step +=
                dt * cell_volume * std::abs(rhs_a0_timelike_from_piw);
            diag_src_timelike_aw_from_pi0_step +=
                dt * cell_volume * rhs_aw_timelike_from_pi0;
            diag_src_timelike_aw_from_pi0_abs_step +=
                dt * cell_volume * std::abs(rhs_aw_timelike_from_pi0);
            diag_src_em_laplacian_abs_step +=
                dt * cell_volume *
                (std::abs(src_a0_lap) + std::abs(src_ax_lap) + std::abs(src_ay_lap) +
                 std::abs(src_az_lap) + std::abs(src_aw_lap));
            diag_src_em_mass_abs_step +=
                dt * cell_volume *
                (std::abs(src_a0_mass) + std::abs(src_ax_mass) + std::abs(src_ay_mass) +
                 std::abs(src_az_mass));
            diag_src_em_current_abs_step +=
                dt * cell_volume *
                (std::abs(src_a0_current) + std::abs(src_ax_current) +
                 std::abs(src_ay_current) + std::abs(src_az_current) +
                 std::abs(src_aw_current));
            diag_src_em_damping_abs_step +=
                dt * cell_volume *
                (std::abs(src_a0_damping) + std::abs(src_ax_damping) +
                 std::abs(src_ay_damping) + std::abs(src_az_damping) +
                 std::abs(src_aw_damping));
            diag_src_em_geometry_shift_abs_step +=
                dt * cell_volume *
                (std::abs(src_a0_geometry_shift) + std::abs(src_ax_geometry_shift) +
                 std::abs(src_ay_geometry_shift) + std::abs(src_az_geometry_shift));
            diag_src_em_gauge_step +=
                dt * cell_volume * src_a0_gauge;
            diag_src_em_gauge_abs_step +=
                dt * cell_volume * std::abs(src_a0_gauge);
            diag_src_em_spatial_mixed_abs_step +=
                dt * cell_volume *
                (std::abs(src_ax_spatial_mixed) + std::abs(src_ay_spatial_mixed) +
                 std::abs(src_az_spatial_mixed) + std::abs(src_aw_spatial_mixed));
            diag_src_em_timelike_abs_step +=
                dt * cell_volume *
                (std::abs(rhs_a0_timelike_from_piw) + std::abs(rhs_aw_timelike_from_pi0));
            if (n == 0) {
              diag_src_em_gauge_mode0_step += dt * cell_volume * src_a0_gauge;
              diag_src_em_gauge_mode0_abs_step += dt * cell_volume * std::abs(src_a0_gauge);
              diag_rhs_a0_mode0_lap_step += dt * cell_volume * src_a0_lap;
              diag_rhs_a0_mode0_mass_step += dt * cell_volume * src_a0_mass;
              diag_rhs_a0_mode0_current_step += dt * cell_volume * src_a0_current;
              diag_rhs_a0_mode0_damping_step += dt * cell_volume * src_a0_damping;
              diag_rhs_a0_mode0_timelike_step +=
                  dt * cell_volume * rhs_a0_timelike_from_piw;
              diag_rhs_a0_mode0_gauge_step += dt * cell_volume * src_a0_gauge;
              diag_rhs_a0_mode0_total_step += dt * cell_volume * rhs_a0;
              diag_rhs_ay_mode0_lap_step += dt * cell_volume * src_ay_lap;
              diag_rhs_ay_mode0_mass_step += dt * cell_volume * src_ay_mass;
              diag_rhs_ay_mode0_current_step += dt * cell_volume * src_ay_current;
              diag_rhs_ay_mode0_damping_step += dt * cell_volume * src_ay_damping;
              diag_rhs_ay_mode0_spatial_mixed_step +=
                  dt * cell_volume * src_ay_spatial_mixed;
              diag_rhs_ay_mode0_even_bridge_step +=
                  dt * cell_volume * src_ay_mode0_even_bridge;
              diag_rhs_ay_mode0_even_bridge_abs_step +=
                  dt * cell_volume * std::abs(src_ay_mode0_even_bridge);
              diag_rhs_ay_mode0_geometry_shift_step +=
                  dt * cell_volume * src_ay_geometry_shift;
              diag_rhs_ay_mode0_geometry_shift_abs_step +=
                  dt * cell_volume * std::abs(src_ay_geometry_shift);
              diag_rhs_ay_mode0_w0_response_step +=
                  dt * cell_volume * src_ay_w0_response;
              diag_rhs_ay_mode0_w0_response_abs_step +=
                  dt * cell_volume * std::abs(src_ay_w0_response);
              const Real ey_mode0 =
                  -pi_old(idx_ay, k, j, i) - gradient_y(a_old, idx_a0, k, j, i);
              const Real bridge_power_mode0 = -src_ay_mode0_even_bridge * ey_mode0;
              bridge_power_mode0_cell = bridge_power_mode0;
              diag_bridge_power_mode0_step += dt * cell_volume * bridge_power_mode0;
              diag_bridge_power_mode0_abs_step +=
                  dt * cell_volume * std::abs(bridge_power_mode0);
              const Real response_bridge_power_mode0 =
                  -src_ay_w0_response * ey_mode0;
              response_bridge_power_mode0_cell = response_bridge_power_mode0;
              diag_response_bridge_power_mode0_step +=
                  dt * cell_volume * response_bridge_power_mode0;
              diag_response_bridge_power_mode0_abs_step +=
                  dt * cell_volume * std::abs(response_bridge_power_mode0);
              const Real geometry_shift_power_mode0 =
                  -src_ay_geometry_shift * ey_mode0;
              diag_geometry_shift_power_mode0_step +=
                  dt * cell_volume * geometry_shift_power_mode0;
              diag_geometry_shift_power_mode0_abs_step +=
                  dt * cell_volume * std::abs(geometry_shift_power_mode0);
              diag_rhs_ay_mode0_total_step += dt * cell_volume * rhs_ay;
              diag_rhs_aw_mode0_lap_step += dt * cell_volume * src_aw_lap;
              diag_rhs_aw_mode0_current_step += dt * cell_volume * src_aw_current;
              diag_rhs_aw_mode0_damping_step += dt * cell_volume * src_aw_damping;
              diag_rhs_aw_mode0_spatial_mixed_step +=
                  dt * cell_volume * src_aw_spatial_mixed;
              diag_rhs_aw_mode0_timelike_step +=
                  dt * cell_volume * rhs_aw_timelike_from_pi0;
              diag_rhs_aw_mode0_total_step += dt * cell_volume * rhs_aw;
            }
            if (n == 2) {
              diag_rhs_aw_mode2_even_bridge_step +=
                  dt * cell_volume * src_aw_mode2_even_bridge;
              diag_rhs_aw_mode2_even_bridge_abs_step +=
                  dt * cell_volume * std::abs(src_aw_mode2_even_bridge);
              const Real ew_mode2 = -pi_old(idx_aw, k, j, i) - tables.ApplyID3Raising(a0_modes, n);
              const Real bridge_power_mode2 = -src_aw_mode2_even_bridge * ew_mode2;
              const Real bridge_power_sum = bridge_power_mode0_cell + bridge_power_mode2;
              diag_bridge_power_mode2_step += dt * cell_volume * bridge_power_mode2;
              diag_bridge_power_mode2_abs_step +=
                  dt * cell_volume * std::abs(bridge_power_mode2);
              diag_bridge_power_sum_step += dt * cell_volume * bridge_power_sum;
              diag_bridge_power_sum_abs_step +=
                  dt * cell_volume * std::abs(bridge_power_sum);
            }
            if (n == 1) {
              diag_rhs_aw_mode1_w0_response_step +=
                  dt * cell_volume * src_aw_w0_response;
              diag_rhs_aw_mode1_w0_response_abs_step +=
                  dt * cell_volume * std::abs(src_aw_w0_response);
              const Real ew_mode1 =
                  -pi_old(idx_aw, k, j, i) - tables.ApplyID3Raising(a0_modes, n);
              const Real response_bridge_power_mode1 =
                  -src_aw_w0_response * ew_mode1;
              const Real response_bridge_power_sum =
                  response_bridge_power_mode0_cell + response_bridge_power_mode1;
              diag_response_bridge_power_mode1_step +=
                  dt * cell_volume * response_bridge_power_mode1;
              diag_response_bridge_power_mode1_abs_step +=
                  dt * cell_volume * std::abs(response_bridge_power_mode1);
              diag_response_bridge_power_sum_step +=
                  dt * cell_volume * response_bridge_power_sum;
              diag_response_bridge_power_sum_abs_step +=
                  dt * cell_volume * std::abs(response_bridge_power_sum);
            }

            if ((response_w0_enable || response_lambda_enable) && (n >= 1)) {
              // Use total higher-mode EM reservoir energy (not just n=1 odd channels)
              // so either response path (w0 or lambda) can be driven in even-seeded runs.
              const Real reservoir_proxy_density =
                  0.5 * ((a_old(idx_a0, k, j, i) * a_old(idx_a0, k, j, i)) +
                         (a_old(idx_ax, k, j, i) * a_old(idx_ax, k, j, i)) +
                         (a_old(idx_ay, k, j, i) * a_old(idx_ay, k, j, i)) +
                         (a_old(idx_az, k, j, i) * a_old(idx_az, k, j, i)) +
                         (a_old(idx_aw, k, j, i) * a_old(idx_aw, k, j, i)) +
                         (pi_old(idx_a0, k, j, i) * pi_old(idx_a0, k, j, i)) +
                         (pi_old(idx_ax, k, j, i) * pi_old(idx_ax, k, j, i)) +
                         (pi_old(idx_ay, k, j, i) * pi_old(idx_ay, k, j, i)) +
                         (pi_old(idx_az, k, j, i) * pi_old(idx_az, k, j, i)) +
                         (pi_old(idx_aw, k, j, i) * pi_old(idx_aw, k, j, i)));
              response_drive_energy_integral_step += cell_volume * reservoir_proxy_density;
              if ((n % 2) == 0) {
                response_drive_even_energy_integral_step +=
                    cell_volume * reservoir_proxy_density;
              } else {
                response_drive_odd_energy_integral_step +=
                    cell_volume * reservoir_proxy_density;
              }
              response_drive_volume_step += cell_volume;
            }

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
          const Real pi0_mode0_new = pi_new(pi0_mode0_idx, k, j, i);
          const Real pix_mode0_new = pi_new(pix_mode0_idx, k, j, i);
          const Real piy_mode0_new = pi_new(piy_mode0_idx, k, j, i);
          const Real piz_mode0_new = pi_new(piz_mode0_idx, k, j, i);
          const Real piw_mode0_new = pi_new(piw_mode0_idx, k, j, i);
          diag_srcpi0_mode0_step += cell_volume * (pi0_mode0_new - pi0_mode0_old);
          diag_srcpix_mode0_step += cell_volume * (pix_mode0_new - pix_mode0_old);
          diag_srcpiy_mode0_step += cell_volume * (piy_mode0_new - piy_mode0_old);
          diag_srcpiz_mode0_step += cell_volume * (piz_mode0_new - piz_mode0_old);
          diag_srcpiw_mode0_step += cell_volume * (piw_mode0_new - piw_mode0_old);

          for (int n = 0; n < n_modes; ++n) {
            a0_modes[n] = a_new(EMIndex(n, kCompA0), k, j, i);
          }
          for (int n = 0; n < n_modes; ++n) {
            ew_modes[n] = -pi_new(EMIndex(n, kCompAW), k, j, i) -
                          tables.ApplyID3Raising(a0_modes, n);
          }

          Real jw_ew_density = 0.0;
          Real ja_ea_density = 0.0;
          for (int q = 0; q < tables.NumQuadrature(); ++q) {
            Real jx_at_node = 0.0;
            Real jy_at_node = 0.0;
            Real jz_at_node = 0.0;
            Real jw_at_node = 0.0;
            Real ew_at_node = 0.0;
            for (int n = 0; n < n_modes; ++n) {
              const Real phi_nq = tables.Phi(n, q);
              jx_at_node += jx_modes[n] * phi_nq;
              jy_at_node += jy_modes[n] * phi_nq;
              jz_at_node += jz_modes[n] * phi_nq;
              jw_at_node += jw_modes[n] * phi_nq;
              ew_at_node += ew_modes[n] * phi_nq;
            }
            const Real ex = ex_nodes[q];
            const Real ey = ey_nodes[q];
            const Real ez = ez_nodes[q];
            ja_ea_density +=
                weights[q] *
                ((jx_at_node * ex) + (jy_at_node * ey) + (jz_at_node * ez));
            jw_ew_density += weights[q] * jw_at_node * ew_at_node;
          }

          const Real s_leak_density =
              (n_modes > 1) ? (-(inv_lambda_root2 * jw_modes[1])) : 0.0;
          diag_ja_ea_step += dt * cell_volume * ja_ea_density;
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

  auto select_bridge_integral = [&](const std::string &channel) -> Real {
    if (channel == "response") {
      return diag_response_bridge_power_sum_step;
    }
    if (channel == "combined") {
      return diag_bridge_power_sum_step + diag_response_bridge_power_sum_step;
    }
    return diag_bridge_power_sum_step;
  };

  const Real response_drive_density =
      (response_drive_volume_step > 0.0)
          ? (response_drive_energy_integral_step / response_drive_volume_step)
          : 0.0;
  const Real response_drive_parity_density =
      (response_drive_volume_step > 0.0)
          ? ((response_drive_odd_energy_integral_step -
              response_drive_even_energy_integral_step) /
             response_drive_volume_step)
          : 0.0;
  const Real response_drive_work_density =
      ((response_drive_volume_step > 0.0) && (dt > 0.0))
          ? (diag_jw_ew_step / (dt * response_drive_volume_step))
          : 0.0;
  const Real response_drive_leak_density =
      ((response_drive_volume_step > 0.0) && (dt > 0.0))
          ? (diag_s_leak_step / (dt * response_drive_volume_step))
          : 0.0;

  if (response_w0_enable) {
    Real response_drive_bridge_integral =
        select_bridge_integral(response_w0_drive_from_bridge_channel);
    if (response_w0_drive_from_bridge_abs) {
      response_drive_bridge_integral = std::abs(response_drive_bridge_integral);
    }
    const Real response_drive_bridge_density =
        ((response_drive_volume_step > 0.0) && (dt > 0.0))
            ? (response_drive_bridge_integral / (dt * response_drive_volume_step))
            : 0.0;
    Real response_drive_work_term = response_drive_work_density;
    if (response_w0_drive_from_work_abs) {
      response_drive_work_term = std::abs(response_drive_work_term);
    }
    Real response_drive_leak_term = response_drive_leak_density;
    if (response_w0_drive_from_leak_abs) {
      response_drive_leak_term = std::abs(response_drive_leak_term);
    }
    response_w0_drive_reservoir =
        response_w0_drive_gain * std::log1p(std::max(response_drive_density, 0.0));
    response_w0_drive_bridge =
        response_w0_drive_from_bridge_gain * response_drive_bridge_density;
    response_w0_drive_parity =
        response_w0_drive_from_parity_gain * response_drive_parity_density;
    response_w0_drive_work =
        response_w0_drive_from_work_gain * response_drive_work_term;
    response_w0_drive_leak =
        response_w0_drive_from_leak_gain * response_drive_leak_term;
    response_w0_drive = response_w0_bias + response_w0_drive_reservoir +
                        response_w0_drive_bridge + response_w0_drive_parity +
                        response_w0_drive_work + response_w0_drive_leak;
    response_w0_force = response_w0_drive - (response_w0_stiffness * response_w0) -
                        (response_w0_damping * response_w0_dot);
    const Real response_w0_dot_old = response_w0_dot;
    const Real response_w0_ddot = response_w0_force / response_w0_mass;
    response_w0_dot += dt * response_w0_ddot;
    const Real response_w0_dot_mid = 0.5 * (response_w0_dot_old + response_w0_dot);
    const Real response_w0_power_drive = response_w0_drive * response_w0_dot_mid;
    const Real response_w0_power_stiffness =
        -(response_w0_stiffness * response_w0) * response_w0_dot_mid;
    const Real response_w0_power_damping =
        -(response_w0_damping * response_w0_dot_mid) * response_w0_dot_mid;
    const Real response_w0_power_net = response_w0_force * response_w0_dot_mid;
    diag_response_w0_power_drive_step = dt * response_w0_power_drive;
    diag_response_w0_power_stiffness_step = dt * response_w0_power_stiffness;
    diag_response_w0_power_damping_step = dt * response_w0_power_damping;
    diag_response_w0_power_net_step = dt * response_w0_power_net;
    response_w0 += dt * response_w0_dot;
    if (response_w0_max_abs > 0.0) {
      response_w0 = std::clamp(response_w0, -response_w0_max_abs, response_w0_max_abs);
    }
    response_w0_energy = 0.5 * ((response_w0_mass * response_w0_dot * response_w0_dot) +
                                (response_w0_stiffness * response_w0 * response_w0));
  } else {
    response_w0 = 0.0;
    response_w0_dot = 0.0;
    response_w0_drive = 0.0;
    response_w0_drive_reservoir = 0.0;
    response_w0_drive_bridge = 0.0;
    response_w0_drive_parity = 0.0;
    response_w0_drive_work = 0.0;
    response_w0_drive_leak = 0.0;
    response_w0_force = 0.0;
    response_w0_energy = 0.0;
    response_w0_initialized = false;
  }

  if (response_lambda_enable) {
    Real response_lambda_bridge_integral =
        select_bridge_integral(response_lambda_drive_from_bridge_channel);
    if (response_lambda_drive_from_bridge_abs) {
      response_lambda_bridge_integral = std::abs(response_lambda_bridge_integral);
    }
    const Real response_lambda_bridge_density =
        ((response_drive_volume_step > 0.0) && (dt > 0.0))
            ? (response_lambda_bridge_integral / (dt * response_drive_volume_step))
            : 0.0;
    Real response_lambda_work_term = response_drive_work_density;
    if (response_lambda_drive_from_work_abs) {
      response_lambda_work_term = std::abs(response_lambda_work_term);
    }
    Real response_lambda_leak_term = response_drive_leak_density;
    if (response_lambda_drive_from_leak_abs) {
      response_lambda_leak_term = std::abs(response_lambda_leak_term);
    }
    response_lambda_drive_reservoir =
        response_lambda_drive_gain * std::log1p(std::max(response_drive_density, 0.0));
    response_lambda_drive_bridge =
        response_lambda_drive_from_bridge_gain * response_lambda_bridge_density;
    response_lambda_drive_work =
        response_lambda_drive_from_work_gain * response_lambda_work_term;
    response_lambda_drive_leak =
        response_lambda_drive_from_leak_gain * response_lambda_leak_term;
    response_lambda_drive =
        response_lambda_bias + response_lambda_drive_reservoir +
        response_lambda_drive_bridge + response_lambda_drive_work +
        response_lambda_drive_leak;
    response_lambda_force =
        response_lambda_drive - (response_lambda_stiffness * response_lambda_fraction) -
        (response_lambda_damping * response_lambda_dot);
    const Real response_lambda_dot_old = response_lambda_dot;
    const Real response_lambda_ddot = response_lambda_force / response_lambda_mass;
    response_lambda_dot += dt * response_lambda_ddot;
    const Real response_lambda_dot_mid =
        0.5 * (response_lambda_dot_old + response_lambda_dot);
    const Real response_lambda_power_drive =
        response_lambda_drive * response_lambda_dot_mid;
    const Real response_lambda_power_stiffness =
        -(response_lambda_stiffness * response_lambda_fraction) * response_lambda_dot_mid;
    const Real response_lambda_power_damping =
        -(response_lambda_damping * response_lambda_dot_mid) * response_lambda_dot_mid;
    const Real response_lambda_power_net = response_lambda_force * response_lambda_dot_mid;
    diag_response_lambda_power_drive_step = dt * response_lambda_power_drive;
    diag_response_lambda_power_stiffness_step = dt * response_lambda_power_stiffness;
    diag_response_lambda_power_damping_step = dt * response_lambda_power_damping;
    diag_response_lambda_power_net_step = dt * response_lambda_power_net;
    response_lambda_fraction += dt * response_lambda_dot;
    if (response_lambda_max_abs_frac > 0.0) {
      response_lambda_fraction = std::clamp(response_lambda_fraction,
                                            -response_lambda_max_abs_frac,
                                            response_lambda_max_abs_frac);
    }
    response_lambda_energy =
        0.5 * ((response_lambda_mass * response_lambda_dot * response_lambda_dot) +
               (response_lambda_stiffness * response_lambda_fraction *
                response_lambda_fraction));
  } else {
    response_lambda_fraction = 0.0;
    response_lambda_dot = 0.0;
    response_lambda_drive = 0.0;
    response_lambda_drive_reservoir = 0.0;
    response_lambda_drive_bridge = 0.0;
    response_lambda_drive_work = 0.0;
    response_lambda_drive_leak = 0.0;
    response_lambda_force = 0.0;
    response_lambda_energy = 0.0;
    response_lambda_initialized = false;
  }

  const Real response_w0_effective_out =
      response_w0 + (response_w0_velocity_bridge_gain * response_w0_dot);
  Real response_w0_lambda_fraction_from_w0_out = 0.0;
  if (response_w0_enable && response_w0_lambda_shift_enable &&
      (response_w0_lambda_shift_gain != 0.0)) {
    response_w0_lambda_fraction_from_w0_out =
        response_w0_lambda_shift_gain * response_w0_effective_out;
    if (response_w0_lambda_shift_max_frac > 0.0) {
      response_w0_lambda_fraction_from_w0_out =
          std::clamp(response_w0_lambda_fraction_from_w0_out,
                     -response_w0_lambda_shift_max_frac,
                     response_w0_lambda_shift_max_frac);
    }
  }
  const Real response_lambda_fraction_total_out =
      response_w0_lambda_fraction_from_w0_out + response_lambda_fraction;
  Real response_w0_lambda_out = lambda * (1.0 + response_lambda_fraction_total_out);
  response_w0_lambda_out = std::max(response_w0_lambda_out, 1.0e-6 * lambda);

  if (response_w0_projection_enable) {
    response_w0_projection_center = response_w0_projection_gain * response_w0;
    if (response_w0_projection_max_abs > 0.0) {
      response_w0_projection_center =
          std::clamp(response_w0_projection_center, -response_w0_projection_max_abs,
                     response_w0_projection_max_abs);
    }
  } else {
    response_w0_projection_center = 0.0;
  }

  auto *diag_jw_ew = modes_pkg->MutableParam<double>("diag/int_jw_ew");
  auto *diag_ja_ea = modes_pkg->MutableParam<double>("diag/int_ja_ea");
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
  auto *diag_srcmomx_mode0 = modes_pkg->MutableParam<double>("diag/int_srcmomx_mode0");
  auto *diag_srcmomy_mode0 = modes_pkg->MutableParam<double>("diag/int_srcmomy_mode0");
  auto *diag_srcmomz_mode0 = modes_pkg->MutableParam<double>("diag/int_srcmomz_mode0");
  auto *diag_srcmomw_mode0 = modes_pkg->MutableParam<double>("diag/int_srcmomw_mode0");
  auto *diag_srcenergy_mode0 =
      modes_pkg->MutableParam<double>("diag/int_srcenergy_mode0");
  auto *diag_srcpi0_mode0 = modes_pkg->MutableParam<double>("diag/int_srcpi0_mode0");
  auto *diag_srcpix_mode0 = modes_pkg->MutableParam<double>("diag/int_srcpix_mode0");
  auto *diag_srcpiy_mode0 = modes_pkg->MutableParam<double>("diag/int_srcpiy_mode0");
  auto *diag_srcpiz_mode0 = modes_pkg->MutableParam<double>("diag/int_srcpiz_mode0");
  auto *diag_srcpiw_mode0 = modes_pkg->MutableParam<double>("diag/int_srcpiw_mode0");
  auto *diag_src_timelike_a0_from_piw =
      modes_pkg->MutableParam<double>("diag/int_src_timelike_a0_from_piw");
  auto *diag_src_timelike_a0_from_piw_abs =
      modes_pkg->MutableParam<double>("diag/int_src_timelike_a0_from_piw_abs");
  auto *diag_src_timelike_aw_from_pi0 =
      modes_pkg->MutableParam<double>("diag/int_src_timelike_aw_from_pi0");
  auto *diag_src_timelike_aw_from_pi0_abs =
      modes_pkg->MutableParam<double>("diag/int_src_timelike_aw_from_pi0_abs");
  auto *diag_src_em_laplacian_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_laplacian_abs");
  auto *diag_src_em_mass_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_mass_abs");
  auto *diag_src_em_current_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_current_abs");
  auto *diag_src_em_damping_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_damping_abs");
  auto *diag_src_em_geometry_shift_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_geometry_shift_abs");
  auto *diag_src_em_spatial_mixed_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_spatial_mixed_abs");
  auto *diag_src_em_timelike_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_timelike_abs");
  auto *diag_src_em_gauge = modes_pkg->MutableParam<double>("diag/int_src_em_gauge");
  auto *diag_src_em_gauge_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_gauge_abs");
  auto *diag_src_em_gauge_mode0 =
      modes_pkg->MutableParam<double>("diag/int_src_em_gauge_mode0");
  auto *diag_src_em_gauge_mode0_abs =
      modes_pkg->MutableParam<double>("diag/int_src_em_gauge_mode0_abs");
  auto *diag_rhs_a0_mode0_lap =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_lap");
  auto *diag_rhs_a0_mode0_mass =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_mass");
  auto *diag_rhs_a0_mode0_current =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_current");
  auto *diag_rhs_a0_mode0_damping =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_damping");
  auto *diag_rhs_a0_mode0_timelike =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_timelike");
  auto *diag_rhs_a0_mode0_gauge =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_gauge");
  auto *diag_rhs_a0_mode0_total =
      modes_pkg->MutableParam<double>("diag/int_rhs_a0_mode0_total");
  auto *diag_rhs_ay_mode0_lap =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_lap");
  auto *diag_rhs_ay_mode0_mass =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_mass");
  auto *diag_rhs_ay_mode0_current =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_current");
  auto *diag_rhs_ay_mode0_damping =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_damping");
  auto *diag_rhs_ay_mode0_spatial_mixed =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_spatial_mixed");
  auto *diag_rhs_ay_mode0_even_bridge =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_even_bridge");
  auto *diag_rhs_ay_mode0_even_bridge_abs =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_even_bridge_abs");
  auto *diag_rhs_ay_mode0_geometry_shift =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_geometry_shift");
  auto *diag_rhs_ay_mode0_geometry_shift_abs =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_geometry_shift_abs");
  auto *diag_bridge_power_mode0 =
      modes_pkg->MutableParam<double>("diag/int_bridge_power_mode0");
  auto *diag_bridge_power_mode0_abs =
      modes_pkg->MutableParam<double>("diag/int_bridge_power_mode0_abs");
  auto *diag_geometry_shift_power_mode0 =
      modes_pkg->MutableParam<double>("diag/int_geometry_shift_power_mode0");
  auto *diag_geometry_shift_power_mode0_abs =
      modes_pkg->MutableParam<double>("diag/int_geometry_shift_power_mode0_abs");
  auto *diag_rhs_aw_mode2_even_bridge =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode2_even_bridge");
  auto *diag_rhs_aw_mode2_even_bridge_abs =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode2_even_bridge_abs");
  auto *diag_bridge_power_mode2 =
      modes_pkg->MutableParam<double>("diag/int_bridge_power_mode2");
  auto *diag_bridge_power_mode2_abs =
      modes_pkg->MutableParam<double>("diag/int_bridge_power_mode2_abs");
  auto *diag_bridge_power_sum =
      modes_pkg->MutableParam<double>("diag/int_bridge_power_sum");
  auto *diag_bridge_power_sum_abs =
      modes_pkg->MutableParam<double>("diag/int_bridge_power_sum_abs");
  auto *diag_response_w0_initialized =
      modes_pkg->MutableParam<bool>("diag/response_w0_initialized");
  auto *diag_response_w0 = modes_pkg->MutableParam<double>("diag/response_w0");
  auto *diag_response_w0_dot = modes_pkg->MutableParam<double>("diag/response_w0_dot");
  auto *diag_response_w0_drive = modes_pkg->MutableParam<double>("diag/response_w0_drive");
  auto *diag_response_w0_drive_reservoir =
      modes_pkg->MutableParam<double>("diag/response_w0_drive_reservoir");
  auto *diag_response_w0_drive_bridge =
      modes_pkg->MutableParam<double>("diag/response_w0_drive_bridge");
  auto *diag_response_w0_drive_parity =
      modes_pkg->MutableParam<double>("diag/response_w0_drive_parity");
  auto *diag_response_w0_drive_work =
      modes_pkg->MutableParam<double>("diag/response_w0_drive_work");
  auto *diag_response_w0_drive_leak =
      modes_pkg->MutableParam<double>("diag/response_w0_drive_leak");
  auto *diag_response_w0_force = modes_pkg->MutableParam<double>("diag/response_w0_force");
  auto *diag_response_w0_energy = modes_pkg->MutableParam<double>("diag/response_w0_energy");
  auto *diag_response_w0_geometry_center =
      modes_pkg->MutableParam<double>("diag/response_w0_geometry_center");
  auto *diag_response_w0_lambda_fraction =
      modes_pkg->MutableParam<double>("diag/response_w0_lambda_fraction");
  auto *diag_response_w0_lambda_eff =
      modes_pkg->MutableParam<double>("diag/response_w0_lambda_eff");
  auto *diag_response_lambda_initialized =
      modes_pkg->MutableParam<bool>("diag/response_lambda_initialized");
  auto *diag_response_lambda_fraction =
      modes_pkg->MutableParam<double>("diag/response_lambda_fraction");
  auto *diag_response_lambda_dot =
      modes_pkg->MutableParam<double>("diag/response_lambda_dot");
  auto *diag_response_lambda_drive =
      modes_pkg->MutableParam<double>("diag/response_lambda_drive");
  auto *diag_response_lambda_drive_reservoir =
      modes_pkg->MutableParam<double>("diag/response_lambda_drive_reservoir");
  auto *diag_response_lambda_drive_bridge =
      modes_pkg->MutableParam<double>("diag/response_lambda_drive_bridge");
  auto *diag_response_lambda_drive_work =
      modes_pkg->MutableParam<double>("diag/response_lambda_drive_work");
  auto *diag_response_lambda_drive_leak =
      modes_pkg->MutableParam<double>("diag/response_lambda_drive_leak");
  auto *diag_response_lambda_force =
      modes_pkg->MutableParam<double>("diag/response_lambda_force");
  auto *diag_response_lambda_energy =
      modes_pkg->MutableParam<double>("diag/response_lambda_energy");
  auto *diag_projection_center_w =
      modes_pkg->MutableParam<double>("diag/projection_center_w");
  auto *diag_response_w0_power_drive =
      modes_pkg->MutableParam<double>("diag/int_response_w0_power_drive");
  auto *diag_response_w0_power_stiffness =
      modes_pkg->MutableParam<double>("diag/int_response_w0_power_stiffness");
  auto *diag_response_w0_power_damping =
      modes_pkg->MutableParam<double>("diag/int_response_w0_power_damping");
  auto *diag_response_w0_power_net =
      modes_pkg->MutableParam<double>("diag/int_response_w0_power_net");
  auto *diag_response_lambda_power_drive =
      modes_pkg->MutableParam<double>("diag/int_response_lambda_power_drive");
  auto *diag_response_lambda_power_stiffness =
      modes_pkg->MutableParam<double>("diag/int_response_lambda_power_stiffness");
  auto *diag_response_lambda_power_damping =
      modes_pkg->MutableParam<double>("diag/int_response_lambda_power_damping");
  auto *diag_response_lambda_power_net =
      modes_pkg->MutableParam<double>("diag/int_response_lambda_power_net");
  auto *diag_rhs_ay_mode0_w0_response =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_w0_response");
  auto *diag_rhs_ay_mode0_w0_response_abs =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_w0_response_abs");
  auto *diag_rhs_aw_mode1_w0_response =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode1_w0_response");
  auto *diag_rhs_aw_mode1_w0_response_abs =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode1_w0_response_abs");
  auto *diag_response_bridge_power_mode0 =
      modes_pkg->MutableParam<double>("diag/int_response_bridge_power_mode0");
  auto *diag_response_bridge_power_mode0_abs =
      modes_pkg->MutableParam<double>("diag/int_response_bridge_power_mode0_abs");
  auto *diag_response_bridge_power_mode1 =
      modes_pkg->MutableParam<double>("diag/int_response_bridge_power_mode1");
  auto *diag_response_bridge_power_mode1_abs =
      modes_pkg->MutableParam<double>("diag/int_response_bridge_power_mode1_abs");
  auto *diag_response_bridge_power_sum =
      modes_pkg->MutableParam<double>("diag/int_response_bridge_power_sum");
  auto *diag_response_bridge_power_sum_abs =
      modes_pkg->MutableParam<double>("diag/int_response_bridge_power_sum_abs");
  auto *diag_rhs_ay_mode0_total =
      modes_pkg->MutableParam<double>("diag/int_rhs_ay_mode0_total");
  auto *diag_rhs_aw_mode0_lap =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode0_lap");
  auto *diag_rhs_aw_mode0_current =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode0_current");
  auto *diag_rhs_aw_mode0_damping =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode0_damping");
  auto *diag_rhs_aw_mode0_spatial_mixed =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode0_spatial_mixed");
  auto *diag_rhs_aw_mode0_timelike =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode0_timelike");
  auto *diag_rhs_aw_mode0_total =
      modes_pkg->MutableParam<double>("diag/int_rhs_aw_mode0_total");
  auto *diag_gauge_l1 = modes_pkg->MutableParam<double>("diag/gauge_l1");
  auto *diag_gauge_l2 = modes_pkg->MutableParam<double>("diag/gauge_l2");
  auto *diag_gauge_max_abs =
      modes_pkg->MutableParam<double>("diag/gauge_max_abs");
  auto *diag_gauge_mode0_l1 = modes_pkg->MutableParam<double>("diag/gauge_mode0_l1");
  auto *diag_gauge_mode0_l2 = modes_pkg->MutableParam<double>("diag/gauge_mode0_l2");
  auto *diag_gauge_mode0_max_abs =
      modes_pkg->MutableParam<double>("diag/gauge_mode0_max_abs");
  *diag_ja_ea += diag_ja_ea_step;
  *diag_jw_ew += diag_jw_ew_step;
  *diag_s_leak += diag_s_leak_step;
  *diag_s_leak_abs += diag_s_leak_abs_step;
  *diag_cont_local_l1 = diag_cont_local_l1_step;
  *diag_cont_local_l2 = diag_cont_local_l2_step;
  *diag_cont_local_max_abs = diag_cont_local_max_abs_step;
  *diag_cont_mode0_l1 = diag_cont_mode0_l1_step;
  *diag_cont_mode0_l2 = diag_cont_mode0_l2_step;
  *diag_cont_mode0_max_abs = diag_cont_mode0_max_abs_step;
  *diag_srcmomx_mode0 += diag_srcmomx_mode0_step;
  *diag_srcmomy_mode0 += diag_srcmomy_mode0_step;
  *diag_srcmomz_mode0 += diag_srcmomz_mode0_step;
  *diag_srcmomw_mode0 += diag_srcmomw_mode0_step;
  *diag_srcenergy_mode0 += diag_srcenergy_mode0_step;
  *diag_srcpi0_mode0 += diag_srcpi0_mode0_step;
  *diag_srcpix_mode0 += diag_srcpix_mode0_step;
  *diag_srcpiy_mode0 += diag_srcpiy_mode0_step;
  *diag_srcpiz_mode0 += diag_srcpiz_mode0_step;
  *diag_srcpiw_mode0 += diag_srcpiw_mode0_step;
  *diag_src_timelike_a0_from_piw += diag_src_timelike_a0_from_piw_step;
  *diag_src_timelike_a0_from_piw_abs += diag_src_timelike_a0_from_piw_abs_step;
  *diag_src_timelike_aw_from_pi0 += diag_src_timelike_aw_from_pi0_step;
  *diag_src_timelike_aw_from_pi0_abs += diag_src_timelike_aw_from_pi0_abs_step;
  *diag_src_em_laplacian_abs += diag_src_em_laplacian_abs_step;
  *diag_src_em_mass_abs += diag_src_em_mass_abs_step;
  *diag_src_em_current_abs += diag_src_em_current_abs_step;
  *diag_src_em_damping_abs += diag_src_em_damping_abs_step;
  *diag_src_em_geometry_shift_abs += diag_src_em_geometry_shift_abs_step;
  *diag_src_em_spatial_mixed_abs += diag_src_em_spatial_mixed_abs_step;
  *diag_src_em_timelike_abs += diag_src_em_timelike_abs_step;
  *diag_src_em_gauge += diag_src_em_gauge_step;
  *diag_src_em_gauge_abs += diag_src_em_gauge_abs_step;
  *diag_src_em_gauge_mode0 += diag_src_em_gauge_mode0_step;
  *diag_src_em_gauge_mode0_abs += diag_src_em_gauge_mode0_abs_step;
  *diag_rhs_a0_mode0_lap += diag_rhs_a0_mode0_lap_step;
  *diag_rhs_a0_mode0_mass += diag_rhs_a0_mode0_mass_step;
  *diag_rhs_a0_mode0_current += diag_rhs_a0_mode0_current_step;
  *diag_rhs_a0_mode0_damping += diag_rhs_a0_mode0_damping_step;
  *diag_rhs_a0_mode0_timelike += diag_rhs_a0_mode0_timelike_step;
  *diag_rhs_a0_mode0_gauge += diag_rhs_a0_mode0_gauge_step;
  *diag_rhs_a0_mode0_total += diag_rhs_a0_mode0_total_step;
  *diag_rhs_ay_mode0_lap += diag_rhs_ay_mode0_lap_step;
  *diag_rhs_ay_mode0_mass += diag_rhs_ay_mode0_mass_step;
  *diag_rhs_ay_mode0_current += diag_rhs_ay_mode0_current_step;
  *diag_rhs_ay_mode0_damping += diag_rhs_ay_mode0_damping_step;
  *diag_rhs_ay_mode0_spatial_mixed += diag_rhs_ay_mode0_spatial_mixed_step;
  *diag_rhs_ay_mode0_even_bridge += diag_rhs_ay_mode0_even_bridge_step;
  *diag_rhs_ay_mode0_even_bridge_abs += diag_rhs_ay_mode0_even_bridge_abs_step;
  *diag_rhs_ay_mode0_geometry_shift += diag_rhs_ay_mode0_geometry_shift_step;
  *diag_rhs_ay_mode0_geometry_shift_abs +=
      diag_rhs_ay_mode0_geometry_shift_abs_step;
  *diag_bridge_power_mode0 += diag_bridge_power_mode0_step;
  *diag_bridge_power_mode0_abs += diag_bridge_power_mode0_abs_step;
  *diag_geometry_shift_power_mode0 += diag_geometry_shift_power_mode0_step;
  *diag_geometry_shift_power_mode0_abs +=
      diag_geometry_shift_power_mode0_abs_step;
  *diag_rhs_aw_mode2_even_bridge += diag_rhs_aw_mode2_even_bridge_step;
  *diag_rhs_aw_mode2_even_bridge_abs += diag_rhs_aw_mode2_even_bridge_abs_step;
  *diag_bridge_power_mode2 += diag_bridge_power_mode2_step;
  *diag_bridge_power_mode2_abs += diag_bridge_power_mode2_abs_step;
  *diag_bridge_power_sum += diag_bridge_power_sum_step;
  *diag_bridge_power_sum_abs += diag_bridge_power_sum_abs_step;
  *diag_response_w0_initialized = response_w0_initialized;
  *diag_response_w0 = response_w0;
  *diag_response_w0_dot = response_w0_dot;
  *diag_response_w0_drive = response_w0_drive;
  *diag_response_w0_drive_reservoir = response_w0_drive_reservoir;
  *diag_response_w0_drive_bridge = response_w0_drive_bridge;
  *diag_response_w0_drive_parity = response_w0_drive_parity;
  *diag_response_w0_drive_work = response_w0_drive_work;
  *diag_response_w0_drive_leak = response_w0_drive_leak;
  *diag_response_w0_force = response_w0_force;
  *diag_response_w0_energy = response_w0_energy;
  *diag_response_w0_geometry_center = response_w0;
  *diag_response_w0_lambda_fraction = response_lambda_fraction_total_out;
  *diag_response_w0_lambda_eff = response_w0_lambda_out;
  *diag_response_lambda_initialized = response_lambda_initialized;
  *diag_response_lambda_fraction = response_lambda_fraction;
  *diag_response_lambda_dot = response_lambda_dot;
  *diag_response_lambda_drive = response_lambda_drive;
  *diag_response_lambda_drive_reservoir = response_lambda_drive_reservoir;
  *diag_response_lambda_drive_bridge = response_lambda_drive_bridge;
  *diag_response_lambda_drive_work = response_lambda_drive_work;
  *diag_response_lambda_drive_leak = response_lambda_drive_leak;
  *diag_response_lambda_force = response_lambda_force;
  *diag_response_lambda_energy = response_lambda_energy;
  *diag_projection_center_w = response_w0_projection_center;
  *diag_response_w0_power_drive += diag_response_w0_power_drive_step;
  *diag_response_w0_power_stiffness += diag_response_w0_power_stiffness_step;
  *diag_response_w0_power_damping += diag_response_w0_power_damping_step;
  *diag_response_w0_power_net += diag_response_w0_power_net_step;
  *diag_response_lambda_power_drive += diag_response_lambda_power_drive_step;
  *diag_response_lambda_power_stiffness += diag_response_lambda_power_stiffness_step;
  *diag_response_lambda_power_damping += diag_response_lambda_power_damping_step;
  *diag_response_lambda_power_net += diag_response_lambda_power_net_step;
  *diag_rhs_ay_mode0_w0_response += diag_rhs_ay_mode0_w0_response_step;
  *diag_rhs_ay_mode0_w0_response_abs += diag_rhs_ay_mode0_w0_response_abs_step;
  *diag_rhs_aw_mode1_w0_response += diag_rhs_aw_mode1_w0_response_step;
  *diag_rhs_aw_mode1_w0_response_abs += diag_rhs_aw_mode1_w0_response_abs_step;
  *diag_response_bridge_power_mode0 += diag_response_bridge_power_mode0_step;
  *diag_response_bridge_power_mode0_abs += diag_response_bridge_power_mode0_abs_step;
  *diag_response_bridge_power_mode1 += diag_response_bridge_power_mode1_step;
  *diag_response_bridge_power_mode1_abs += diag_response_bridge_power_mode1_abs_step;
  *diag_response_bridge_power_sum += diag_response_bridge_power_sum_step;
  *diag_response_bridge_power_sum_abs += diag_response_bridge_power_sum_abs_step;
  *diag_rhs_ay_mode0_total += diag_rhs_ay_mode0_total_step;
  *diag_rhs_aw_mode0_lap += diag_rhs_aw_mode0_lap_step;
  *diag_rhs_aw_mode0_current += diag_rhs_aw_mode0_current_step;
  *diag_rhs_aw_mode0_damping += diag_rhs_aw_mode0_damping_step;
  *diag_rhs_aw_mode0_spatial_mixed += diag_rhs_aw_mode0_spatial_mixed_step;
  *diag_rhs_aw_mode0_timelike += diag_rhs_aw_mode0_timelike_step;
  *diag_rhs_aw_mode0_total += diag_rhs_aw_mode0_total_step;
  *diag_gauge_l1 = diag_gauge_l1_step;
  *diag_gauge_l2 = diag_gauge_l2_step;
  *diag_gauge_max_abs = diag_gauge_max_abs_step;
  *diag_gauge_mode0_l1 = diag_gauge_mode0_l1_step;
  *diag_gauge_mode0_l2 = diag_gauge_mode0_l2_step;
  *diag_gauge_mode0_max_abs = diag_gauge_mode0_max_abs_step;
}

} // namespace Modes4D
