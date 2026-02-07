#include "package4d_modes.hpp"

#include <algorithm>

#include "diagnostics4d.hpp"
#include "em4d_modes.hpp"
#include "mode_tables.hpp"
#include "plasma4d_modes.hpp"
#include "utils/error_checking.hpp"

namespace Modes4D {

std::shared_ptr<parthenon::StateDescriptor> Initialize(parthenon::ParameterInput *pin) {
  auto pkg = std::make_shared<parthenon::StateDescriptor>("modes4d");

  const bool enabled = pin->GetOrAddBoolean("modes4d", "enabled", false);
  pkg->AddParam<bool>("enabled", enabled);
  if (!enabled) {
    return pkg;
  }

  const int n_modes = pin->GetOrAddInteger("modes4d", "n_modes", 4);
  const int n_quadrature =
      pin->GetOrAddInteger("modes4d", "n_quadrature", std::max(n_modes, n_modes + 2));
  const double lambda = pin->GetOrAddReal("modes4d", "lambda", 1.0);
  const double em_c_wave = pin->GetOrAddReal("modes4d", "em_c_wave", 1.0);
  const double em_damping = pin->GetOrAddReal("modes4d", "em_damping", 0.0);
  const double em_mu0 = pin->GetOrAddReal("modes4d", "em_mu0", 1.0);
  const bool em_conservative_transport =
      pin->GetOrAddBoolean("modes4d", "em_conservative_transport", false);
  const double plasma_qom_ion = pin->GetOrAddReal("modes4d", "plasma_qom_ion", 1.0);
  const double plasma_qom_electron =
      pin->GetOrAddReal("modes4d", "plasma_qom_electron", -1.0);
  const double hydro_gamma = pin->GetOrAddReal("hydro", "gamma", 5.0 / 3.0);
  const double plasma_gamma = pin->GetOrAddReal("modes4d", "plasma_gamma", hydro_gamma);
  const double plasma_force_source_gain =
      pin->GetOrAddReal("modes4d", "plasma_force_source_gain", 1.0);
  const double plasma_momw_source_gain =
      pin->GetOrAddReal("modes4d", "plasma_momw_source_gain", 1.0);
  const double plasma_momw_pressure_source_gain =
      pin->GetOrAddReal("modes4d", "plasma_momw_pressure_source_gain", 0.0);
  const double plasma_momw_damping =
      pin->GetOrAddReal("modes4d", "plasma_momw_damping", 0.0);
  const double plasma_w_flux_source_gain =
      pin->GetOrAddReal("modes4d", "plasma_w_flux_source_gain", 0.0);
  const double plasma_rho_divj_gain =
      pin->GetOrAddReal("modes4d", "plasma_rho_divj_gain", 0.0);
  const double plasma_rho_floor = pin->GetOrAddReal("modes4d", "plasma_rho_floor", 1.0e-12);
  const double plasma_energy_source_gain =
      pin->GetOrAddReal("modes4d", "plasma_energy_source_gain", 1.0);
  const double plasma_energy_floor =
      pin->GetOrAddReal("modes4d", "plasma_energy_floor", 0.0);
  const double plasma_pressure_transport_gain =
      pin->GetOrAddReal("modes4d", "plasma_pressure_transport_gain", 0.0);
  const double plasma_pressure_rusanov_gain =
      pin->GetOrAddReal("modes4d", "plasma_pressure_rusanov_gain", 1.0);
  const double plasma_pressure_signal_speed_cap =
      pin->GetOrAddReal("modes4d", "plasma_pressure_signal_speed_cap", 10.0);
  const double plasma_pressure_flux_relative_cap =
      pin->GetOrAddReal("modes4d", "plasma_pressure_flux_relative_cap", 50.0);
  const int plasma_pressure_transport_max_mode =
      pin->GetOrAddInteger("modes4d", "plasma_pressure_transport_max_mode", 0);
  const bool plasma_pseudospectral_transport =
      pin->GetOrAddBoolean("modes4d", "plasma_pseudospectral_transport", false);
  const double plasma_pressure_floor =
      pin->GetOrAddReal("modes4d", "plasma_pressure_floor", 0.0);
  const double pulse_amp = pin->GetOrAddReal("problem/em4d_pulse", "amplitude", 1.0e-3);
  const double pulse_sigma = pin->GetOrAddReal("problem/em4d_pulse", "sigma", 0.08);
  const double pulse_x0 = pin->GetOrAddReal("problem/em4d_pulse", "x0", 0.0);

  PARTHENON_REQUIRE(n_modes > 0, "modes4d/n_modes must be > 0");
  PARTHENON_REQUIRE(n_quadrature >= n_modes,
                    "modes4d/n_quadrature must be >= modes4d/n_modes");
  PARTHENON_REQUIRE(lambda > 0.0, "modes4d/lambda must be > 0");
  PARTHENON_REQUIRE(plasma_gamma > 1.0, "modes4d/plasma_gamma must be > 1");
  PARTHENON_REQUIRE(plasma_pressure_rusanov_gain >= 0.0,
                    "modes4d/plasma_pressure_rusanov_gain must be >= 0");
  PARTHENON_REQUIRE(plasma_pressure_signal_speed_cap > 0.0,
                    "modes4d/plasma_pressure_signal_speed_cap must be > 0");
  PARTHENON_REQUIRE(plasma_pressure_flux_relative_cap > 0.0,
                    "modes4d/plasma_pressure_flux_relative_cap must be > 0");
  PARTHENON_REQUIRE(plasma_pressure_transport_max_mode >= 0,
                    "modes4d/plasma_pressure_transport_max_mode must be >= 0");

  ModeConfig config;
  config.n_modes = n_modes;
  config.n_quadrature = n_quadrature;
  config.lambda = lambda;

  pkg->AddParam<int>("n_modes", n_modes);
  pkg->AddParam<int>("n_quadrature", n_quadrature);
  pkg->AddParam<double>("lambda", lambda);
  pkg->AddParam<double>("em4d/c_wave", em_c_wave);
  pkg->AddParam<double>("em4d/damping", em_damping);
  pkg->AddParam<double>("em4d/mu0", em_mu0);
  pkg->AddParam<bool>("em4d/use_conservative_transport", em_conservative_transport);
  pkg->AddParam<double>("plasma4d/qom_ion", plasma_qom_ion);
  pkg->AddParam<double>("plasma4d/qom_electron", plasma_qom_electron);
  pkg->AddParam<double>("plasma4d/gamma", plasma_gamma);
  pkg->AddParam<double>("plasma4d/force_source_gain", plasma_force_source_gain);
  pkg->AddParam<double>("plasma4d/momw_source_gain", plasma_momw_source_gain);
  pkg->AddParam<double>("plasma4d/momw_pressure_source_gain",
                        plasma_momw_pressure_source_gain);
  pkg->AddParam<double>("plasma4d/momw_damping", plasma_momw_damping);
  pkg->AddParam<double>("plasma4d/w_flux_source_gain", plasma_w_flux_source_gain);
  pkg->AddParam<double>("plasma4d/rho_divj_gain", plasma_rho_divj_gain);
  pkg->AddParam<double>("plasma4d/rho_floor", plasma_rho_floor);
  pkg->AddParam<double>("plasma4d/energy_source_gain", plasma_energy_source_gain);
  pkg->AddParam<double>("plasma4d/energy_floor", plasma_energy_floor);
  pkg->AddParam<double>("plasma4d/pressure_transport_gain", plasma_pressure_transport_gain);
  pkg->AddParam<double>("plasma4d/pressure_rusanov_gain", plasma_pressure_rusanov_gain);
  pkg->AddParam<double>("plasma4d/pressure_signal_speed_cap",
                        plasma_pressure_signal_speed_cap);
  pkg->AddParam<double>("plasma4d/pressure_flux_relative_cap",
                        plasma_pressure_flux_relative_cap);
  pkg->AddParam<int>("plasma4d/pressure_transport_max_mode",
                     plasma_pressure_transport_max_mode);
  pkg->AddParam<bool>("plasma4d/use_pseudospectral_transport",
                      plasma_pseudospectral_transport);
  pkg->AddParam<double>("plasma4d/pressure_floor", plasma_pressure_floor);
  pkg->AddParam<double>("em4d_pulse/amplitude", pulse_amp);
  pkg->AddParam<double>("em4d_pulse/sigma", pulse_sigma);
  pkg->AddParam<double>("em4d_pulse/x0", pulse_x0);
  pkg->AddParam<ModeTables>("mode_tables", ModeTables(config));

  RegisterEMVariables(pkg.get(), n_modes, em_conservative_transport);
  RegisterPlasmaVariables(pkg.get(), n_modes);
  RegisterDiagnostics(pkg.get());

  return pkg;
}

} // namespace Modes4D
