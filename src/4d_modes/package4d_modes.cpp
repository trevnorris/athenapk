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

  PARTHENON_REQUIRE(n_modes > 0, "modes4d/n_modes must be > 0");
  PARTHENON_REQUIRE(n_quadrature >= n_modes,
                    "modes4d/n_quadrature must be >= modes4d/n_modes");
  PARTHENON_REQUIRE(lambda > 0.0, "modes4d/lambda must be > 0");

  ModeConfig config;
  config.n_modes = n_modes;
  config.n_quadrature = n_quadrature;
  config.lambda = lambda;

  pkg->AddParam<int>("n_modes", n_modes);
  pkg->AddParam<int>("n_quadrature", n_quadrature);
  pkg->AddParam<double>("lambda", lambda);
  pkg->AddParam<ModeTables>("mode_tables", ModeTables(config));

  RegisterEMVariables(pkg.get(), n_modes);
  RegisterPlasmaVariables(pkg.get(), n_modes);
  RegisterDiagnostics(pkg.get());

  return pkg;
}

} // namespace Modes4D
