#include "em4d_modes.hpp"

#include <string>
#include <vector>

namespace Modes4D {
using namespace parthenon::package::prelude;

void RegisterEMVariables(parthenon::StateDescriptor *pkg, const int n_modes) {
  std::vector<std::string> labels(5 * n_modes);
  for (int n = 0; n < n_modes; ++n) {
    labels[5 * n + 0] = "a0_mode_" + std::to_string(n);
    labels[5 * n + 1] = "ax_mode_" + std::to_string(n);
    labels[5 * n + 2] = "ay_mode_" + std::to_string(n);
    labels[5 * n + 3] = "az_mode_" + std::to_string(n);
    labels[5 * n + 4] = "aw_mode_" + std::to_string(n);
  }

  Metadata m({Metadata::Cell, Metadata::Independent, Metadata::FillGhost,
              Metadata::WithFluxes},
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

  m = Metadata({Metadata::Cell, Metadata::Independent, Metadata::FillGhost,
                Metadata::WithFluxes},
               std::vector<int>({5 * n_modes}), pi_labels);
  pkg->AddField("em4d_pi", m);
}

} // namespace Modes4D
