#include "plasma4d_modes.hpp"

#include <string>
#include <vector>

namespace Modes4D {
using namespace parthenon::package::prelude;

void RegisterPlasmaVariables(parthenon::StateDescriptor *pkg, const int n_modes) {
  constexpr int kSpecies = 2; // ions + electrons in the initial design.
  constexpr int kVarsPerModePerSpecies = 6; // rho, mom_x/y/z/w, energy

  std::vector<std::string> labels(kSpecies * kVarsPerModePerSpecies * n_modes);
  int idx = 0;
  for (int s = 0; s < kSpecies; ++s) {
    const std::string prefix = (s == 0) ? "ion" : "electron";
    for (int n = 0; n < n_modes; ++n) {
      labels[idx++] = prefix + "_rho_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momx_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momy_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momz_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_momw_mode_" + std::to_string(n);
      labels[idx++] = prefix + "_energy_mode_" + std::to_string(n);
    }
  }

  Metadata m({Metadata::Cell, Metadata::Independent, Metadata::FillGhost,
              Metadata::WithFluxes},
             std::vector<int>({static_cast<int>(labels.size())}), labels);
  pkg->AddField("plasma4d_cons", m);
}

} // namespace Modes4D
