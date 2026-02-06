#include "pgen_harris4d.hpp"

#include "utils/error_checking.hpp"

namespace Modes4D {

void InitializeHarrisModes(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {
  (void)pmb;
  (void)pin;
  PARTHENON_FAIL(
      "harris_4d setup is scaffolded but not implemented yet. Next step is zero-mode "
      "Harris initialization plus A_y^(0) perturbation.");
}

} // namespace Modes4D
