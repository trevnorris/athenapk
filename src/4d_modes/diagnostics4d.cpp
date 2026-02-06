#include "diagnostics4d.hpp"

namespace Modes4D {
using namespace parthenon::package::prelude;

void RegisterDiagnostics(parthenon::StateDescriptor *pkg) {
  // Placeholder accumulators for the required energy/leakage ledger channels.
  pkg->AddParam<double>("diag/int_jw_ew", 0.0);
  pkg->AddParam<double>("diag/int_s_leak", 0.0);
}

} // namespace Modes4D
