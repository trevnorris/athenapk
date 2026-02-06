#ifndef MODES4D_EM4D_MODES_HPP_
#define MODES4D_EM4D_MODES_HPP_

#include <parthenon/package.hpp>

namespace Modes4D {

void RegisterEMVariables(parthenon::StateDescriptor *pkg, int n_modes);

} // namespace Modes4D

#endif // MODES4D_EM4D_MODES_HPP_
