#ifndef MODES4D_PACKAGE4D_MODES_HPP_
#define MODES4D_PACKAGE4D_MODES_HPP_

#include <memory>

#include <parthenon/package.hpp>

namespace Modes4D {

std::shared_ptr<parthenon::StateDescriptor> Initialize(parthenon::ParameterInput *pin);

} // namespace Modes4D

#endif // MODES4D_PACKAGE4D_MODES_HPP_
