#ifndef MODES4D_MODE_TRANSFORM_HPP_
#define MODES4D_MODE_TRANSFORM_HPP_

#include <vector>

#include "mode_tables.hpp"

namespace Modes4D {

// Thin wrappers around ModeTables transforms so call sites have a dedicated API.
std::vector<double> ForwardTransform(const ModeTables &tables,
                                     const std::vector<double> &mode_values);
std::vector<double> InverseTransform(const ModeTables &tables,
                                     const std::vector<double> &node_values);

} // namespace Modes4D

#endif // MODES4D_MODE_TRANSFORM_HPP_
