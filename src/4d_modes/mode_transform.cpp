#include "mode_transform.hpp"

namespace Modes4D {

std::vector<double> ForwardTransform(const ModeTables &tables,
                                     const std::vector<double> &mode_values) {
  return tables.ModesToNodes(mode_values);
}

std::vector<double> InverseTransform(const ModeTables &tables,
                                     const std::vector<double> &node_values) {
  return tables.NodesToModes(node_values);
}

} // namespace Modes4D
