#ifndef MODES4D_EM4D_MODES_HPP_
#define MODES4D_EM4D_MODES_HPP_

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

namespace Modes4D {

void RegisterEMVariables(parthenon::StateDescriptor *pkg, int n_modes,
                         bool enable_conservative_transport);
parthenon::TaskStatus AddEMTransportFluxes(parthenon::MeshData<parthenon::Real> *md);
void SourceUnsplit(parthenon::MeshData<parthenon::Real> *md,
                   const parthenon::SimTime &tm, const parthenon::Real dt);

} // namespace Modes4D

#endif // MODES4D_EM4D_MODES_HPP_
