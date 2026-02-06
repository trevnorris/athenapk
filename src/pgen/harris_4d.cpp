//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2026, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include "pgen.hpp"

#if ATHENAPK_ENABLE_4D_MODES
#include "../4d_modes/em4d_modes.hpp"
#include "../4d_modes/pgen_harris4d.hpp"
#endif
#include "utils/error_checking.hpp"

namespace harris_4d {
using namespace parthenon::package::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
#if ATHENAPK_ENABLE_4D_MODES
  PARTHENON_REQUIRE(pin->GetString("hydro", "fluid") == "glmmhd",
                    "harris_4d requires hydro/fluid=glmmhd.");
  auto modes_pkg = pmb->packages.Get("modes4d");
  PARTHENON_REQUIRE(modes_pkg->Param<bool>("enabled"),
                    "harris_4d requires modes4d/enabled=true in the input file.");

  Modes4D::InitializeHarrisModes(pmb, pin);
#else
  (void)pmb;
  (void)pin;
  PARTHENON_FAIL(
      "harris_4d requested but AthenaPK_ENABLE_4D_MODES=OFF at configure time.");
#endif
}

void SourceUnsplit(MeshData<Real> *md, const parthenon::SimTime &tm, const Real dt) {
#if ATHENAPK_ENABLE_4D_MODES
  Modes4D::SourceUnsplit(md, tm, dt);
#else
  (void)md;
  (void)tm;
  (void)dt;
#endif
}

} // namespace harris_4d
