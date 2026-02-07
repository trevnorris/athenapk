// AthenaPK - a performance portable block structured AMR MHD code
// Copyright (c) 2020-2026, Athena Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE")

// Parthenon headers
#include "mesh/mesh.hpp"
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../main.hpp"

namespace brio_wu {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  auto jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  auto kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  const Real rho_l = pin->GetOrAddReal("problem/brio_wu", "rho_l", 1.0);
  const Real pres_l = pin->GetOrAddReal("problem/brio_wu", "pres_l", 1.0);
  const Real vx_l = pin->GetOrAddReal("problem/brio_wu", "vx_l", 0.0);
  const Real vy_l = pin->GetOrAddReal("problem/brio_wu", "vy_l", 0.0);
  const Real vz_l = pin->GetOrAddReal("problem/brio_wu", "vz_l", 0.0);
  const Real bx_l = pin->GetOrAddReal("problem/brio_wu", "bx_l", 0.75);
  const Real by_l = pin->GetOrAddReal("problem/brio_wu", "by_l", 1.0);
  const Real bz_l = pin->GetOrAddReal("problem/brio_wu", "bz_l", 0.0);

  const Real rho_r = pin->GetOrAddReal("problem/brio_wu", "rho_r", 0.125);
  const Real pres_r = pin->GetOrAddReal("problem/brio_wu", "pres_r", 0.1);
  const Real vx_r = pin->GetOrAddReal("problem/brio_wu", "vx_r", 0.0);
  const Real vy_r = pin->GetOrAddReal("problem/brio_wu", "vy_r", 0.0);
  const Real vz_r = pin->GetOrAddReal("problem/brio_wu", "vz_r", 0.0);
  const Real bx_r = pin->GetOrAddReal("problem/brio_wu", "bx_r", 0.75);
  const Real by_r = pin->GetOrAddReal("problem/brio_wu", "by_r", -1.0);
  const Real bz_r = pin->GetOrAddReal("problem/brio_wu", "bz_r", 0.0);

  const Real x_discont = pin->GetOrAddReal("problem/brio_wu", "x_discont", 0.0);
  const Real gamma = pin->GetReal("hydro", "gamma");
  const Real gm1_inv = 1.0 / (gamma - 1.0);

  auto &mbd = pmb->meshblock_data.Get();
  auto &cons = mbd->Get("cons").data;
  auto &coords = pmb->coords;

  pmb->par_for(
      "Init brio_wu", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const bool left = coords.Xc<1>(i) < x_discont;

        const Real rho = left ? rho_l : rho_r;
        const Real pres = left ? pres_l : pres_r;
        const Real vx = left ? vx_l : vx_r;
        const Real vy = left ? vy_l : vy_r;
        const Real vz = left ? vz_l : vz_r;
        const Real bx = left ? bx_l : bx_r;
        const Real by = left ? by_l : by_r;
        const Real bz = left ? bz_l : bz_r;

        const Real mx = rho * vx;
        const Real my = rho * vy;
        const Real mz = rho * vz;
        const Real ekin = 0.5 * (mx * mx + my * my + mz * mz) / rho;
        const Real emag = 0.5 * (bx * bx + by * by + bz * bz);

        cons(IDN, k, j, i) = rho;
        cons(IM1, k, j, i) = mx;
        cons(IM2, k, j, i) = my;
        cons(IM3, k, j, i) = mz;
        cons(IB1, k, j, i) = bx;
        cons(IB2, k, j, i) = by;
        cons(IB3, k, j, i) = bz;
        cons(IEN, k, j, i) = pres * gm1_inv + ekin + emag;
      });
}
} // namespace brio_wu
