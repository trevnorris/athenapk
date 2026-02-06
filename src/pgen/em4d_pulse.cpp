//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2026, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include "pgen.hpp"

#include <cmath>

#include "../main.hpp"
#include "utils/error_checking.hpp"

namespace em4d_pulse {
using namespace parthenon::package::prelude;

namespace {
} // namespace

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  PARTHENON_REQUIRE(pin->GetString("hydro", "fluid") == "glmmhd",
                    "em4d_pulse requires hydro/fluid=glmmhd.");

  auto modes_pkg = pmb->packages.Get("modes4d");
  PARTHENON_REQUIRE(modes_pkg->Param<bool>("enabled"),
                    "em4d_pulse requires modes4d/enabled=true in the input file.");

  const Real rho0 = pin->GetOrAddReal("problem/em4d_pulse", "rho0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem/em4d_pulse", "p0", 1.0);
  const Real amp = modes_pkg->Param<double>("em4d_pulse/amplitude");
  const Real sigma = modes_pkg->Param<double>("em4d_pulse/sigma");
  const Real x0 = modes_pkg->Param<double>("em4d_pulse/x0");
  const Real gamma = pin->GetOrAddReal("hydro", "gamma", 5.0 / 3.0);

  PARTHENON_REQUIRE(sigma > 0.0, "problem/em4d_pulse/sigma must be > 0.");

  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto &rc = pmb->meshblock_data.Get();
  auto &cons_dev = rc->Get("cons").data;
  auto &em_a_dev = rc->Get("em4d_a").data;
  auto &em_pi_dev = rc->Get("em4d_pi").data;
  auto &plasma_dev = rc->Get("plasma4d_cons").data;

  auto cons = cons_dev.GetHostMirrorAndCopy();
  auto em_a = em_a_dev.GetHostMirrorAndCopy();
  auto em_pi = em_pi_dev.GetHostMirrorAndCopy();
  auto plasma = plasma_dev.GetHostMirrorAndCopy();
  auto &coords = pmb->coords;

  const Real gm1 = gamma - 1.0;
  const int n_modes = modes_pkg->Param<int>("n_modes");
  const int n_em_vars = em_a.GetDim(4);
  const int n_plasma_vars = plasma.GetDim(4);

  for (int k = kb.s; k <= kb.e; ++k) {
    for (int j = jb.s; j <= jb.e; ++j) {
      for (int i = ib.s; i <= ib.e; ++i) {
        const Real x = coords.Xc<1>(i);
        const Real pulse = amp * std::exp(-0.5 * std::pow((x - x0) / sigma, 2));

        cons(IDN, k, j, i) = rho0;
        cons(IM1, k, j, i) = 0.0;
        cons(IM2, k, j, i) = 0.0;
        cons(IM3, k, j, i) = 0.0;
        cons(IEN, k, j, i) = p0 / gm1;
        cons(IB1, k, j, i) = 0.0;
        cons(IB2, k, j, i) = 0.0;
        cons(IB3, k, j, i) = 0.0;
        cons(IPS, k, j, i) = 0.0;

        for (int v = 0; v < n_em_vars; ++v) {
          em_a(v, k, j, i) = 0.0;
          em_pi(v, k, j, i) = 0.0;
        }
        em_a(2, k, j, i) = pulse; // mode 0, A_y component

        for (int v = 0; v < n_plasma_vars; ++v) {
          plasma(v, k, j, i) = 0.0;
        }
        for (int s = 0; s < 2; ++s) {
          const int species_offset = s * 6 * n_modes;
          plasma(species_offset + 0, k, j, i) = 0.5 * rho0;
          plasma(species_offset + 5, k, j, i) = 0.5 * p0 / gm1;
        }
      }
    }
  }

  cons_dev.DeepCopy(cons);
  em_a_dev.DeepCopy(em_a);
  em_pi_dev.DeepCopy(em_pi);
  plasma_dev.DeepCopy(plasma);
}

void SourceUnsplit(MeshData<Real> *md, const parthenon::SimTime &tm, const Real dt) {
  auto modes_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("modes4d");

  const Real c_wave = modes_pkg->Param<double>("em4d/c_wave");
  const Real pulse_amp = modes_pkg->Param<double>("em4d_pulse/amplitude");
  const Real pulse_sigma = modes_pkg->Param<double>("em4d_pulse/sigma");
  const Real pulse_x0 = modes_pkg->Param<double>("em4d_pulse/x0");
  const Real t_new = tm.time + dt;

  const int num_blocks = md->NumBlocks();
  for (int b = 0; b < num_blocks; ++b) {
    auto &bd = md->GetBlockData(b);
    auto *pmb = bd->GetBlockPointer();

    auto &a_dev = bd->Get("em4d_a").data;
    auto &pi_dev = bd->Get("em4d_pi").data;
    auto a = a_dev.GetHostMirrorAndCopy();
    auto pi = pi_dev.GetHostMirrorAndCopy();

    IndexRange ib = bd->GetBoundsI(IndexDomain::interior);
    IndexRange jb = bd->GetBoundsJ(IndexDomain::interior);
    IndexRange kb = bd->GetBoundsK(IndexDomain::interior);

    auto &coords = pmb->coords;
    const Real x1min = pmb->pmy_mesh->mesh_size.xmin(X1DIR);
    const Real x1max = pmb->pmy_mesh->mesh_size.xmax(X1DIR);
    const Real lx = x1max - x1min;
    const Real sigma2 = pulse_sigma * pulse_sigma;
    const Real xc_raw = pulse_x0 + (c_wave * t_new);
    const Real xc = x1min + std::fmod((xc_raw - x1min) + 1000.0 * lx, lx);

    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          for (int v = 0; v < a.GetDim(4); ++v) {
            a(v, k, j, i) = 0.0;
            pi(v, k, j, i) = 0.0;
          }

          const Real x = coords.Xc<1>(i);
          Real dx = x - xc;
          if (dx > 0.5 * lx) dx -= lx;
          if (dx < -0.5 * lx) dx += lx;

          const Real ay = pulse_amp * std::exp(-0.5 * (dx * dx) / sigma2);
          const Real piy = c_wave * (dx / sigma2) * ay;

          a(2, k, j, i) = ay;  // mode 0, A_y
          pi(2, k, j, i) = piy;
        }
      }
    }

    a_dev.DeepCopy(a);
    pi_dev.DeepCopy(pi);
  }
}

} // namespace em4d_pulse
