// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// system.h contains a list of properties needed to define a system of protein
// hard spheres

#ifndef PEN_SYSTEM_H
#define PEN_SYSTEM_H

#include "vec3.h"
#include <vector>

struct System {
  // positions of all beads at time t
  std::vector<Vec3> pos;
  // velocities of all beads at time t
  std::vector<Vec3> vel;
  // bead clocks
  std::vector<double> times;
  // collision counters of beads
  std::vector<uint64_t> counter;

  System(unsigned int nbeads)
      : pos(nbeads), vel(nbeads), times(nbeads), counter(nbeads) {}
};

#endif
