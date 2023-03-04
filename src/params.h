// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PEN_PARAMS_H
#define PEN_PARAMS_H

#include "config.h"

struct Param {
  // mass of each bead
  double m;
  // mass of each solvent particle
  double m_sol;
  // diameter of each bead
  double sigma_bb;
  // radius of each bead
  double sigma_sb;
  // squared diameter
  double sigma_bb2;
  double sigma_sb2;
  // distance constraints for nearest neighbors
  double near_min;
  double near_max;
  // smallest bond length squared between nearest neighbors
  double near_min2;
  // largest bond length squared between nearest neighbors
  double near_max2;
  // distance constraints for next-nearest neighbors
  double nnear_min;
  double nnear_max;
  // smallest bond length squared between next-nearest neighbors
  double nnear_min2;
  // largest bond length squared between next-nearest neighbors
  double nnear_max2;
  // distance constraints for nonlocal beads
  double rh;
  double rc;
  // smallest bond length squared between nonlocal beads
  double rh2;
  // largest bond length squared between nonlocal beads
  double rc2;
  // nonlocal bond energy
  double eps;

  // vector of indices of beads which form bonds
  NonlocalBonds nonlocal_bonds;
  // vector of indices of beads which form bonds that can be broken
  NonlocalBonds transient_bonds;
  // transient_bonds placeholder for equilibrium run
  NonlocalBonds transient_bonds_eq;
  // vector of indices of beads which form bonds that cannot be broken
  NonlocalBonds permanent_bonds;

  // number of tries to place beads in box
  uint64_t tries;
  // number of beads
  unsigned int nbeads;
  // number of solvent particles
  unsigned int nsol;
  // length of box
  double length;
  // number of cells in each dimension of the box
  unsigned int ncell;
  // number of trajectories
  unsigned int total_iter;
  // number of time intervals
  unsigned int nsteps;
  unsigned int nsteps_eq;
  // step at which to write output to file
  unsigned int write_step;
  // increment time
  double del_t;
  // frequency of solvent collisions
  double del_t_coll;
  // seeds for random number generator
  std::vector<unsigned int> seeds;
  // Boltzmann constant
  double k_B;
  // temperature of the system
  double temp;
  // max number of nonlocal bonds that can form
  unsigned int max_nbonds;
  // solvent rotation angle
  unsigned int angle;
};

#endif
