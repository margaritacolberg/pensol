// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// hardspheres.cc initializes the positions and velocities of a hard sphere
// protein, and handles the initialization and processing of collision and cell
// crossing events of the hard spheres in a priority queue

#include "hardspheres.h"
#include <algorithm>
#include <cassert>
#include <cmath>

// calculate the x-, y- and z-coordinates of a point that lies on
// the surface of a unit sphere;
// sphere point picking: http://mathworld.wolfram.com/SpherePointPicking.html;
// see also Marsaglia, 1972
void unit_sphere(Random &mt, double &x, double &y, double &z) {
  std::uniform_real_distribution<> uniform_vw(-1.0, 1.0);
  double v, w, s;

  // reject points for which v^2 + w^2 >= 1
  do {
    v = uniform_vw(mt);
    w = uniform_vw(mt);
    s = v * v + w * w;
  } while (!(s < 1));

  // using the random numbers v and w, calculate the x-, y- and
  // z-coordinates
  x = 2 * v * std::sqrt(1 - v * v - w * w);
  y = 2 * w * std::sqrt(1 - v * v - w * w);
  z = 1 - 2 * (v * v + w * w);
}

// initial positions of all beads
void init_pos(std::vector<Vec3> &pos, const Box &box, Random &mt,
              const Param &p) {
  const unsigned int nbeads = pos.size();
  std::uniform_real_distribution<> uniform_near(0.0, 1.0);

  if (pos.size() < 1)
    return;

  // location of the first bead
  pos[0].x = 0.0;
  pos[0].y = 0.0;
  pos[0].z = 0.0;

  uint64_t tries = 0;
  unsigned int i = 1;

  const double near_min3 = std::pow(p.near_min, 3.0);
  const double near_max3 = std::pow(p.near_max, 3.0);

  while (i < nbeads) {
    if (!(tries++ < p.tries)) {
      throw std::logic_error(
          "number of tries exceeded while placing bead in box");
    }

    // calculate the x-, y- and z-coordinates of a point that
    // lies on the surface of a unit sphere
    double near_dx, near_dy, near_dz;
    unit_sphere(mt, near_dx, near_dy, near_dz);

    // pick a scaling factor that corresponds to the distance
    // constraint between the nearest beads
    const double val = uniform_near(mt);
    const double near = std::cbrt(near_min3 + (val * (near_max3 - near_min3)));

    // center the point at the position of the previous
    // point, and scale according to the distance constraint
    // between the nearest beads
    double x = near * near_dx + pos[i - 1].x;
    double y = near * near_dy + pos[i - 1].y;
    double z = near * near_dz + pos[i - 1].z;

    // if the next-nearest neighbor exists,
    if (i >= 2) {
      // calculate the distance between the new bead
      // and the next-nearest neighbor
      double nnear_dx = x - pos[i - 2].x;
      double nnear_dy = y - pos[i - 2].y;
      double nnear_dz = z - pos[i - 2].z;
      box.mindist(nnear_dx, nnear_dy, nnear_dz);

      // check if the distance constraint between the next-
      // nearest beads is satisfied
      double nnear_dist2 =
          nnear_dx * nnear_dx + nnear_dy * nnear_dy + nnear_dz * nnear_dz;
      if (nnear_dist2 < p.nnear_min2 || nnear_dist2 > p.nnear_max2)
        continue;
    }

    // declare a variable to check if there is an overlap of two beads
    bool overlap = false;

    for (unsigned int j = 0; j < i; j++) {
      double dx = x - pos[j].x;
      double dy = y - pos[j].y;
      double dz = z - pos[j].z;
      box.mindist(dx, dy, dz);

      const double dist2 = dx * dx + dy * dy + dz * dz;

      // check if the nonlocal beads overlap
      if (i >= j + 3) {
        if (dist2 < p.rh2) {
          // if there is an overlap,
          overlap = true;
          // exit without checking the other distances
          break;
        }
      }

      const Config bond_mask = p.transient_bonds.get_bond_mask(j, i);

      // check if there is a bond between two nonlocal beads
      if (bond_mask && dist2 < p.rc2) {
        // if there is a bond,
        overlap = true;
        // exit without checking the other distances
        break;
      }
    }

    if (overlap)
      continue;

    // store only the x-, y- and z-coordinates whose bead
    // is the allowable distance from its next-nearest neighbor
    pos[i].x = x;
    pos[i].y = y;
    pos[i].z = z;
    i++;
    tries = 0;
  }
}

// if two nonlocal beads which can form a transient bond happened to reach
// bonding distance at the end of a simulation in the equilibrium trajectory,
// and the configuration at the end of this simulation is being used as input to
// the bond forming trajectory, turn on the transient bond
void init_update_config(std::vector<Vec3> &pos, UpdateConfig &update_config,
                        const Box &box, const double rc2,
                        const NonlocalBonds &transient_bonds) {
  const unsigned int nbeads = pos.size();

  // reset bonds
  update_config = UpdateConfig();

  for (unsigned int i = 0; i + 3 < nbeads; i++) {
    for (unsigned int j = i + 3; j < nbeads; j++) {
      double dx = pos[j].x - pos[i].x;
      double dy = pos[j].y - pos[i].y;
      double dz = pos[j].z - pos[i].z;
      box.mindist(dx, dy, dz);

      const double dist2 = dx * dx + dy * dy + dz * dz;

      // determine if a bond between i and j will form
      const Config t_bond_mask = transient_bonds.get_bond_mask(i, j);

      if (t_bond_mask && dist2 < rc2) {
        assert(update_config.non_bonded(t_bond_mask));
        update_config.flip_bond(t_bond_mask);
        assert(update_config.bonded(t_bond_mask));
      }
    }
  }
}

void shift_vel_CM_to_0(std::vector<Vec3> &vel) {
  const unsigned int nbeads = vel.size();

  double vel_x_tot = 0.0, vel_y_tot = 0.0, vel_z_tot = 0.0;
  for (unsigned int i = 0; i < nbeads; i++) {
    vel_x_tot += vel[i].x;
    vel_y_tot += vel[i].y;
    vel_z_tot += vel[i].z;
  }

  double vel_x_CM = vel_x_tot / double(nbeads);
  double vel_y_CM = vel_y_tot / double(nbeads);
  double vel_z_CM = vel_z_tot / double(nbeads);

  // to account for finite size of system, shift velocities
  // to yield a center of mass velocity of zero
  for (unsigned int i = 0; i < nbeads; i++) {
    vel[i].x -= vel_x_CM;
    vel[i].y -= vel_y_CM;
    vel[i].z -= vel_z_CM;
  }
}

// initial velocities of all beads from a Maxwell-Boltzmann distribution
void init_vel(std::vector<Vec3> &vel, Random &mt, const double temp,
              const double m) {
  const unsigned int nbeads = vel.size();

  const double mean = 0.0;
  const double sigma = std::sqrt(temp / m);
  std::normal_distribution<> gaussian{mean, sigma};

  for (unsigned int i = 0; i < nbeads; i++) {
    vel[i].x = gaussian(mt);
    vel[i].y = gaussian(mt);
    vel[i].z = gaussian(mt);
  }

  shift_vel_CM_to_0(vel);
}

// if an inner collision occurs (two beads reach the minimum
// allowable separation), the function returns true and returns the
// time until collision in t; if a collision does not occur, the
// function returns false
bool t_until_inner_coll(double xij, double yij, double zij, double vxij,
                        double vyij, double vzij, const double sigma2,
                        double &t) {
  double rvij = xij * vxij + yij * vyij + zij * vzij;

  // ! accounts for comparison with NaN values
  if (!(rvij < 0)) {
    return false;
  }

  double vvij = vxij * vxij + vyij * vyij + vzij * vzij;
  double rrij = xij * xij + yij * yij + zij * zij;
  double rad1 = vvij * (rrij - sigma2);
  double rvij2 = rvij * rvij;

  // ! accounts for comparison with NaN values
  if (!(rvij2 > rad1)) {
    return false;
  }

  t += (-rvij - std::sqrt(rvij2 - rad1)) / vvij;

  return true;
}

// if an outer collision occurs (two beads reach the maximum allowable
// separation), the function returns the time until collision in t
void t_until_outer_coll(double xij, double yij, double zij, double vxij,
                        double vyij, double vzij, const double sigma2,
                        double &t) {
  double rvij = xij * vxij + yij * vyij + zij * vzij;
  double vvij = vxij * vxij + vyij * vyij + vzij * vzij;
  double rrij = xij * xij + yij * yij + zij * zij;
  double rad1 = vvij * (rrij - sigma2);
  double rvij2 = rvij * rvij;

  t += (-rvij + std::sqrt(rvij2 - rad1)) / vvij;
}

// if a collision occurs, the function returns the velocities
// of the colliding beads post-collision
bool v_after_coll(double xij, double yij, double zij, const double sigma2,
                  double &vxi, double &vyi, double &vzi, double &vxj,
                  double &vyj, double &vzj, const std::optional<double> &dU,
                  const double m) {
  double vxij = vxj - vxi;
  double vyij = vyj - vyi;
  double vzij = vzj - vzi;
  const double mu = 0.5 * m;
  double rvij = (xij * vxij + yij * vyij + zij * vzij) / sigma2;
  const double rvij2 = rvij * rvij;
  double dU_norm;

  bool is_transmission;
  if (dU && rvij2 > (dU_norm = 2 * (*dU) / (sigma2 * mu))) {
    // transmission with change in energy
    rvij -= std::copysign(std::sqrt(rvij2 - dU_norm), rvij);
    LOG_DEBUG("event is a transmission with change in energy");
    is_transmission = true;
  } else {
    // reflection without change in energy
    rvij += rvij;
    LOG_DEBUG("event is a reflection without change in energy");
    is_transmission = false;
  }

  double pxij = rvij * xij;
  double pyij = rvij * yij;
  double pzij = rvij * zij;

  vxi += pxij / 2;
  vyi += pyij / 2;
  vzi += pzij / 2;
  vxj -= pxij / 2;
  vyj -= pyij / 2;
  vzj -= pzij / 2;

  return is_transmission;
}

// split the box into smaller cells, and store beads in each cell
void init_cells(const std::vector<Vec3> &pos, const Box &box, Cells &cells) {
  const unsigned int nbeads = pos.size();

  // clear old particle indices
  for (auto &cell : cells.cells) {
    cell.clear();
  }

  // number the beads
  for (unsigned int j = 0; j < nbeads; j++) {
    // calculate the 3D indices of the cells, based on the x-, y- and
    // z-coordinates of the beads;
    // clamp cell index to allowed interval [0, ncell - 1] to account
    // for floating-point rounding, eg. x = -0.00001 is icell_x = 0,
    // not -1
    double x = pos[j].x;
    double y = pos[j].y;
    double z = pos[j].z;

    box.minpos(x, y, z);
    unsigned int icell_x =
        std::clamp(int(cells.inv_lcell * x), 0, int(cells.ncell - 1));
    unsigned int icell_y =
        std::clamp(int(cells.inv_lcell * y), 0, int(cells.ncell - 1));
    unsigned int icell_z =
        std::clamp(int(cells.inv_lcell * z), 0, int(cells.ncell - 1));
    // convert the 3D indices to a single index for each cell
    unsigned int icell =
        icell_x + icell_y * cells.ncell + icell_z * cells.ncell * cells.ncell;

    // store the 3D indices of the cell
    cells.cell_x.emplace_back(icell_x);
    cells.cell_y.emplace_back(icell_y);
    cells.cell_z.emplace_back(icell_z);
    LOG_DEBUG("particle " << j << " at " << pos[j].x << ", " << pos[j].y << ", "
                          << pos[j].z << " in cell " << icell_x << ", "
                          << icell_y << ", " << icell_z);
    // store the 1D index for each cell, and the beads
    // contained within, in a vector of vectors
    cells.cells[icell].emplace_back(j);
  }
}

// calculate the change in clocks, position, and velocity of
// two beads
void delta_tpv(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
               const Box &box, const std::vector<double> &times, unsigned int i,
               unsigned int j, double &t, double &dx, double &dy, double &dz,
               double &dvx, double &dvy, double &dvz) {
  if (times[i] < times[j]) {
    std::swap(i, j);
  }

  t = times[i];
  double dt = t - times[j];
  assert(dt >= 0);

  dx = pos[j].x + vel[j].x * dt - pos[i].x;
  dy = pos[j].y + vel[j].y * dt - pos[i].y;
  dz = pos[j].z + vel[j].z * dt - pos[i].z;
  box.mindist(dx, dy, dz);

  dvx = vel[j].x - vel[i].x;
  dvy = vel[j].y - vel[i].y;
  dvz = vel[j].z - vel[i].z;
}

// if a collision between two beads occurs, calculate the absolute
// time of that collision, and add an event to the priority queue which
// contains the absolute time of the collision, indices of the colliding
// particles, and their current collision counters
void if_coll(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
             const double rh2, const double rc2, const Box &box,
             const std::vector<uint64_t> &counter, EventQueue &event_queue,
             const std::vector<double> &times, unsigned int i, unsigned int j,
             const NonlocalBonds &transient_bonds,
             const NonlocalBonds &permanent_bonds,
             const UpdateConfig &update_config, const unsigned int max_nbonds) {
  double t, dx, dy, dz, dvx, dvy, dvz;

  if (i > j)
    std::swap(i, j);

  if (j == i || j == i + 1 || j == i + 2)
    return;

  delta_tpv(pos, vel, box, times, i, j, t, dx, dy, dz, dvx, dvy, dvz);

  // determine if a bond between i and j will form
  const Config t_bond_mask = transient_bonds.get_bond_mask(i, j);
  const Config p_bond_mask = permanent_bonds.get_bond_mask(i, j);

  // count the number of transient bonds already present
  const unsigned int nbonds = update_config.count_bonds();
  assert(nbonds <= max_nbonds);

  // if two nonlocal beads collide and form a bond
  // (i and j must initially not be bonded),
  if (t_bond_mask && (nbonds < max_nbonds) &&
      update_config.non_bonded(t_bond_mask) &&
      t_until_inner_coll(dx, dy, dz, dvx, dvy, dvz, rc2, t)) {
    MaxNonlocalInnerEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
    // else if two nonlocal beads collide elastically,
  } else if (t_until_inner_coll(dx, dy, dz, dvx, dvy, dvz, rh2, t)) {
    MinNonlocalInnerEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
    // else if two nonlocal bonded beads reach rc,
  } else if (update_config.bonded(t_bond_mask) || p_bond_mask) {
    t_until_outer_coll(dx, dy, dz, dvx, dvy, dvz, rc2, t);
    MaxNonlocalOuterEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
  }
}

// iterate collisions over all particle pairs
void iterate_coll(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
                  const double rh2, const double rc2, const Box &box,
                  const std::vector<uint64_t> &counter, EventQueue &event_queue,
                  const std::vector<double> &times, const Cells &cells,
                  unsigned int icell, unsigned int i,
                  const NonlocalBonds &transient_bonds,
                  const NonlocalBonds &permanent_bonds,
                  const UpdateConfig &update_config,
                  const unsigned int max_nbonds) {
  for (unsigned int j : cells.cells[icell]) {
    if_coll(pos, vel, rh2, rc2, box, counter, event_queue, times, i, j,
            transient_bonds, permanent_bonds, update_config, max_nbonds);
  }
}

void walls_of_cell(unsigned int &zmin, unsigned int &zmax, unsigned int &ymin,
                   unsigned int &ymax, unsigned int &xmin, unsigned int &xmax,
                   const Cells &cells, unsigned int i) {
  zmin = cells.cell_z[i] - 1 + cells.ncell;
  zmax = cells.cell_z[i] + 1 + cells.ncell;
  ymin = cells.cell_y[i] - 1 + cells.ncell;
  ymax = cells.cell_y[i] + 1 + cells.ncell;
  xmin = cells.cell_x[i] - 1 + cells.ncell;
  xmax = cells.cell_x[i] + 1 + cells.ncell;
}

// calculate the collisions occurring in each cell and its neighboring
// cells, but without any cell crossings
void add_events_for_one_bead(
    const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
    const double rh2, const double rc2, const Box &box,
    const std::vector<uint64_t> &counter, EventQueue &event_queue,
    const std::vector<double> &times, const Cells &cells, unsigned int i,
    const NonlocalBonds &transient_bonds, const NonlocalBonds &permanent_bonds,
    const UpdateConfig &update_config, const unsigned int max_nbonds) {
  unsigned int zmin, zmax, ymin, ymax, xmin, xmax;
  walls_of_cell(zmin, zmax, ymin, ymax, xmin, xmax, cells, i);

  // calculate the 3D indices of a cell and its neighbors
  for (unsigned int z = zmin; z <= zmax; z++) {
    for (unsigned int y = ymin; y <= ymax; y++) {
      for (unsigned int x = xmin; x <= xmax; x++) {
        unsigned int icell_x = x % cells.ncell;
        unsigned int icell_y = y % cells.ncell;
        unsigned int icell_z = z % cells.ncell;

        // convert the 3D indices to a 1D index
        unsigned int icell = icell_x + icell_y * cells.ncell +
                             icell_z * cells.ncell * cells.ncell;
        // iterate collisions over the cell
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }
}

// fill priority queue with collisions of all particle pairs in
// each cell and its neighboring cells
void add_events_for_all_beads(
    const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
    unsigned int nbeads, const double rh2, const double rc2, const Box &box,
    const std::vector<uint64_t> &counter, EventQueue &event_queue,
    const std::vector<double> &times, const Cells &cells,
    const NonlocalBonds &transient_bonds, const NonlocalBonds &permanent_bonds,
    const UpdateConfig &update_config, const unsigned int max_nbonds) {
  assert(pos.size() == nbeads);
  // for a collision between two beads,
  for (unsigned int i = 0; i < nbeads; i++) {
    add_events_for_one_bead(pos, vel, rh2, rc2, box, counter, event_queue,
                            times, cells, i, transient_bonds, permanent_bonds,
                            update_config, max_nbonds);
  }
}

// if a cell crossing occurs, the function returns true and returns the
// time until cell crossing in t, the wall crossed, and the x-, y- and
// z-indices of the new cell; if a cell crossing does not occur, the
// function returns false
bool t_until_cell(double xi, double yi, double zi, const Box &box,
                  const Cells &cells, double vxi, double vyi, double vzi,
                  unsigned int ixo, unsigned int iyo, unsigned int izo,
                  double &dt, BeadCellEvent::Wall &wall, unsigned int &ixn,
                  unsigned int &iyn, unsigned int &izn) {
  bool crossing = false;

  // position of left cell wall
  double x_wall = cells.lcell * ixo;
  // position of bottom cell wall
  double y_wall = cells.lcell * iyo;
  // position of back cell wall
  double z_wall = cells.lcell * izo;

  double dx = x_wall - xi;
  double dy = y_wall - yi;
  double dz = z_wall - zi;
  box.mindist(dx, dy, dz);

  // if the velocity of the bead is positive in x-direction,
  if (vxi > 0) {
    // calculate the time needed to reach the right wall
    double txr = (dx + cells.lcell) / vxi;

    // if a cell crossing occurs,
    if (!crossing || txr < dt) {
      crossing = true;
      dt = txr;
      wall = BeadCellEvent::xpos;
      ixn = (ixo + 1) % cells.ncell;
      iyn = iyo;
      izn = izo;
    }
  }

  // else if the velocity of the bead is negative in x-direction,
  else if (vxi < 0) {
    // calculate the time needed to reach the left wall
    double txl = dx / vxi;

    if (!crossing || txl < dt) {
      crossing = true;
      dt = txl;
      wall = BeadCellEvent::xneg;
      ixn = (ixo - 1 + cells.ncell) % cells.ncell;
      iyn = iyo;
      izn = izo;
    }
  }

  // if the velocity of the bead is positive in y-direction,
  if (vyi > 0) {
    // calculate the time needed to reach the top wall
    double tyr = (dy + cells.lcell) / vyi;

    if (!crossing || tyr < dt) {
      crossing = true;
      dt = tyr;
      wall = BeadCellEvent::ypos;
      ixn = ixo;
      iyn = (iyo + 1) % cells.ncell;
      izn = izo;
    }
  }

  // else if the velocity of the bead is negative in y-direction,
  else if (vyi < 0) {
    // calculate the time needed to reach the bottom wall
    double tyl = dy / vyi;

    if (!crossing || tyl < dt) {
      crossing = true;
      dt = tyl;
      wall = BeadCellEvent::yneg;
      ixn = ixo;
      iyn = (iyo - 1 + cells.ncell) % cells.ncell;
      izn = izo;
    }
  }

  // if the velocity of the bead is positive in z-direction,
  if (vzi > 0) {
    // calculate the time needed to reach the front wall
    double tzr = (dz + cells.lcell) / vzi;

    if (!crossing || tzr < dt) {
      crossing = true;
      dt = tzr;
      wall = BeadCellEvent::zpos;
      ixn = ixo;
      iyn = iyo;
      izn = (izo + 1) % cells.ncell;
    }
  }

  // else if the velocity of the bead is negative in z-direction,
  else if (vzi < 0) {
    // calculate the time needed to reach the back wall
    double tzl = dz / vzi;

    if (!crossing || tzl < dt) {
      crossing = true;
      dt = tzl;
      wall = BeadCellEvent::zneg;
      ixn = ixo;
      iyn = iyo;
      izn = (izo - 1 + cells.ncell) % cells.ncell;
    }
  }

  return crossing;
}

// check if a cell crossing of a bead will occur, and if it does,
// calculate the absolute time of that cell crossing, and add an event
// to the priority queue which contains the absolute time of the cell
// crossing, index of cell crossing particle, its current collision
// counter, indices of the old and new cells, and the cell wall that was
// crossed
void if_cell(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
             const Box &box, const std::vector<uint64_t> &counter,
             EventQueue &event_queue, const std::vector<double> &times,
             unsigned int i, const Cells &cells) {
  double dt;
  BeadCellEvent::Wall wall;
  unsigned int ixn, iyn, izn;

  if (t_until_cell(pos[i].x, pos[i].y, pos[i].z, box, cells, vel[i].x, vel[i].y,
                   vel[i].z, cells.cell_x[i], cells.cell_y[i], cells.cell_z[i],
                   dt, wall, ixn, iyn, izn)) {
    BeadCellEvent ev{times[i] + dt, i, counter[i], ixn, iyn, izn, wall};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
  }
}

// initialize priority queue with cell crossings of all particles
void init_cell_events(const std::vector<Vec3> &pos,
                      const std::vector<Vec3> &vel, unsigned int nbeads,
                      const Box &box, const std::vector<uint64_t> &counter,
                      EventQueue &event_queue, const std::vector<double> &times,
                      const Cells &cells) {
  for (unsigned int i = 0; i < nbeads; i++) {
    if_cell(pos, vel, box, counter, event_queue, times, i, cells);
  }
}

// move a cell crossing bead from the old cell to the new cell
void move_to_new_cell(Cells &cells, unsigned int i, unsigned int ixn,
                      unsigned int iyn, unsigned int izn) {
  // calculate the 1D index of the old cell from the 2D indices
  unsigned int old_cell = cells.cell_x[i] + cells.cell_y[i] * cells.ncell +
                          cells.cell_z[i] * cells.ncell * cells.ncell;
  assert(old_cell < cells.cells.size());

  // calculate the 1D index of the new cell from the 2D indices
  unsigned int new_cell =
      ixn + iyn * cells.ncell + izn * cells.ncell * cells.ncell;
  assert(new_cell < cells.cells.size());

  LOG_DEBUG("move particle " << i << " from " << cells.cell_x[i] << ", "
                             << cells.cell_y[i] << ", " << cells.cell_z[i]
                             << " to " << ixn << ", " << iyn << ", " << izn);
  auto ind =
      std::find(cells.cells[old_cell].begin(), cells.cells[old_cell].end(), i);
  assert(ind != cells.cells[old_cell].end());

  // move the cell crossing bead from the old cell to the
  // new cell
  cells.cells[new_cell].emplace_back(i);
  // erase the cell crossing bead from the list of indices
  // of the beads stored in the old cell
  std::swap(*ind, cells.cells[old_cell].back());
  cells.cells[old_cell].pop_back();

  cells.cell_x[i] = ixn;
  cells.cell_y[i] = iyn;
  cells.cell_z[i] = izn;
}

// update priority queue with collisions of all particle pairs involving
// a particle which recently crossed into a new cell
void add_events_for_bead_after_crossing(
    const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
    const double rh2, const double rc2, const Box &box,
    const std::vector<uint64_t> &counter, EventQueue &event_queue,
    const std::vector<double> &times, const Cells &cells_bead, unsigned int i,
    BeadCellEvent::Wall wall, const NonlocalBonds &transient_bonds,
    const NonlocalBonds &permanent_bonds, const UpdateConfig &update_config,
    const unsigned int max_nbonds) {
  unsigned int zmin, zmax, ymin, ymax, xmin, xmax;
  walls_of_cell(zmin, zmax, ymin, ymax, xmin, xmax, cells_bead, i);

  // if the wall in the positive x-direction is crossed,
  if (wall == BeadCellEvent::xpos) {
    // calculate the 3D indices of the 9 new neighbor cells
    unsigned int icell_x = xmax % cells_bead.ncell;

    for (unsigned int z = zmin; z <= zmax; z++) {
      for (unsigned int y = ymin; y <= ymax; y++) {
        unsigned int icell_y = y % cells_bead.ncell;
        unsigned int icell_z = z % cells_bead.ncell;

        // convert the 3D indices to a 1D index
        unsigned int icell = icell_x + icell_y * cells_bead.ncell +
                             icell_z * cells_bead.ncell * cells_bead.ncell;
        // iterate collisions over the cell
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells_bead, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }

  if (wall == BeadCellEvent::xneg) {
    unsigned int icell_x = xmin % cells_bead.ncell;

    for (unsigned int z = zmin; z <= zmax; z++) {
      for (unsigned int y = ymin; y <= ymax; y++) {
        unsigned int icell_y = y % cells_bead.ncell;
        unsigned int icell_z = z % cells_bead.ncell;

        unsigned int icell = icell_x + icell_y * cells_bead.ncell +
                             icell_z * cells_bead.ncell * cells_bead.ncell;
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells_bead, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }

  if (wall == BeadCellEvent::ypos) {
    unsigned int icell_y = ymax % cells_bead.ncell;

    for (unsigned int z = zmin; z <= zmax; z++) {
      for (unsigned int x = xmin; x <= xmax; x++) {
        unsigned int icell_x = x % cells_bead.ncell;
        unsigned int icell_z = z % cells_bead.ncell;

        unsigned int icell = icell_x + icell_y * cells_bead.ncell +
                             icell_z * cells_bead.ncell * cells_bead.ncell;
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells_bead, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }

  if (wall == BeadCellEvent::yneg) {
    unsigned int icell_y = ymin % cells_bead.ncell;

    for (unsigned int z = zmin; z <= zmax; z++) {
      for (unsigned int x = xmin; x <= xmax; x++) {
        unsigned int icell_x = x % cells_bead.ncell;
        unsigned int icell_z = z % cells_bead.ncell;

        unsigned int icell = icell_x + icell_y * cells_bead.ncell +
                             icell_z * cells_bead.ncell * cells_bead.ncell;
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells_bead, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }

  if (wall == BeadCellEvent::zpos) {
    unsigned int icell_z = zmax % cells_bead.ncell;

    for (unsigned int y = ymin; y <= ymax; y++) {
      for (unsigned int x = xmin; x <= xmax; x++) {
        unsigned int icell_x = x % cells_bead.ncell;
        unsigned int icell_y = y % cells_bead.ncell;

        unsigned int icell = icell_x + icell_y * cells_bead.ncell +
                             icell_z * cells_bead.ncell * cells_bead.ncell;
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells_bead, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }

  if (wall == BeadCellEvent::zneg) {
    unsigned int icell_z = zmin % cells_bead.ncell;

    for (unsigned int y = ymin; y <= ymax; y++) {
      for (unsigned int x = xmin; x <= xmax; x++) {
        unsigned int icell_x = x % cells_bead.ncell;
        unsigned int icell_y = y % cells_bead.ncell;

        unsigned int icell = icell_x + icell_y * cells_bead.ncell +
                             icell_z * cells_bead.ncell * cells_bead.ncell;
        iterate_coll(pos, vel, rh2, rc2, box, counter, event_queue, times,
                     cells_bead, icell, i, transient_bonds, permanent_bonds,
                     update_config, max_nbonds);
      }
    }
  }
}

// if two nearest neighboring beads within a bond reach minimum or maximum
// allowable separation, calculate the absolute time of that collision, and add
// an event to the priority queue which contains the absolute time of the
// collision, indices of the colliding particles, their current collision
// counters, and the type of collision they experienced
void if_nearest_bond(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
                     const Box &box, const std::vector<uint64_t> &counter,
                     EventQueue &event_queue, const std::vector<double> &times,
                     unsigned int i, unsigned int j, double near_min2,
                     double near_max2) {
  double t, dx, dy, dz, dvx, dvy, dvz;
  delta_tpv(pos, vel, box, times, i, j, t, dx, dy, dz, dvx, dvy, dvz);
  assert(check_bond(near_min2, near_max2, dx, dy, dz));

  if (t_until_inner_coll(dx, dy, dz, dvx, dvy, dvz, near_min2, t)) {
    MinNearestEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
  } else {
    t_until_outer_coll(dx, dy, dz, dvx, dvy, dvz, near_max2, t);
    MaxNearestEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
  }
}

// fill priority queue with nearest bond events between all particle
// pairs
void init_nearest_bond_events(const std::vector<Vec3> &pos,
                              const std::vector<Vec3> &vel, unsigned int nbeads,
                              const Box &box,
                              const std::vector<uint64_t> &counter,
                              EventQueue &event_queue,
                              const std::vector<double> &times,
                              double near_min2, double near_max2) {
  for (unsigned int i = 0; i + 1 < nbeads; i++) {
    if_nearest_bond(pos, vel, box, counter, event_queue, times, i, i + 1,
                    near_min2, near_max2);
  }
}

// if two next-nearest neighboring beads within a bond reach minimum or maximum
// allowable separation, calculate the absolute time of that collision, and add
// an event to the priority queue which contains the absolute time of the
// collision, indices of the colliding particles, their current collision
// counters, and the type of collision they experienced
void if_nnearest_bond(const std::vector<Vec3> &pos,
                      const std::vector<Vec3> &vel, const Box &box,
                      const std::vector<uint64_t> &counter,
                      EventQueue &event_queue, const std::vector<double> &times,
                      unsigned int i, unsigned int j, double nnear_min2,
                      double nnear_max2) {
  double t, dx, dy, dz, dvx, dvy, dvz;
  delta_tpv(pos, vel, box, times, i, j, t, dx, dy, dz, dvx, dvy, dvz);
  assert(check_bond(nnear_min2, nnear_max2, dx, dy, dz));

  if (t_until_inner_coll(dx, dy, dz, dvx, dvy, dvz, nnear_min2, t)) {
    MinNNearestEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
  } else {
    t_until_outer_coll(dx, dy, dz, dvx, dvy, dvz, nnear_max2, t);
    MaxNNearestEvent ev{t, i, j, counter[i], counter[j]};
    LOG_DEBUG("queueing " << ev);
    event_queue.emplace(ev);
  }
}

// fill priority queue with next-nearest bond events between all
// particle pairs
void init_nnearest_bond_events(const std::vector<Vec3> &pos,
                               const std::vector<Vec3> &vel,
                               unsigned int nbeads, const Box &box,
                               const std::vector<uint64_t> &counter,
                               EventQueue &event_queue,
                               const std::vector<double> &times,
                               double nnear_min2, double nnear_max2) {
  for (unsigned int i = 0; i + 2 < nbeads; i++) {
    if_nnearest_bond(pos, vel, box, counter, event_queue, times, i, i + 2,
                     nnear_min2, nnear_max2);
  }
}

// check if there is an overlap of two beads
bool check_overlap(const unsigned int i, const std::vector<Vec3> &pos,
                   const std::vector<Vec3> &vel,
                   const std::vector<double> &times, const double sigma2,
                   const Box &box) {
  const unsigned int nbeads = pos.size();

  const double tol = 100000 * std::numeric_limits<double>::epsilon();
  const double sigma2_tol = sigma2 * (1 - tol);

  for (unsigned int j = 0; j < nbeads; j++) {
    if (j == i || j == i + 1 || j == i - 1 || j == i + 2 || j == i - 2)
      continue;

    const double dt = times[i] - times[j];
    double xi = pos[i].x;
    double yi = pos[i].y;
    double zi = pos[i].z;
    double xj = pos[j].x + vel[j].x * dt;
    double yj = pos[j].y + vel[j].y * dt;
    double zj = pos[j].z + vel[j].z * dt;
    double dx = xj - xi;
    double dy = yj - yi;
    double dz = zj - zi;
    box.mindist(dx, dy, dz);

    double dist2 = dx * dx + dy * dy + dz * dz;
    // ! accounts for comparison with NaN values
    if (!(dist2 > sigma2_tol)) {
      // print coordinates of centers of overlapping beads
      LOG_DEBUG("particle " << i << " at " << xi << ", " << yi << ", " << zi);
      LOG_DEBUG("particle " << j << " at " << xj << ", " << yj << ", " << zj);
      LOG_DEBUG("distance " << std::setprecision(15) << std::sqrt(dist2));
      return false;
    }
  }

  return true;
}

// check if bond conditions are satisfied
bool check_bond(const double min_dist2, const double max_dist2, double dx,
                double dy, double dz) {
  // after the positions of two colliding beads are updated, the bond
  // distance may be slightly smaller / larger than the minimum / maximum
  // bond distance due to floating point rounding; to account for this, we
  // compare the bond distance by including a small tolerance that is a
  // multiple of the machine epsilon
  //
  // to compare the squared minimum/maximum bond distance, we must also
  // square the tolerance: (1 ± ε')^2 = 1 ± 2 * ε' + ε'^2 ≈ 1 ± 2 * ε' for
  // ε' << 1, e.g., if ε' = 5 * ε, tol is 10 * ε
  const double tol = 100000000 * std::numeric_limits<double>::epsilon();
  const double min_dist2_tol = min_dist2 * (1 - tol);
  const double max_dist2_tol = max_dist2 * (1 + tol);

  double dist2 = dx * dx + dy * dy + dz * dz;
  // ! accounts for comparison with NaN values
  if (!(dist2 > min_dist2_tol && dist2 < max_dist2_tol)) {
    LOG_DEBUG("diff b/w dist and tol: min = "
              << std::sqrt(dist2) - std::sqrt(min_dist2_tol)
              << " and max = " << std::sqrt(max_dist2_tol) - std::sqrt(dist2));
    return false;
  }

  return true;
}

bool check_local_dist(const std::vector<Vec3> &pos, const Box &box,
                      const double near_min2, const double near_max2,
                      const double nnear_min2, const double nnear_max2) {
  const unsigned int nbeads = pos.size();

  for (unsigned int i = 0; i + 1 < nbeads; i++) {
    double dx = pos[i + 1].x - pos[i].x;
    double dy = pos[i + 1].y - pos[i].y;
    double dz = pos[i + 1].z - pos[i].z;
    box.mindist(dx, dy, dz);

    if (!check_bond(near_min2, near_max2, dx, dy, dz)) {
      LOG_DEBUG("nearest neighbors " << i << " and " << i + 1 << " overlap");
      return false;
    }
  }

  for (unsigned int i = 0; i + 2 < nbeads; i++) {
    double dx = pos[i + 2].x - pos[i].x;
    double dy = pos[i + 2].y - pos[i].y;
    double dz = pos[i + 2].z - pos[i].z;
    box.mindist(dx, dy, dz);

    if (!check_bond(nnear_min2, nnear_max2, dx, dy, dz)) {
      LOG_DEBUG("next-nearest neighbors " << i << " and " << i + 2
                                          << " overlap");
      return false;
    }
  }

  return true;
}

bool check_nonlocal_dist(const std::vector<Vec3> &pos, const Box &box,
                         const double rc2, const double rh2,
                         const NonlocalBonds &permanent_bonds) {
  const unsigned int nbeads = pos.size();

  for (unsigned int i = 0; i < nbeads - 3; i++) {
    for (unsigned int j = i + 3; j < nbeads; j++) {
      double dx = pos[i].x - pos[j].x;
      double dy = pos[i].y - pos[j].y;
      double dz = pos[i].z - pos[j].z;
      box.mindist(dx, dy, dz);

      const double dist2 = dx * dx + dy * dy + dz * dz;

      if (!(dist2 > rh2)) {
        LOG_DEBUG("nonlocal beads " << i << " and " << j
                                    << " overlap with distance "
                                    << std::sqrt(dist2));
        return false;
      }

      const Config p_bond_mask = permanent_bonds.get_bond_mask(i, j);

      if (p_bond_mask && !(dist2 < rc2)) {
        LOG_DEBUG("permanent bond between beads "
                  << i << " and " << j << " separated by distance "
                  << std::sqrt(dist2) << " has been broken");
        return false;
      }
    }
  }

  return true;
}

// update positions of the colliding beads post-collision
void update_pos(Vec3 &pos, Vec3 &vel, double &time, double update_time) {
  const double dt = update_time - time;
  pos.x += vel.x * dt;
  pos.y += vel.y * dt;
  pos.z += vel.z * dt;

  time = update_time;
}

// update velocities of the colliding beads post-collision
bool final_vel(const std::vector<Vec3> &pos, std::vector<Vec3> &vel,
               const double sigma2, unsigned int i, unsigned int j,
               const Box &box, const std::optional<double> &dU,
               const double m) {
  double dx = pos[j].x - pos[i].x;
  double dy = pos[j].y - pos[i].y;
  double dz = pos[j].z - pos[i].z;
  box.mindist(dx, dy, dz);

  return v_after_coll(dx, dy, dz, sigma2, vel[i].x, vel[i].y, vel[i].z,
                      vel[j].x, vel[j].y, vel[j].z, dU, m);
}

double total_k_E(const std::vector<Vec3> &vel, const double m) {
  const unsigned int nbeads = vel.size();

  double k_E = 0.0;
  for (unsigned int i = 0; i < nbeads; i++) {
    double vx2 = vel[i].x * vel[i].x;
    double vy2 = vel[i].y * vel[i].y;
    double vz2 = vel[i].z * vel[i].z;

    k_E += 0.5 * m * (vx2 + vy2 + vz2);
  }
  std::cout << k_E << " kinetic energy of beads" << std::endl;

  return k_E;
}

double compute_hamiltonian(const std::vector<Vec3> &vel, const double m,
                           UpdateConfig &update_config) {
  return total_k_E(vel, m) + update_config.U;
}

void dist_between_nonlocal_beads(const std::vector<Vec3> &pos, const Box &box,
                                 const NonlocalBonds &nonlocal_bonds,
                                 std::vector<double> &dist) {
  const unsigned int nbeads = pos.size();

  dist.clear();

  for (unsigned int i = 0; i < nbeads - 3; i++) {
    for (unsigned int j = i + 3; j < nbeads; j++) {
      if (nonlocal_bonds.get_bond_mask(i, j)) {
        double dx = pos[i].x - pos[j].x;
        double dy = pos[i].y - pos[j].y;
        double dz = pos[i].z - pos[j].z;
        box.mindist(dx, dy, dz);

        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        dist.emplace_back(d);
      }
    }
  }
}

bool process_event(const MinNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config, UpdateConfigWriter &) {
  const unsigned int nbeads = sys.pos.size();

  // if the actual (counter[ev.i, j]) and predicted (ni, nj)
  // collision counters of the colliding beads match,
  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  final_vel(sys.pos, sys.vel, p.near_min2, ev.i, ev.j, box, {}, p.m);

  // predict future collisions of beads which are nearest and next-nearest
  // to colliding beads
  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j - 2, ev.j, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  assert(ev.i + 1 == ev.j);

  if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                  ev.i, ev.j, p.near_min2, p.near_max2);

  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.i + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i, ev.i + 2, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  // fill priority queue with collisions of all particle
  // pairs involving a recently-collided particle (executed
  // many times to update the priority queue and remove
  // invalid entries)
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const MaxNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config, UpdateConfigWriter &) {
  const unsigned int nbeads = sys.pos.size();

  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  final_vel(sys.pos, sys.vel, p.near_max2, ev.i, ev.j, box, {}, p.m);

  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j - 2, ev.j, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  assert(ev.i + 1 == ev.j);

  if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                  ev.i, ev.j, p.near_min2, p.near_max2);

  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.i + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i, ev.i + 2, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const MinNNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config, UpdateConfigWriter &) {
  const unsigned int nbeads = sys.pos.size();

  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  final_vel(sys.pos, sys.vel, p.nnear_min2, ev.i, ev.j, box, {}, p.m);

  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  assert(ev.i + 2 == ev.j);

  if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                  ev.i, ev.i + 1, p.near_min2, p.near_max2);
  if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                   ev.i, ev.j, p.nnear_min2, p.nnear_max2);
  if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                  ev.j - 1, ev.j, p.near_min2, p.near_max2);

  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const MaxNNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config, UpdateConfigWriter &) {
  const unsigned int nbeads = sys.pos.size();

  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  final_vel(sys.pos, sys.vel, p.nnear_max2, ev.i, ev.j, box, {}, p.m);

  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  assert(ev.i + 2 == ev.j);

  if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                  ev.i, ev.i + 1, p.near_min2, p.near_max2);
  if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                   ev.i, ev.j, p.nnear_min2, p.nnear_max2);
  if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                  ev.j - 1, ev.j, p.near_min2, p.near_max2);

  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const MinNonlocalInnerEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config, UpdateConfigWriter &) {
  const unsigned int nbeads = sys.pos.size();

  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  final_vel(sys.pos, sys.vel, p.rh2, ev.i, ev.j, box, {}, p.m);

  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  if (ev.i + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i, ev.i + 1, p.near_min2, p.near_max2);
  }
  if (ev.i + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i, ev.i + 2, p.nnear_min2, p.nnear_max2);
  }

  if (ev.j > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j - 2, ev.j, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j - 1, ev.j, p.near_min2, p.near_max2);
  }
  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const MaxNonlocalInnerEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config,
                   UpdateConfigWriter &update_config_writer) {
  const unsigned int nbeads = sys.pos.size();

  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  const Config t_bond_mask = p.transient_bonds.get_bond_mask(ev.i, ev.j);

  // flip bit to form bond
  if (final_vel(sys.pos, sys.vel, p.rc2, ev.i, ev.j, box, -p.eps, p.m)) {
    assert(update_config.non_bonded(t_bond_mask));
    update_config.flip_bond(t_bond_mask);
    assert(update_config.bonded(t_bond_mask));

    update_config.U -= p.eps;

    // store the integer of the configuration and the time
    // of the event
    LOG_DEBUG("bond is formed " << update_config.config);
    update_config_writer.config_int.emplace_back(update_config.config);
    update_config_writer.config_time.emplace_back(ev.t);
  }

  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  if (ev.i + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i, ev.i + 1, p.near_min2, p.near_max2);
  }
  if (ev.i + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i, ev.i + 2, p.nnear_min2, p.nnear_max2);
  }

  if (ev.j > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j - 2, ev.j, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j - 1, ev.j, p.near_min2, p.near_max2);
  }
  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const MaxNonlocalOuterEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config,
                   UpdateConfigWriter &update_config_writer) {
  const unsigned int nbeads = sys.pos.size();

  if (sys.counter[ev.i] != ev.ni || sys.counter[ev.j] != ev.nj)
    return false;

  LOG_DEBUG("processing " << ev);

  sys.counter[ev.i]++;
  sys.counter[ev.j]++;

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);
  update_pos(sys.pos[ev.j], sys.vel[ev.j], sys.times[ev.j], ev.t);
  assert(check_overlap(ev.i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  assert(check_overlap(ev.j, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));

  const Config t_bond_mask = p.transient_bonds.get_bond_mask(ev.i, ev.j);

  std::optional<double> dU;
  // if i, j belong to set of transient bonds,
  if (t_bond_mask && update_config.bonded(t_bond_mask)) {
    dU = p.eps;
  }

  // if the bond is broken,
  if (final_vel(sys.pos, sys.vel, p.rc2, ev.i, ev.j, box, dU, p.m)) {
    assert(update_config.bonded(t_bond_mask));
    update_config.flip_bond(t_bond_mask);
    assert(update_config.non_bonded(t_bond_mask));

    update_config.U += p.eps;

    LOG_DEBUG("bond is broken " << update_config.config);
    update_config_writer.config_int.emplace_back(update_config.config);
    update_config_writer.config_time.emplace_back(ev.t);
  }

  if (ev.i > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i - 2, ev.i, p.nnear_min2, p.nnear_max2);
  }
  if (ev.i > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i - 1, ev.i, p.near_min2, p.near_max2);
  }
  if (ev.i + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.i, ev.i + 1, p.near_min2, p.near_max2);
  }
  if (ev.i + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.i, ev.i + 2, p.nnear_min2, p.nnear_max2);
  }

  if (ev.j > 1) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j - 2, ev.j, p.nnear_min2, p.nnear_max2);
  }
  if (ev.j > 0) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j - 1, ev.j, p.near_min2, p.near_max2);
  }
  if (ev.j + 1 < nbeads) {
    if_nearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                    ev.j, ev.j + 1, p.near_min2, p.near_max2);
  }
  if (ev.j + 2 < nbeads) {
    if_nnearest_bond(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times,
                     ev.j, ev.j + 2, p.nnear_min2, p.nnear_max2);
  }

  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.i,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);
  add_events_for_one_bead(sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter,
                          event_queue, sys.times, cells_bead, ev.j,
                          p.transient_bonds, p.permanent_bonds, update_config,
                          p.max_nbonds);

  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.j,
          cells_bead);

  return true;
}

bool process_event(const BeadCellEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells_bead,
                   UpdateConfig &update_config, UpdateConfigWriter &) {
  if (sys.counter[ev.i] != ev.ni)
    return false;

  LOG_DEBUG("processing " << ev);

  move_to_new_cell(cells_bead, ev.i, ev.ixn, ev.iyn, ev.izn);

  update_pos(sys.pos[ev.i], sys.vel[ev.i], sys.times[ev.i], ev.t);

  add_events_for_bead_after_crossing(
      sys.pos, sys.vel, p.rh2, p.rc2, box, sys.counter, event_queue, sys.times,
      cells_bead, ev.i, ev.wall, p.transient_bonds, p.permanent_bonds,
      update_config, p.max_nbonds);

  // add next cell crossing event, based on currently processed
  // cell crossing event, assuming no collisions
  if_cell(sys.pos, sys.vel, box, sys.counter, event_queue, sys.times, ev.i,
          cells_bead);

  return true;
}
