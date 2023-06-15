// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// solvent.cc simulates the penetrating solvent model, in which the solvent
// particles are simulated implicitly by specifying only the random forces
// that occur on the surface of the protein to mimic solvent-protein collisions

#include "solvent.h"

// update velocities of the beads after collisions with solvent particles
void collide_sol(System &sys, Random &mt, const Param &p) {
  const double Vc = 1.0;
  const double rho = p.nsol / (p.length * p.length * p.length);

  const unsigned int nbeads = sys.pos.size();

  const double a = p.angle * M_PI / 180;
  unsigned int cos_a = cos(a);
  unsigned int sin_a = sin(a);

  for (unsigned int i = 0; i < nbeads; i++) {
    double u_x, u_y, u_z;
    unit_sphere(mt, u_x, u_y, u_z);

    // rotation matrix with angle "a" degrees
    const double Rxx = cos_a + u_x * u_x * (1 - cos_a);
    const double Rxy = u_x * u_y * (1 - cos_a) - u_z * sin_a;
    const double Rxz = u_x * u_z * (1 - cos_a) + u_y * sin_a;
    const double Ryx = u_y * u_x * (1 - cos_a) + u_z * sin_a;
    const double Ryy = cos_a + u_y * u_y * (1 - cos_a);
    const double Ryz = u_y * u_z * (1 - cos_a) - u_x * sin_a;
    const double Rzx = u_z * u_x * (1 - cos_a) - u_y * sin_a;
    const double Rzy = u_z * u_y * (1 - cos_a) + u_x * sin_a;
    const double Rzz = cos_a + u_z * u_z * (1 - cos_a);

    // number of solvent particles in a cell
    std::poisson_distribution<unsigned int> poisson(Vc * rho);
    unsigned int Ns;

    do {
      Ns = poisson(mt);
    } while (Ns == 0);

    // total mass of solvent particles in a cell
    double mass = Ns * p.m_sol;
    // total mass of all particles
    double M = p.m + mass;

    const double mean = 0.0;
    double sigma = std::sqrt(p.temp / mass);
    std::normal_distribution<> gaussian{mean, sigma};

    double v1 = gaussian(mt);
    double v2 = gaussian(mt);
    double v3 = gaussian(mt);

    // calculate the center of mass velocity
    const double v_cm_x = ((p.m / M) * sys.vel[i].x) + ((mass / M) * v1);
    const double v_cm_y = ((p.m / M) * sys.vel[i].y) + ((mass / M) * v2);
    const double v_cm_z = ((p.m / M) * sys.vel[i].z) + ((mass / M) * v3);

    double dv_x = sys.vel[i].x - v_cm_x;
    double dv_y = sys.vel[i].y - v_cm_y;
    double dv_z = sys.vel[i].z - v_cm_z;

    sys.vel[i].x = v_cm_x + dv_x * Rxx + dv_y * Rxy + dv_z * Rxz;
    sys.vel[i].y = v_cm_y + dv_x * Ryx + dv_y * Ryy + dv_z * Ryz;
    sys.vel[i].z = v_cm_z + dv_x * Rzx + dv_y * Rzy + dv_z * Rzz;
  }
}
