#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# calc_diff.py calculates the diffusion coefficient using both old and new sets
# of parameters; results are printed to terminal
#
# example of how to run:
# python calc_diff.py ../examples/crambin.json 8.87 0.013

import argparse
import json
import math
import scipy.special as sc


def main(args):
    with open(args.json, 'r') as input_json:
        json_data = json.load(input_json)

    m = json_data['m']
    m_sol = json_data['m_sol']
    mu = m / m_sol
    V = json_data['length']**3
    rho = json_data['nsol'] / V
    new_rho = 20
    c_gamma = 1/3
    del_t = json_data['del_t_coll']
    new_del_t = 0.018

    # see Schofield et al., J. Chem. Phys., 2012, 136
    #
    # Kummer's function of the first kind
    kummer = sc.hyp1f1(1, 2 + mu, -rho)
    new_kummer = sc.hyp1f1(1, 2 + mu, -new_rho)
    gamma = ((1 - c_gamma) / (1 + mu)) * rho * kummer
    new_gamma = ((1 - c_gamma) / (1 + mu)) * new_rho * new_kummer

    # original (old parameter) D
    D = (1 / mu) * del_t * ((2 - gamma) / (2 * gamma))
    print('old parameter D in simulation units:', D)
    # adjust rho and del_t so that D matches with Stokes-Einstein D, if
    # original D is not identical to Stokes-Einstein D
    new_D = (1 / mu) * new_del_t * ((2 - new_gamma) / (2 * new_gamma))
    print('new parameter D in simulation units:', new_D)

    # length of covalent bond between two amino acids (from Bayat et al., J.
    # Chem. Phys., 2012, 136)
    a = 3.84 * 10**(-10)
    amu = 1.6605 * 10**(-27)
    # mass of one water molecule
    m = 18.0 * amu
    kT = (1.3806 * 10**(-23)) * 310.15
    # see Howard et al., Curr. Opin. Chem. Eng., 2019, 23
    tau = math.sqrt((m * a**2) / kT)
    print('old parameter D in m^2/s:', D * (a**2 / tau))
    print('new parameter D in m^2/s:', new_D * (a**2 / tau))

    # viscosity of water at body temperatures (310.15 K), in Pa s
    nu = 0.0006913
    # radius of generic amino acid
    R = 0.5 * 10**(-9)
    D_se = kT / (6 * math.pi * nu * R)
    print('Stokes-Einstein D in m^2/s:', D_se)
    print('Stokes-Einstein D in simulation units:', D_se * (tau / a**2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('new_rho', help='rho for new parameters')
    parser.add_argument('new_del_t', help='del_t for new parameters')
    args = parser.parse_args()

    main(args)
