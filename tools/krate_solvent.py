#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# krate_solvent.py uses the configurational probability obtained at each step
# of a trajectory, averaged over an ensemble of trajectories, to construct a
# decay curve of configurational probability vs. t; the fit of this curve
# produces the krate, Pb, and the standard error for the transition whose
# initial and final states differ by one bond in the presence of a solvent
#
# example of how to run:
# python ../../tools/krate_solvent.py hardspheres_4_0100101100_1100101100_eps_3.0.json
#
# note that the dir cannot contain any .h5.tmp files; if any exist, they must
# be deleted before running krate_solvent.py
#
# note also that to run krate_solvent.py, the user must first enter the
# transition specific directory (for example,
# hardspheres_4_0100101100_1100101100_eps_3.0), to run the above command

import matrix_element

import argparse
import csv
import h5py
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import scipy.integrate
import scipy.optimize
from operator import add
from sklearn import utils


def main(args):
    with open(args.json, 'r') as input_json:
        data_json = json.load(input_json)

    bits_from_json = re.search(r'_[0-9]+_([01]+)_([01]+)_.*\.json$', args.json)
    bits_i, bits_j = bits_from_json.groups()

    eps = data_json['eps']
    del_t = data_json['del_t']
    traj_num = 0

    nsteps = data_json['nsteps']
    max_t = (nsteps - 1) * del_t

    config_prob = []
    config_prob_store = []
    for file_path in os.listdir('.'):
        if re.search('hardspheres_[0-9]+_([01]+)_([01]+)_[0-9]+\.h5',
                file_path):
            traj_num += 1
            with h5py.File(file_path, 'r') as f:
                config_prob_i = f['unique_config']['config_prob'][:][0]

                config_prob_store.append(config_prob_i)

                if len(config_prob) == 0:
                    config_prob = config_prob_i
                else:
                    config_prob = list(map(add, config_prob, config_prob_i))

    config_prob_store = list(map(list, zip(*config_prob_store)))

    print('number of trajectories:', traj_num)
    config_prob = np.array(config_prob) / traj_num
    config_prob = [1.0 - i for i in config_prob]

    t = []
    for i in range(len(config_prob)):
        t.append(i * del_t)

    t = np.array(t)

    Pb, krate = Pb_krate(eps, config_prob, t)
    nboot = 300
    se = Pb_krate_se(config_prob_store, t, eps, traj_num, nboot)

    print('parameters: Pb: {}, krate: {}'.format(Pb, krate))
    print('standard error: {}'.format(se))

    csv_name = '../krate_solvent.csv'
    with open(csv_name, 'a') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows([[bits_i, bits_j, eps, Pb, krate, se]])

    plt.plot(t, config_prob, label='config prob data')
    plt.plot(t, P(t, Pb, krate), '--', label='fit')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('Probability')
    plt.savefig('probability_{}.pdf'.format(eps), format='pdf')


def P(t, Pb, krate):
    return Pb - ((Pb - 1.0) * np.exp(-krate * t))


def Pb_krate(eps, config_prob, t):
    # Pb is obtained from the plateau of the fit over the configurational
    # probability vs. t decay curve; Pb is defined as p_i / (p_i + p_j), where
    # p_i is the probability of being in state i, the initial state, and p_j is
    # the probability of being in state j, the final state with one bond more
    # than state i
    Pb = 0.0
    # krate is defined as K_ij + K_ji, where K_ij the rate of forming a bond by
    # transitioning from state i to j, and K_ji is the rate of breaking a bond
    # by transitioning from state j to state i
    krate = 0.0
    # for eps = 1.0, data is particularly noisy so curve fitting uses different
    # method to produce closer fit
    if eps == 1.0:
        # number of points averaged over to produce Pb is arbitrarily chosen
        Pb = np.sum(config_prob[len(config_prob)-400:len(config_prob)-1]) / 399
        krate = (1 / scipy.integrate.simps(config_prob - Pb, t)) * (1.0 - Pb)
    else:
        params, cv = scipy.optimize.curve_fit(P, t, config_prob)
        Pb, krate = params

    return Pb, krate


def Pb_krate_se(config_prob_store, t, eps, traj_num, nboot):
    Pb = []
    krate = []
    se = []
    for i in range(nboot):
        boot = []
        for j in range(len(config_prob_store)):
            boot_j = utils.resample(config_prob_store[j],
                    n_samples=len(config_prob_store[j]), random_state=None)
            boot.append(1.0 - np.mean(boot_j))

        Pb_i, krate_i = Pb_krate(eps, boot, t)
        Pb.append(Pb_i)
        krate.append(krate_i)

        residuals = boot - P(t, Pb_i, krate_i)
        ss_res = np.sum(residuals**2)
        n_val = len(boot)
        se.append(np.sqrt(ss_res / (n_val - 2)))

    return np.mean(se)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    args = parser.parse_args()

    main(args)
