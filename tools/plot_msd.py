#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# plot_msd.py calculates the diffusion coefficient from the mean squared
# displacement (msd) and plots the msd vs. t and the linear fit over a portion
# of the msd vs. t data from which the diffusion coefficient is obtained
#
# example of how to run:
# python ../../tools/plot_msd.py hardspheres_4_0100101100_1100101100_eps_0.0.json 40 60 10
#
# note that for crambin, plot_msd.py is run only in dir for which eps was set
# to 0.0
#
# note also that to run plot_msd.py, the user must first enter the transition
# specific directory (for example,
# hardspheres_4_0100101100_1100101100_eps_0.0), to run the above command
#
# in the example above, the values of t_i, t_j, and the gap are set for the old
# parameter runs; the new parameter runs require the values 150, 300, and 50

import argparse
import csv
import json
import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import re
import scipy.stats
from sklearn import utils


matplotlib.rcParams.update({
    # figsize: width, height
    "figure.figsize": [2.75, 2.5],
    "figure.subplot.left": 0.18,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.99,
    "pgf.texsystem": "lualatex",
    "pgf.rcfonts": False,
    "pgf.preamble": [
        "\\usepackage{amsmath}",
        "\\usepackage{unicode-math}",
        "\\setmainfont{TeX Gyre Pagella}",
        "\\setmathfont{TeX Gyre Pagella Math}",
        "\\setmathfont[range={cal,bfcal},StylisticSet=1]{XITS Math}"
        ],
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 11,
    "legend.labelspacing": 0.2,
    "legend.loc": "lower right",
})


def main(args):
    with open(args.json_in, 'r') as input_json:
        data = json.load(input_json)

    del_t = data['del_t']
    nsteps = data['nsteps']
    write_step = data['write_step']
    n_pos = int(nsteps / write_step)

    t = []
    for i in range(n_pos):
        t.append(i * float(del_t))

    t_i = int(args.t_i)
    t_j = int(args.t_j)
    gap = int(args.gap)

    traj_num = 0
    msd = np.zeros(n_pos)
    msd_store = []
    for file_path in os.listdir('.'):
        if re.search('hardspheres_[0-9]+_([01]+)_([01]+)_[0-9]+_msd\.csv',
                file_path):
            print(file_path)
            traj_num += 1
            with open(file_path, 'r') as input_csv:
                data = csv.reader(input_csv, delimiter=',')

                count = 0
                msd_i = []
                for row in data:
                    msd[count] += float(row[0])
                    msd_i.append(float(row[0]))
                    count += 1

                msd_store.append(msd_i)

    msd_store = np.array(msd_store)

    print('number of trajectories:', traj_num)
    msd /= traj_num

    fit = np.polyfit(t[t_i:t_j], msd[t_i:t_j], 1)
    D = fit[0] / 6.0

    nboot = 300
    l_D_err, u_D_err = D_err(nboot, traj_num, msd_store, t, t_i, t_j, gap)

    print('first t is', t[t_i])
    print('last t is', t[t_j])
    print('the diffusion coefficient is', D)
    print('the CI of the diffusion coefficient is [{}, {}]'.format(l_D_err, u_D_err))

    with open('diff_coeff.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerow([D, l_D_err, u_D_err])

    line = []
    t_fit = []
    for i in range(len(t)):
        if t[i] > t[t_i] and t[i] < t[t_j]:
            line_i = (fit[0] * t[i]) + fit[1]
            line.append(line_i)
            t_fit.append(t[i])

    plt.plot(t[30:70], msd[30:70], label='Data')
    plt.plot(t_fit, line, label='Fit')
    plt.legend()
    plt.xlabel('$t$')
    plt.ylabel('MSD')
    plt.savefig('msd.pdf', format='pdf')
    plt.clf()


def D_err(nboot, traj_num, msd, t, t_min, t_max, gap):
    D = []
    for i in range(nboot):
        t_i, t_j = get_times(t_min, t_max, gap)

        boot_i = utils.resample(msd, n_samples=len(msd),
                random_state=None)

        boot_i = np.sum(boot_i, axis=0) / traj_num
        fit = np.polyfit(t[t_i:t_j], boot_i[t_i:t_j], 1)
        D.append(fit[0] / 6.0)

    D = np.sort(D)

    lower_ind = int(0.025 * nboot)
    upper_ind = int(0.975 * nboot)

    lower_D = D[lower_ind]
    upper_D = D[upper_ind]

    return lower_D, upper_D


def get_times(t_min, t_max, gap):
    t_i = random.randint(t_min, t_max-gap)
    t_j = random.randint(t_i+gap, t_max)

    return t_i, t_j


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json_in', help='json input file')
    parser.add_argument('t_i',
            help='index of first t in range of t for fitting')
    parser.add_argument('t_j',
            help='index of last t in range of t for fitting')
    parser.add_argument('gap',
            help='minimum gap in indices between first t and last t')
    args = parser.parse_args()

    main(args)
