#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# plot_msd.py calculates the diffusion coefficient from the mean squared
# displacement (msd) and plots the msd vs. t and the linear fit over a portion
# of the msd vs. t data from which the diffusion coefficient is obtained
#
# example of how to run:
# python ../../tools/plot_msd.py hardspheres_4_0100101100_1100101100_eps_0.0.json 8 16
#
# note that for crambin, plot_msd.py is run only in dir for which eps was set
# to 0.0
#
# note also that to run plot_msd.py, the user must first enter the transition
# specific directory (for example,
# hardspheres_4_0100101100_1100101100_eps_0.0), to run the above command

import argparse
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import scipy.stats


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

    traj_num = 0
    msd = np.zeros(n_pos)
    boot = []
    for file_path in os.listdir('.'):
        if re.search('hardspheres_[0-9]+_([01]+)_([01]+)_[0-9]+_msd\.csv',
                file_path):
            print(file_path)
            traj_num += 1
            boot_i = []
            with open(file_path, 'r') as input_csv:
                data = csv.reader(input_csv, delimiter=',')

                count = 0
                for row in data:
                    msd[count] += float(row[0])
                    boot_i.append(float(row[0]))
                    count += 1

            fit_i = np.polyfit(t[t_i:t_j], boot_i[t_i:t_j], 1)
            boot.append(fit_i[0] / 6.0)

    sem = scipy.stats.sem(boot, axis=0)
    print('number of trajectories', traj_num)
    msd /= traj_num

    fit = np.polyfit(t[t_i:t_j], msd[t_i:t_j], 1)
    D = fit[0] / 6.0

    print('first t is', t[t_i])
    print('last t is', t[t_j])
    print('the diffusion coefficient is', D)
    print('the error of the diffusion coefficient is', sem)

    with open('diff_coeff.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerow([D])

    line = []
    t_fit = []
    for i in range(len(t)):
        if t[i] > t[t_i] and t[i] < t[t_j]:
            line_i = (fit[0] * t[i]) + fit[1]
            line.append(line_i)
            t_fit.append(t[i])

    plt.plot(t[0:150], msd[0:150], label='MSD data')
    plt.plot(t_fit, line, label='fit')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('MSD')
    plt.savefig('msd.png')
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json_in', help='json input file')
    parser.add_argument('t_i',
            help='index of first t in range of t for fitting')
    parser.add_argument('t_j',
            help='index of last t in range of t for fitting')
    args = parser.parse_args()

    main(args)
