#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# plot_vacf.py calculates a set of diffusion coefficients from the velocity
# autocorrelation function (vacf) and plots the vacf vs. t and diffusion
# coefficient vs. t data; plot_vacf.py is best used for a system of a single
# Brownian particle, as the diffusion coefficient it calculates is unreliable
# for larger systems (in such cases, use plot_msd.py)
#
# example of how to run:
# python ../../tools/plot_vacf.py hardspheres_4_0100101100_1100101100_eps_0.0.json 200 250
#
# note that for crambin, plot_vacf.py is run only in dir for which eps was set
# to 0.0
#
# note also that to run plot_vacf.py, the user must first enter the transition
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


def main(args):
    with open(args.json_in, 'r') as input_json:
        data = json.load(input_json)

    del_t = data['del_t']
    traj_num = data['traj_num']
    nsteps = data['nsteps']
    write_step = data['write_step']
    n_vel = int(nsteps / write_step)

    t = []
    for i in range(n_vel):
        t.append(i * float(del_t))

    t_i = int(args.t_i)
    t_j = int(args.t_j)

    vacf = np.zeros(n_vel)
    for file_path in os.listdir('.'):
        if re.search('hardspheres_[0-9]+_([01]+)_([01]+)_[0-9]+_vacf\.csv',
                file_path):
            print(file_path)
            with open(file_path, 'r') as input_csv:
                data = csv.reader(input_csv, delimiter=',')

                count = 0
                for row in data:
                    vacf[count] += float(row[0])
                    count += 1

    vacf /= traj_num
    print('initial val of vacf is {}'.format(vacf[0]))

    # exclude oscillating data at high values of t
    new_vacf = vacf[:800]
    new_t = t[:800]

    plt.plot(t[:100], vacf[:100])
    plt.xlabel('t')
    plt.ylabel('VACF')
    plt.savefig('vacf.pdf', format='pdf')
    plt.clf()

    diff, time = calculate_diff(new_vacf, new_t)

    D = np.mean(diff[t_i:t_j])
    print('the diffusion coefficient is', D)

    with open('diff_coeff.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerow([D])

    plt.plot(time, diff)
    plt.xlabel('t')
    plt.ylabel('Diffusion coefficient')
    plt.savefig('diff.png')
    plt.clf()


def calculate_diff(vacf, t):
    diff = []
    time = []

    # composite trapezoidal rule for integration
    for i in range(1, len(vacf)):
        trap_sum = vacf[0] + vacf[i]

        for j in range(1, i-1, 1):
            trap_sum += (2.0 * vacf[j])

        trap_sum *= (t[i] / (2.0 * i))

        diff.append(trap_sum.tolist())
        time.append(t[i])

    return diff, time


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json_in', help='json input file')
    parser.add_argument('t_i',
            help='index of first t in range of t for fitting')
    parser.add_argument('t_j',
            help='index of last t in range of t for fitting')
    args = parser.parse_args()

    main(args)
