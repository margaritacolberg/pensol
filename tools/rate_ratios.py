#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# rate_ratios.py calculates the ratio of two krates of the hybridmc and solvent
# models for two different eps; this program assumes the krates are obtained
# only for one transition at a range of eps; rate ratios eliminate the need to
# evaluate the diffusion coefficient for hybridmc krates
#
# example of how to run:
# python ../tools/rate_ratios.py eps_Pu_krate.csv avg_s_bias.csv hybridmc_0_0_1.csv

import matrix_element

import argparse
import csv
import math


def main(args):
    _, _, t_ind, int_i, int_j = matrix_element.get_state_data(args.t_csv)

    beta = 1.0

    s_bias = get_s_bias(args.s_csv, int_i, int_j)
    eps, Pu, krate = get_eps_Pu_krate(args.csv_in)

    K_ji = []
    K_ij = []
    for i in range(len(eps)):
        K_elem_ji, K_elem_ij = matrix_element.K_elem(args.t_csv, beta, eps[i],
                0.0, s_bias[0], s_bias[1], t_ind)
        K_ji.append(K_elem_ji)
        K_ij.append(K_elem_ij)

    output_1 = []
    for i in range(len(krate)):
        for j in range(len(krate)):
            if i != j:
                # hybridmc rate ratio
                ratio_1 = (K_ji[i] + K_ij[i]) / (K_ji[j] + K_ij[j])
                # solvent rate ratio
                ratio_2 = krate[i] / krate[j]
                output_1.append([eps[i], eps[j], ratio_1, ratio_2])

    with open('rate_ratios.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output_1)

    output_2 = []
    for i in range(len(Pu)):
        P = P_ratio(beta, eps[i], 0.0, max(s_bias), min(s_bias))
        output_2.append([eps[i], P, Pu[i]])

    with open('prob_compare.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output_2)


def get_s_bias(csv_in, int_i, int_j):
    s_bias = []
    with open(csv_in, 'r') as input_csv:
        s_data = csv.reader(input_csv, delimiter=',')

        for row in s_data:
            state = int(row[0], base=2)
            if (state == int_i) or (state == int_j):
                s_bias.append(float(row[1]))

    return s_bias


def get_eps_Pu_krate(csv_in):
    eps = []
    Pu = []
    krate = []
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            eps.append(float(row[0]))
            Pu.append(float(row[1]))
            krate.append(float(row[2]))

    return eps, Pu, krate


def P_ratio(beta, e_i, e_j, s_i, s_j):
    # the final state j has one more bond than the initial state i
    s_bias_diff = s_i - s_j
    eps = e_i - e_j

    P_i_div_P_j = math.exp((-beta*eps) + s_bias_diff)

    return 1 / (1 + (1 / P_i_div_P_j))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_in', help='eps, Pu, krate csv input file')
    parser.add_argument('s_csv', help='hybridmc avg_s_bias csv input file')
    parser.add_argument('t_csv', help='hybridmc mfpt csv input file')
    args = parser.parse_args()

    main(args)
