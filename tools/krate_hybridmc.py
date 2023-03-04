#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# krate_hybridmc.py uses the diffusion coefficient obtained from the solvent
# model to calculate the krate and Pb for the hybridmc model for the transition
# whose initial and final states differ by one bond
#
# example of how to run:
# python ../../tools/krate_hybridmc.py hardspheres_4_0100101100_1100101100_eps_0.0.json ../hybridmc_4_0100101100_1100101100/hybridmc_4_0100101100_1100101100.csv ../hybridmc_4_0100101100_1100101100/hybridmc_4_0100101100_1100101100_s_bias_error.csv 3.0
#
# note that the value of eps specified on the command line must match the eps
# chosen for the solvent krate calculation when comparing hybridmc and solvent
# models
#
# note also that to run krate_hybridmc.py, the user must first enter the
# transition specific directory (for example,
# hardspheres_4_0100101100_1100101100_eps_0.0), to run the above command

import matrix_element

import argparse
import csv
import json
import math
import numpy as np


def main(args):
    with open(args.json, 'r') as input_json:
        json_data = json.load(input_json)

    beta = 1 / json_data['temp']
    eps = float(args.eps)
    del_t = json_data['del_t']

    s_bias_diff = np.array(get_s_bias(args.csv_s))
    # the final state j has one more bond than the initial state i
    bits_i, bits_j, t_ind, int_i, int_j = matrix_element.get_state_data(args.csv_t)

    P_i_div_P_j = math.exp((-beta*eps) + s_bias_diff)
    Pb = 1 / (1 + (1 / P_i_div_P_j))

    inner_fpt, outer_fpt = matrix_element.get_fpt(args.csv_t, t_ind)

    D = get_D('diff_coeff.csv')

    K_ji_inv = P_i_div_P_j * (inner_fpt / D) + (outer_fpt / D)
    # K_ij the rate of forming a bond by transitioning from state i to j, and
    # K_ji is the rate of breaking a bond by transitioning from state j to
    # state i
    K_ji = 1 / K_ji_inv
    K_ij = K_ji * P_i_div_P_j

    krate = K_ij + K_ji
    print('rate from hybridmc is', krate)

    csv_name = '../krate_hybridmc.csv'
    with open(csv_name, 'a') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows([[bits_i, bits_j, eps, D, Pb, krate]])


def get_D(csv_in):
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            D = float(row[0])

    return D


def get_s_bias(csv_in):
    s_bias = []
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            s_bias.append(float(row[0]))

    return s_bias


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('csv_t', help='csv mfpt file')
    parser.add_argument('csv_s', help='csv s_bias file')
    parser.add_argument('eps', help='epsilon input value')
    args = parser.parse_args()

    main(args)
