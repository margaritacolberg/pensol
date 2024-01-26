#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# krate_hybridmc.py uses the diffusion coefficient obtained from the solvent
# model to calculate the krate and Pu for the hybridmc model for the transition
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

    exp_s_bias_diff, rel_e_s_bias = get_exp_s_bias(args.csv_s)
    # the final state j has one more bond than the initial state i
    bits_i, bits_j, t_ind, int_i, int_j = matrix_element.get_state_data(args.csv_t)

    # percent relative error for P_i_div_P_j is the same as rel_e_s_bias
    P_i_div_P_j = math.exp(-beta * eps) * exp_s_bias_diff
    P_j_div_P_i = 1 / P_i_div_P_j
    Pu = 1 / (1 + P_j_div_P_i)
    Pu_rel_e = (rel_e_s_bias * P_j_div_P_i) / (1 + P_j_div_P_i)
    Pu_err = (Pu_rel_e / 100) * Pu
    print('Pu for hybridmc is {}, with error {}'.format(Pu, Pu_err))

    inner_fpt, outer_fpt = matrix_element.get_fpt(args.csv_t, t_ind)

    D, l_D_rel_e, u_D_rel_e = get_D('diff_coeff.csv')

    term_1 = P_i_div_P_j * (inner_fpt / D)
    term_2 = outer_fpt / D
    K_ji_inv = term_1 + term_2

    # K_ij the rate of forming a bond by transitioning from state i to j, and
    # K_ji is the rate of breaking a bond by transitioning from state j to
    # state i
    K_ji = 1 / K_ji_inv
    K_ij = K_ji * P_i_div_P_j

    l_krate_err = calculate_K_err(l_D_rel_e, term_1, term_2, rel_e_s_bias, K_ji_inv, K_ji, K_ij)
    u_krate_err = calculate_K_err(u_D_rel_e, term_1, term_2, rel_e_s_bias, K_ji_inv, K_ji, K_ij)

    krate = K_ij + K_ji
    print('krate for hybridmc is {}, with error [{}, {}]'.format(krate, l_krate_err, u_krate_err))

    csv_name = '../krate_hybridmc.csv'
    with open(csv_name, 'a') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows([[bits_i, bits_j, eps, D, Pu, Pu_err, krate, l_krate_err, u_krate_err]])


def get_D(csv_in):
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            D = float(row[0])
            l_D_abs_e = D - float(row[1])
            l_D_rel_e = (l_D_abs_e / D) * 100
            u_D_abs_e = float(row[2]) - D
            u_D_rel_e = (u_D_abs_e / D) * 100

    return D, l_D_rel_e, u_D_rel_e


def get_exp_s_bias(csv_in):
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            exp_s_bias = float(row[0])
            rel_e_s_bias = float(row[2])

    return exp_s_bias, rel_e_s_bias


def calculate_K_err(D_rel_e, term_1, term_2, rel_e_s_bias, K_ji_inv, K_ji, K_ij):
    term_1_rel_e = math.sqrt(rel_e_s_bias**2 + 1.0**2 + D_rel_e**2)
    term_1_abs_e = (term_1_rel_e / 100) * term_1
    term_2_rel_e = math.sqrt(5.0**2 + D_rel_e**2)
    term_2_abs_e = (term_2_rel_e / 100) * term_2
    K_ji_inv_abs_e = math.sqrt(term_1_abs_e**2 + term_2_abs_e**2)

    K_ji_rel_e = (K_ji_inv_abs_e / K_ji_inv) * 100
    K_ji_abs_e = (K_ji_rel_e / 100) * K_ji
    K_ij_rel_e = math.sqrt(K_ji_rel_e**2 + rel_e_s_bias**2)
    K_ij_abs_e = (K_ij_rel_e / 100) * K_ij

    krate_err = math.sqrt(K_ji_abs_e**2 + K_ij_abs_e**2)

    return krate_err


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('csv_t', help='csv mfpt file')
    parser.add_argument('csv_s', help='csv s_bias file')
    parser.add_argument('eps', help='epsilon input value')
    args = parser.parse_args()

    main(args)
