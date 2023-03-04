#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# example of how to run:
# python ../tools/run_krate_hybridmc.py 8 16 3.0

import argparse
import glob
import re
import subprocess


def main(args):
    dirs = 'hardspheres_*'

    diff_file = '../../tools/plot_msd.py'
    krate_file = '../../tools/krate_hybridmc.py'

    for dir_name in sorted(glob.glob(dirs)):
        # only dir with eps = 0.0 can be used to calculate diffusion
        # coefficient for hybridmc model
        if not dir_name.endswith('0.0'):
            continue

        json_name = dir_name + '.json'
        rm_eps = re.sub('_eps_0.0$', '', dir_name)

        base_name = 'hybridmc' + rm_eps.lstrip('hardspheres')
        hybridmc_name = base_name + '.csv'
        hybridmc_path = '../{}/{}'.format(base_name, hybridmc_name)
        s_bias_name = base_name + '_s_bias_error.csv'
        s_bias_path = '../{}/{}'.format(base_name, s_bias_name)

        diff_command = ['python', diff_file, json_name, args.t_i, args.t_j]
        subprocess.check_call(diff_command, cwd=dir_name)
        krate_command = ['python', krate_file, json_name, hybridmc_path,
                s_bias_path, args.eps]
        subprocess.check_call(krate_command, cwd=dir_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('t_i',
            help='index of first t in range of t for fitting')
    parser.add_argument('t_j',
            help='index of last t in range of t for fitting')
    parser.add_argument('eps', help='epsilon input value')
    args = parser.parse_args()

    main(args)
