#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# run_vacf_or_msd.py calculates the velocity autocorrelation function (vacf) or
# the mean squared displacement (msd) between a pair of beads involved in the
# transient bond; run_vacf_or_msd.py is carried out on a computer cluster only
#
# example of how to run:
# python ../../tools/run_vacf_or_msd.py ../hardspheres_4_0100101100_1100101100.json ../../tools/msd.py hardspheres_4_0100101100_1100101100
#
# note that to run run_vacf_or_msd.py, the user must first enter the transition
# specific directory (for example,
# hardspheres_4_0100101100_1100101100_eps_3.0), to run the above command

import argparse
import json
import subprocess


def main(args):
    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    traj_num = data['traj_num']

    command = ['sbatch', '--account=def-jmschofi', '--time=1-0',
            '--cpus-per-task=1', '--mem-per-cpu=2G',
            '--array=1-{}'.format(traj_num), args.file_name, args.json,
            args.base_name]

    subprocess.check_call(command)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('file_name', help='vacf or msd python file to run')
    parser.add_argument('base_name',
            help='base name of hardspheres hdf5 input file')
    args = parser.parse_args()

    main(args)
