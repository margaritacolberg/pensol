#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# run_diff_eps.py runs an ensemble of trajectories in parallel for one
# transition between an initial state of a folding protein, and a final state
# which contains one additional bond; the states can be two intermediate
# states; run_diff_eps.py is carried out on a computer cluster only
#
# example of how to run:
# python ../tools/run_diff_eps.py hardspheres_4_0100101100_1100101100.json --eps 3.0 --input_hdf5 ../hybridmc_4_0100101100_1100101100/hybridmc_3_0000101100_0100101100.h5
#
# note that run_diff_eps.py creates a new dir which it enters to generate the
# output files, thus, the path of the json input file, init_json.py, and the
# executable must be from this newly created dir, and not from the dir from
# which the simulation is initiated
#
# note also that the higher the value of non-zero eps, the longer the
# simulation will run and the more accurate the results will be

import argparse
import json
import os
import shutil
import subprocess


def main(args):
    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    traj_num = data['traj_num']

    output_json = new_eps(data, args.json, args.eps)
    file_name = '../../tools/init_json.py'
    exe = '../../release/hardspheres'

    dir_name = os.path.splitext(output_json)[0]
    os.makedirs(dir_name, exist_ok=True)

    new_path = os.path.join(dir_name, output_json)
    shutil.move(output_json, new_path)

    command = ['sbatch', '--account=def-jmschofi', '--time=1-0',
            '--cpus-per-task=1', '--mem-per-cpu=256M',
            '--array=1-{}'.format(traj_num), file_name, output_json, exe]

    if args.input_hdf5 is not None:
        command += ['--input_hdf5', args.input_hdf5]

    if args.dry_run:
        print('would run: {}'.format(command, cwd=dir_name))
    else:
        subprocess.check_call(command, cwd=dir_name)


def new_eps(data, input_file, eps):
    data['eps'] = eps

    input_name = os.path.splitext(input_file)[0]

    print('create', input_name)

    json_name = input_name + '_eps_' + str(eps) + '.json'
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    return json_name


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('--eps', type=float, help='epsilon input value')
    parser.add_argument('--input_hdf5', help='hdf5 input file')
    parser.add_argument('-n', '--dry_run', action='store_true', default=False,
            help='output commands without running simulations')
    args = parser.parse_args()

    main(args)
