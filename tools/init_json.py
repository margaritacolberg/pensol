#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import json
import os
import subprocess
import sys


def main(args):
    traj_num = os.getenv('SLURM_ARRAY_TASK_ID', default=1)

    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    nonlocal_bonds = sorted_bonds(data['nonlocal_bonds'])
    transient_bonds = sorted_bonds(data['transient_bonds'])
    permanent_bonds = sorted_bonds(data['permanent_bonds'])

    bits_in = ''
    bits_out = ''
    layer = 0

    for i in nonlocal_bonds:
        if i in permanent_bonds:
            bits_in += '1'
            bits_out += '1'
            layer += 1
        elif i in transient_bonds:
            bits_in += '0'
            bits_out += '1'
        else:
            bits_in += '0'
            bits_out += '0'

    output_name = 'hardspheres_{}_{}_{}_{}'.format(layer, bits_in, bits_out,
            traj_num)

    hdf5_name = '{}.h5'.format(output_name)

    data['config_in'] = int(bits_in, 2)
    data['config_out'] = int(bits_out, 2)
    data['seeds'] = [hash(os.urandom(32)), hash(traj_num)]

    if os.path.exists(hdf5_name):
        return

    json_name = '{}.json'.format(output_name)
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # for layer = 1 or greater,
    command = [args.hardspheres, json_name, hdf5_name]
    if args.input_hdf5 is not None:
        command += ['--input-file', args.input_hdf5]

    print(command)
    sys.stdout.flush()

    log_name = '{}.log'.format(output_name)
    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                stderr=subprocess.STDOUT)


def sorted_bonds(bond_list):
    bond_list = [sorted(el) for el in bond_list]
    bond_list.sort()
    return bond_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('hardspheres', help='hardspheres executable')
    parser.add_argument('--input_hdf5',
            help='hdf5 with initial state with 1+ permanent bond')
    args = parser.parse_args()

    main(args)
