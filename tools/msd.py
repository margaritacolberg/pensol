#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import h5py
import json
import numpy as np
import os


def main(args):
    traj_num = os.getenv('SLURM_ARRAY_TASK_ID', default=1)

    with open(args.json_in, 'r') as input_json:
        data = json.load(input_json)

    del_t = data['del_t']
    nsteps = data['nsteps']
    total_iter = data['total_iter']
    write_step = data['write_step']
    t_bond = data['transient_bonds'][0]
    n_pos = int(nsteps / write_step)

    name = '{}_{}'.format(args.base_name, traj_num)
    csv_out = name + '_msd.csv'
    if os.path.exists(csv_out):
        return

    hdf5_in = name + '.h5'
    with h5py.File(hdf5_in, 'r') as f:
        pos = f['pos']

        msd = np.zeros(n_pos)
        for i in range(total_iter):
            print('iteration', i)
            store_pos = []
            for j in range(i * n_pos, (i + 1) * n_pos):
                store_pos.append(pos[j, t_bond[1]] - pos[j, t_bond[0]])
            msd += calculate_msd(store_pos, del_t)

    msd /= total_iter

    np.savetxt(csv_out, np.transpose(msd), delimiter=',', newline='\n')


def calculate_msd(pos, del_t):
    n_pos = len(pos)
    msd = np.zeros(n_pos)
    ntpoints = np.zeros(n_pos, dtype=int)

    for i in range(1, n_pos):
        print('round', i)
        for j in range(n_pos - i):
            dxt = pos[j+i][0] - pos[j][0]
            dyt = pos[j+i][1] - pos[j][1]
            dzt = pos[j+i][2] - pos[j][2]

            msd[i] += (dxt**2 + dyt**2 + dzt**2)
            ntpoints[i] += 1

        if ntpoints[i] > 0:
            msd[i] /= ntpoints[i]

    return msd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json_in', help='json input file')
    parser.add_argument('base_name', help='base name of hardspheres hdf5 file')
    args = parser.parse_args()

    main(args)
