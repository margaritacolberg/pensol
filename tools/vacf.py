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

    nsteps = data['nsteps']
    total_iter = data['total_iter']
    write_step = data['write_step']
    t_bond = data['transient_bonds'][0]
    n_vel = int(nsteps / write_step)

    name = '{}_{}'.format(args.base_name, traj_num)
    csv_out = name + '_vacf.csv'

    hdf5_in = name + '.h5'
    with h5py.File(hdf5_in, 'r') as f:
        vel = f['unique_config']['vel']

        vacf = np.zeros(n_vel)
        for i in range(total_iter):
            print('iteration', i)
            store_vel = []
            for j in range(i * n_vel, (i + 1) * n_vel):
                store_vel.append(vel[j, t_bond[1]] - vel[j, t_bond[0]])
            vacf += calculate_vacf(store_vel)

    vacf /= total_iter

    np.savetxt(csv_out, np.transpose(vacf), delimiter=',', newline='\n')


def calculate_vacf(vel):
    n_vel = len(vel)
    vacf = np.zeros(n_vel)
    ntpoints = np.zeros(n_vel, dtype=int)

    for i in range(n_vel):
        print('round', i)
        for j in range(n_vel - i):
            vxt = vel[j+i][0] * vel[j][0]
            vyt = vel[j+i][1] * vel[j][1]
            vzt = vel[j+i][2] * vel[j][2]

            vacf[i] += (vxt + vyt + vzt) / 3.0
            ntpoints[i] += 1

        if ntpoints[i] > 0:
            vacf[i] /= ntpoints[i]

    return vacf


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json_in', help='json input file')
    parser.add_argument('base_name', help='base name of hardspheres hdf5 file')
    args = parser.parse_args()

    main(args)
