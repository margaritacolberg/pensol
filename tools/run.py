#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# run.py runs an ensemble of trajectories in series for one transition between
# an initial state of a folding protein, and a final state which contains one
# additional bond; the states can be two intermediate states; run.py is carried
# out on a local computer only; recommended for testing and debugging purposes
# only, using one or a small number of short trajectories, since the simulation
# is resource intensive and is best reserved for computer clusters
#
# example of how to run:
# python ../tools/run.py ../hardspheres_4_0100101100_1100101100.json --input_hdf5 ../hybridmc_4_0100101100_1100101100/hybridmc_3_0000101100_0100101100.h5
#
# note that run.py creates a new dir which it enters to generate the output
# files, thus, the path of the json input file, init_json.py and the executable
# must be from this newly created dir, and not from the dir from which the
# simulation is initiated

import argparse
import os
import subprocess


def main(args):
    file_name = os.path.basename(args.json)
    dir_name = os.path.splitext(file_name)[0]

    if os.path.isdir(dir_name):
        return

    tmp_dir_name = '{}.tmp'.format(dir_name)

    file_name = '../../tools/init_json.py'
    exe = '../../release/hardspheres'

    init_json_args_list = ['python', file_name, args.json, exe]

    if args.input_hdf5 is not None:
        init_json_args_list += ['--input_hdf5', args.input_hdf5]

    if not os.path.isdir(tmp_dir_name):
        os.mkdir(tmp_dir_name)

    subprocess.run(init_json_args_list, cwd=tmp_dir_name)

    os.rename(src=tmp_dir_name, dst=dir_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('--input_hdf5', help='hdf5 input file')
    args = parser.parse_args()

    main(args)
