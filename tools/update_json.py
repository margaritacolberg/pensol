#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# update_json.py creates a json input for each solvent simulation based on a
# master json file contained in the examples dir and the corresponding hybridmc
# transition contained in the crambin dir; note that the newly created
# hardspheres json files, and the hybridmc json and hdf5 files are both needed
# for run_diff_eps.py, so the crambin dir can be used for data collection
# without moving the json files elsewhere
#
# example of how to run:
# python ../tools/update_json.py ../examples/crambin.json

import argparse
import glob
import json
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='hardspheres master json input file')
    args = parser.parse_args()

    output = run(args.json)


def run(json_ha):
    src = 'hybridmc_*/hybridmc_*.json'

    for json_hy in glob.glob(src):
        json_name = re.search(r'_([0-9]+)_([01]+)_([01]+)\.json$', json_hy)

        if json_name is None:
            continue

        layer, bits_in, bits_out = json_name.groups()

        output_name = 'hardspheres_{}_{}_{}'.format(layer, bits_in, bits_out)

        with open(json_hy, 'r') as input_json:
            data_hy = json.load(input_json)

        with open(json_ha, 'r') as input_json:
            data_ha = json.load(input_json)

        data_ha['transient_bonds'] = data_hy['transient_bonds']
        data_ha['permanent_bonds'] = data_hy['permanent_bonds']

        json_name = '{}.json'.format(output_name)
        with open(json_name, 'w') as output_json:
            json.dump(data_ha, output_json)


if __name__ == '__main__':
    main()
