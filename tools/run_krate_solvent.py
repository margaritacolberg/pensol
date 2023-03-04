#!/usr/bin/env python3
# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import glob
import re
import subprocess


def main(args):
    dirs = 'hardspheres_*'

    krate_file = '../../tools/krate_solvent.py'

    for dir_name in sorted(glob.glob(dirs)):
        if re.search(r'hardspheres_[0-9]+_([01]+)_([01]+)_eps_.*' + args.eps,
                dir_name):
            json_name = dir_name + '.json'

            command = ['python', krate_file, json_name]
            subprocess.check_call(command, cwd=dir_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('eps', help='epsilon input value')
    args = parser.parse_args()

    main(args)
