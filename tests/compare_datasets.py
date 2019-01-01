#!/usr/bin/env python
"""Compare two RSF files to determine if their datasets are equal.

Comparison is done using expression:

    abs(a - b) <= (atol + rtol * abs(b))

where a and b are datasets, atol is the absolute tolerance (default 1e-8),
rtol the relative tolerance (default 1e-5).

The result of the comparison is provided as an exit code of this script.
Exit codes are:
0 - datasets are close to each other within tolerance
1 - datasets are not close to each other within tolerance
2 - datasets have different dimensions
"""

import argparse
import json
import os
import sys

import numpy as np
import rsf.api as rsf


def main(argv=None):
    if argv is None:
        argv = sys.argv

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('filename1', help='Filename of the 1st file for comparison')
    p.add_argument('filename2', help='Filename of the 2nd file for comparison')
    p.add_argument(
        '--rtol', help='Relative tolerance', type=float, default=1e-5)
    p.add_argument(
        '--atol', help='Absolute tolerance', type=float, default=1e-15)
    p.add_argument(
        '--verbose', '-v', help='Print results to stdout', action='store_true')

    args = p.parse_args()

    filename1, filename2 = args.filename1, args.filename2
    rtol, atol = args.rtol, args.atol
    verbose = args.verbose

    data1, data2 = get_datasets(filename1, filename2)

    try:
        result = np.allclose(data1, data2, rtol=rtol, atol=atol)
    except ValueError:
        if verbose:
            print('Shapes do not match')
        # `np.allclose` returns `ValueError`, when shapes do not match.
        return 2

    if result:
        if verbose:
            print('Datasets are elementwise close within tolerance')
        return 0
    else:
        if verbose:
            print('Datasets are not elementwise close within tolerance')
        return 1


def get_datasets(filename1, filename2):
    rsf_format = '.rsf' in filename1 and '.rsf' in filename2
    bin_format = '.bin' in filename1 and '.bin' in filename2

    if rsf_format:
        file1 = rsf.Input(filename1)
        file2 = rsf.Input(filename2)

        data1 = file1.getalldata()
        data2 = file2.getalldata()
    elif bin_format:
        top_dir = filename1.split('/')[0]
        json_filename = os.path.join(top_dir, 'in_and_out', 'asofi3D.json')
        params = read_asofi3d_json(json_filename)

        TSNAP1, TSNAP2 = params['TSNAP1'], params['TSNAP2']
        TSNAPINC = params['TSNAPINC']
        nsnap = int(1 + np.floor((TSNAP2 - TSNAP1) / TSNAPINC))

        NX, NY, NZ = params['NX'], params['NY'], params['NZ']
        IDX, IDY, IDZ = params['IDX'], params['IDY'], params['IDZ']

        data1 = np.fromfile(filename1, dtype=np.float32)
        data2 = np.fromfile(filename2, dtype=np.float32)
        data1.reshape((NY / IDY, NX / IDX, NZ / IDZ, nsnap))
        data2.reshape((NY / IDY, NX / IDX, NZ / IDZ, nsnap))
    else:
        raise Exception('Cannot determine dataset format!')

    return data1, data2


def read_asofi3d_json(json_filename):
    with open(json_filename) as fp:
        line_list = fp.readlines()

    clean_lines = []
    for line in line_list:
        line_strip = line.strip()
        if line_strip.startswith('#') or line_strip.startswith('//'):
            continue
        else:
            clean_lines.append(line)

    clean_json_string = ''.join(clean_lines)
    params = json.loads(clean_json_string)

    for key in ['NX', 'NY', 'NZ', 'IDX', 'IDY', 'IDZ']:
        params[key] = int(params[key])

    for key in ['TIME', 'TSNAP1', 'TSNAP2', 'TSNAPINC']:
        params[key] = float(params[key])

    if params['TSNAP2'] > params['TIME']:
        params['TSNAP2'] = params['TIME']

    return params


if __name__ == '__main__':
    sys.exit(main())
