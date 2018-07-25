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
import sys

import numpy as np
import rsf.api as rsf


def main(argv=None):
    if argv is None:
        argv = sys.argv

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('filename1', help='Filename of the 1st file for comparison')
    p.add_argument('filename2', help='Filename of the 2nd file for comparison')
    p.add_argument('--rtol', help='Relative tolerance',
                   type=float, default=1e-5)
    p.add_argument('--atol', help='Absolute tolerance',
                   type=float, default=1e-15)
    p.add_argument('--verbose', '-v', help='Print results to stdout',
                   action='store_true')

    args = p.parse_args()

    filename1, filename2 = args.filename1, args.filename2
    rtol, atol = args.rtol, args.atol
    verbose = args.verbose

    file1 = rsf.Input(filename1)
    file2 = rsf.Input(filename2)

    data1 = file1.getalldata()
    data2 = file2.getalldata()

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


if __name__ == '__main__':
    sys.exit(main())
