#!/usr/bin/env python
# Run dot-product test on the generated data.
# Dot-product test checks self-adjointness of the code.
# Precisely, we run two simulations in which source term is generated
# randomly on time step 1 and then the pressure wavefield is recorded on
# time step 2.
# Then we check that the dot products:
# s_1 .dot. p_2   and   s_2 .dot. p_1
# are equal to each other, where s, p are source and pressure, respectively,
# while subscripts 1 and 2 are the number of simulation.
import argparse
import os
import sys

import numpy as np


def main(argv=None):
    if argv is None:
        argv = sys.argv

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', help='Directory')
    p.add_argument('--rtol', help='Relative tolerance',
                   type=float, default=1e-8)
    p.add_argument('--atol', help='Absolute tolerance',
                   type=float, default=1e-15)
    p.add_argument('--verbose', '-v', help='Print results to stdout',
                   action='store_true')

    args = p.parse_args()

    dirname = args.directory
    rtol = args.rtol
    atol = args.atol

    p_data_1 = _get_data(dirname, 'snap_1/test.bin.p')
    # p_data_1 = np.reshape(p_data_1, (2, -1))

    s_data_1 = _get_data(dirname, 'source_field_1/1.bin')

    p_data_2 = _get_data(dirname, 'snap_2/test.bin.p')
    # p_data_2 = np.reshape(p_data_2, (2, -1))

    s_data_2 = _get_data(dirname, 'source_field_2/1.bin')

    result_1 = np.dot(s_data_1, p_data_2)
    result_2 = np.dot(s_data_2, p_data_1)

    if np.allclose(result_1, result_2, rtol=rtol, atol=atol):
        return 0
    else:
        return 1


def _get_data(dirname, filename):
    full_path = os.path.join(dirname, filename)
    return np.fromfile(full_path, dtype=np.float32)


if __name__ == '__main__':
    sys.exit(main())
