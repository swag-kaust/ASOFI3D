## About

CircleCI status:

[![CircleCI](https://circleci.com/gh/swag-kaust/TD.svg?style=svg&circle-token=2bb57e94a999ba7f33afb12bf091751af6bda219)](https://circleci.com/gh/swag-kaust/TD)

![](doc/tex/eps/invisible_gamma_1.gif)

ASOFI3D stands for Anisotropic Seismic mOdeling with FInite differences in 3D.
This code is a modification of
[SOFI3D](https://git.scc.kit.edu/GPIAG-Software/SOFI3D/wikis/home) code
to accommodate orthorhombic anisotropy.


## Obtaining the code

Get the code from this repo:

    git clone git@github.com:swag-kaust/ASOFI3D.git

Then switch to the created directory:

    cd ASOFI3D


## Building the code

The prerequisites for ASOFI3D are:

* C compiler (for example, `gcc`) with the support of C11 standard
* MPI library (for example, [OpenMPI](https://www.open-mpi.org/))
* [GNU Make](https://www.gnu.org/software/make/) build system


## Compile with gcc and OpenMPI on Ubuntu

On recent Ubuntu versions such as 14.04, 16.04, or 18.04 all prerequisites
can be obtained by the following commands:

    sudo apt-get install gcc make libopenmpi-dev

Then while in the root directory of the code, build the code via command

    make

which compiles the solver and several auxiliary utilities.

## Example usage

After successful compilation, you can run the code via command

    ./run_asofi3D.sh np dirname

where `np` is a number of MPI processes you want to use and `dirname` is the
directory that contains configuration of the problem to solve.
Parameter `dirname` is optional and defaults to `par`, so that the main
configuration file of the solver is `par/in_and_out/asofi3D.json`.


## Running the tests

To run the tests, Madagascar is an additional prerequisite.
Tests are run via the command

    make test
