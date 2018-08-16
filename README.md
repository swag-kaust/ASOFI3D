# About

[![CircleCI](https://circleci.com/gh/swag-kaust/TD/tree/master.svg?style=svg)](https://circleci.com/gh/swag-kaust/TD/tree/master)

ASOFI3D stands for Anisotropic Seismic mOdeling with FInite differences.
This code is an extension of 
[SOFI3D](https://git.scc.kit.edu/GPIAG-Software/SOFI3D/wikis/home)

# Obtaining the code

Get the code from this repo:

    git clone git@github.com:swag-kaust/TD.git

Then switch to the cloned repo directory:

    cd TD

# Building the code

You need a C compiler and an MPI library to compile the code.
Popular combination is C compiler `gcc` and MPI library
[OpenMPI](https://www.open-mpi.org/).
On recent Ubuntu versions such as 14.04, 16.04, or 18.04 all prerequisites
can be obtained by the following commands:

    sudo apt-get install build-essential libopenmpi-dev

Then while in the root directory of the code, command

    make

will compile all codes (solver and auxiliary utilities).
Currently, compilation emits multiple warnings that can be safely ignored.

# Running the code

After successful compilation, you can run the code

    ./run_ASOFI3D.sh np dirname

where `np` is a number of MPI processes you want to use and `dirname` is the
directory that contains configuration of the problem to solve.
Parameter `dirname` is optional and defaults to `par`, so that the main
configuration file of the solver is `par/in_and_out/sofi3D.json`.

## Building code using Intel compiler and MPI library (if you're in SWAG)

You need a C compiler and an MPI library to compile the code.
Known working configuration (tested on `kw13806` workstation) is Intel C
compiler and the Intel MPI library.
To load them on your KAUST-provided machine, execute the following commands:

    module load intel/2016
    module load intelmpi/2016/intel-2016


# SOFI3D instructions

## What is SOFI3D?

SOFI3D stands for Seismic mOdeling with FInite differences and denotes our 3-D
viscoelastic time domain massive parallel modeling code.    The manual and a
reference paper is included in the download archive or can be downloaded
[here](https://git.scc.kit.edu/GPIAG-Software/SOFI3D/wikis/home)

## Download and Newsletter

You can Download the [latest stable Release]
(https://git.scc.kit.edu/GPIAG-Software/SOFI3D/tree/Release)
or the current [beta-version]
(https://git.scc.kit.edu/GPIAG-Software/SOFI3D/tree/master).
Also a SOFI3D branch with additional benchmarks is available:
[overnightbuilt]
(https://git.scc.kit.edu/GPIAG-Software/SOFI3D/tree/overnightbuilt)

To receive news and updates please
[register](https://www.gpi.kit.edu/Software-WS.php) on the email list
sofi@lists.kit.edu.  Please use this list also to ask questions on using the
software or to report problems or bugs.

