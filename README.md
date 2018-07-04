# Time Domain Group

Extension of the code SOFI3D to anisotropic wave propagation (ASOFI3D).

## Obtaining the code

Clone the code:

    git clone git@github.com:swag-kaust/TD.git

Then change your current directory:

    cd TD

and finally switch to the `refactoring` branch:

    git checkout -b refactoring


## Building the code

You need a C compiler and and an MPI library to compile the code.
Known working configuration (tested on `kw13806` workstation) is Intel C
compiler and the Intel MPI library.
To load them on your KAUST-provided machine, execute the following commands:

    module load intel/2016
    module load intelmpi/2016/intel-2016

Then, to compile the code, change the directory:

    cd ASOFI3D/src

and then:

    make all

Currently, compilation emits multiple warnings that can be safely ignored.

After successful compilation, you can run the code, by going back to the root
directory:

    cd ../..

and run:

    ./run_ASOFI3D.sh np

where `np` is a number of MPI processes you want.


#OLD INSTRUCTIONS
in order to compile the code on a KAUST workstation you need to initialize mpi
e.g. run source mpiInit.sh


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

