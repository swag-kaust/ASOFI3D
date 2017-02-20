#!/bin/bash

# Just compile + run + prepare snapshots for visualization

# Load MPI modules in case if you are using KAUST workstation
source mpiInit.sh

# Go to directory with stock-given scripts
cd ASOFI3D/par

# Compile the whole code
echo "Compilation"
./compileSOFI3D.sh
echo "OK"

cd madagascar
scons
cd -

# Run the code
echo "Run code"
./startSOFI3D.sh 8
echo "OK"

# Merge snapshots into visualizable ones
echo "Prepare snapshots"
../bin/snapmerge in_and_out/sofi3D.json
echo "OK"

cd ../..

echo "Done."
echo "Run MATLAB script ./ASOFI3D/mfiles/snap3D_allplanes.m, to see the wavefield"
