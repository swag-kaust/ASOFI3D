#!/bin/bash

# Just compile + run + prepare snapshots for visualization

# Load MPI modules in case if you are using KAUST workstation
source mpiInit.sh

# Go to directory with stock-given scripts
cd SOFI3D_org/par

# Compile the whole code
echo "Compilation"
./compileSOFI3D.sh
echo "OK"

# Run the code
echo "Run code"
./startSOFI3D.sh
echo "OK"

# Merge snapshots into visualizable ones
echo "Prepare snapshots"
../bin/snapmerge in_and_out/sofi3D.json
echo "OK"

cd ../..

echo "Done."
echo "Run MATLAB script ./SOFI3D_org/mfiles/snap3D_allplanes.m, to see the wavefield"