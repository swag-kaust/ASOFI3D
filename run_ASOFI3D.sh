#!/usr/bin/env bash
# Compile and run ASOFI3D solver and then prepare snapshots for visualization.

# Correctly determine the exit status of a pipe.
set -o pipefail

# Separator line to separate output of different stages.
sep="************************************************************************"

# Number of the MPI processes. Default value is adapted for Macbook Pro laptop.
nmpiprocs=4
if [[ $# -gt 0 ]]; then
    nmpiprocs="$1"
fi

# Go to directory with stock-given scripts
cd ASOFI3D/par

# Compile the whole code
#echo "Compilation"
#./compileSOFI3D.sh
#echo "OK"

#cd madagascar
#scons
#cd -

# Run the code
printf "Run code\n"
./startSOFI3D.sh "$nmpiprocs"
if [ $? -eq 0 ]; then
    printf "OK\n"
else
    printf "ERROR: running ASOFI3D solver has failed. Check the output above\n"
    exit 1
fi

# Merge snapshots made by individual MPI processes for visualization.
printf "%s\n" "$sep"
printf "Prepare snapshots\n"
../bin/snapmerge in_and_out/sofi3D.json
if [ $? -eq 0 ]; then
    printf "OK\n"
else
    printf "ERROR: merging snapshots has failed\n"
    exit 1
fi

cd ../..

printf "%s\n" "$sep"
printf "Done.\n"
script="./ASOFI3D/mfiles/snap3D_allplanes.m"
printf "Run MATLAB script %s to see the wavefield.\n" "$script"
