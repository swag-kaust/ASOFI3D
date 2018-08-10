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

# Working directory for simulation input and output.
sim_dir="par/"
if [[ $# -gt 1 ]]; then
    sim_dir="$2"
fi

# Add `bin` directory to $PATH to be able to execute programs
# (sofi3d, snapmerge, etc.) without specifying the file path.
script_dir="$(cd "$(dirname $0)" && pwd)"
export PATH=${script_dir}/bin:$PATH

# Go to the working directory with simulation input files.
pushd "$sim_dir" > /dev/null || exit 1

# Compile the whole code.
echo "$sep"
echo "Compilation"
make asofi3D
if [ "$?" -eq 0 ]; then
    echo "Compilation successful"
else
    echo "ERROR: Compilation failed" > /dev/stderr
    exit 1
fi

# Run the code
config_file="in_and_out/sofi3D.json"
printf "Run code\n"
mpirun -n $nmpiprocs nice -19 sofi3D $config_file | tee in_and_out/sofi3D.jout
if [ $? -eq 0 ]; then
    echo "OK"
else
    echo "ERROR: running ASOFI3D solver has failed. Check the output above" > /dev/stderr
    exit 1
fi

# Merge snapshots made by individual MPI processes for visualization.
echo "$sep"
echo "Prepare snapshots"
snapmerge $config_file
if [ $? -eq 0 ]; then
    echo "OK"
else
    echo "ERROR: merging snapshots has failed" > /dev/stderr
    exit 1
fi

popd > /dev/null || exit 1

echo "$sep"
echo "Done."
script="./ASOFI3D/mfiles/snap3D_allplanes.m"
printf "Run MATLAB script %s to see the wavefield.\n" "$script"
