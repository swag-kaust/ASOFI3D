#!/bin/bash
# Run solver with MPI writing the output to the screen and a file.
# Note that the following directives (starting with #SBATCH)
# are processed by the Slurm job-queueing system.

#SBATCH --account=k1056
#SBATCH --job-name=ASOFI
#SBATCH --time=00-00:01:00      # Time format is DD-HH:MM:SS
#SBATCH --hint=nomultithread
#SBATCH --ntasks=32             # Total number of MPI processes
#SBATCH --ntasks-per-node=32    # Number of MPI processes per node (default 1)
#SBATCH --nodes=1               # Number of nodes


# Note that the `srun` command is recommended instead of the `mpirun`
# when running under Slurm. Number of MPI processes is determined by the above
# directive `--ntasks`.
echo "SLURM_NTASKS = $SLURM_NTASKS" 
srun ../bin/sofi3D ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout

#----excecute with LAMMPI
#lamboot -v mpihosts
#lamboot

#mpirun -np 8 nice -19 ../bin/sofi3D_acoustic ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout

#mpirun -np 8 nice -19 ../bin/sofi3D ./in_and_out/sofi3D.inp | tee ./in_and_out/sofi3D.out

# profiling with nvprof
#mpirun -np 8 nice -19 nvprof --cpu-profiling on ../bin/sofi3D ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout

#merge snapshots
#../bin/snapmerge ./in_and_out/sofi3D.json


