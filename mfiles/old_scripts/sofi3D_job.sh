#!/bin/bash
#
#PBS -l nodes=27:ppn=1,walltime=48:00:00
#PBS -q long
#

cd $WORK/sofi3D/par
mpirun -np `cat $PBS_NODEFILE | wc -l` -hostfile $PBS_NODEFILE ../bin/sofi3D_acoustic ./par/in_and_out/sofi3D.inp > ./in_and_out/sofi3D.out
