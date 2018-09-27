#!/bin/bash
#
#PBS -l nodes=100:ppn=1,walltime=48:00:00
#PBS -q long
#
pbsdsh "rm /tmp/mpd2.console_koed" 
mpdboot -n `sort -u $PBS_NODEFILE | wc -l` -f $PBS_NODEFILE
cd $WORK/sofi3D/par
mpdrun -np `cat $PBS_NODEFILE | wc -l` ../bin/sofi3D_acoustic ./in_and_out/sofi3D.inp > ./in_and_out/sofi3D.out
mpdallexit
