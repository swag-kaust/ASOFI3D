#!/bin/bash
#
#PBS -l nodes=100:ppn=1,walltime=48:00:00
#PBS -q long
#
pbsdsh rm /tmp/mpd2.console_koed 
#mpdboot -n `sort -u $PBS_NODEFILE | wc -l` -f $PBS_NODEFILE
#mpdallexit

