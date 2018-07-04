#shell script for starting a benchmark run with $1 processors
#$1 is a parameter that has to be passed through by e.g. the command line "./startBENCH.sh 8"

#----excecute with OPENMPI
mpirun -np $1 nice -19 ../bin/sofi3D ./in_and_out/benchmark/speedup_np$1.json | tee ./in_and_out/benchmark/speedup_np$1.jout
#mpirun -np $1 nice -19 ../bin/sofi3D ./in_and_out/benchmark/scaleup_np$1.json | tee ./in_and_out/benchmark/scaleup_np$1.jout
#mpirun --hostfile mpihosts -np $1 nice -19 ../bin/sofi3D ./in_and_out/benchmark/speedup_np$1.json | tee ./in_and_out/benchmark/speedup_np$1.jout
#mpirun --hostfile mpihosts -np $1 nice -19 ../bin/sofi3D ./in_and_out/benchmark/scaleup_np$1.json | tee ./in_and_out/benchmark/scaleup_np$1.jout


#----excecute with LAMMPI
#lamboot -v mpihosts
#lamboot 
#mpirun -np $1 nice -19 ../bin/sofi3D ./in_and_out/benchmark/speedup_np$1.json | tee ./in_and_out/benchmark/speedup_np$1.jout
#mpirun -np $1 nice -19 ../bin/sofi3D ./in_and_out/benchmark/scaleup_np$1.json | tee ./in_and_out/benchmark/scaleup_np$1.jout

