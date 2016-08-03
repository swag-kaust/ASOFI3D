#----excecute with LAMMPI
#lamboot -v mpihosts
#lamboot

#mpirun -np 8 nice -19 ../bin/sofi3D_acoustic ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout
#mpirun -np 8 nice -19 ../bin/sofi3D ./in_and_out/sofi3D.inp | tee ./in_and_out/sofi3D.out
mpirun -np 8 nice -19 ../bin/sofi3D ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout

#----execute with OPENMPI2
#mpirun -np 8 nice -19 ../bin/sofi3D ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout
#mpirun --hostfile mpihosts -np 8 nice -19 ../bin/sofi3D ./in_and_out/sofi3D.json | tee ./in_and_out/sofi3D.jout

#merge snapshots
#../bin/snapmerge ./in_and_out/sofi3D.json
