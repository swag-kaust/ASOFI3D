# Below you can set C compiler and necessary compiler flags by hand.
# Modify only if you know what you're doing.

# On Linux cluster running LAM
# CC = hcc
# LDLIBS = -lm -lmpi
# CFLAGS = -g -Wall -O4

# On CRAY T3E
# CC=cc

# On Linux Cluster using OpenMPI
# CC=mpicc
# -ta=tesla:managed -acc -Minfo=accel
# CC = mpicc -ta=multicore -acc -Minfo=accel
# LDLIBS = -lm
# CFLAGS = -Wall -O3

# On InstitutsCluster2 using intelMPI AND ITAC (no optimization)
# CC = mpicc
# LDLIBS = -tcollect
# CFLAGS = -tcollect -trace

# On InstitutsCluster2 using Open MPI and scalasca (no optimization)
# CC = skin mpicc
# LDLIBS =
# CFLAGS =-O3

# On InstitutsCluster2 using Open MPI and DDT (no optimization)
# CC = mpicc
# LDLIBS = -g
# CFLAGS = -g -O0 -i_dynamic

# On JUROPA cluster
# CC=mpicc
# LDLIBS=-lm
# CFLAGS=-O3 -ipo

# On HLRN system
# CC=mpcc
# LDLIBS=-lm

# On Linux cluster ALTIX.RZ.UNI-KIEL:DE
# CC=icc
# CFLAGS=-mp -O3 -ip0q
# LDLIBS=-lmpi -lm

# On Linux cluster ALTIX.HRZ.TU-FREIBERG.DE
# CC=icc
# CFLAGS=-mp -O3 -ip0
# LDLIBS=-lmpi -lm -i-static
