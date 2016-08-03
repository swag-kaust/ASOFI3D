/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI3D.
 * 
 * SOFI3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI3D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   globvar.h - global variables used in SOFI3D
 *
 *  ----------------------------------------------------------------------*/

/* definition of global variables used in the finite difference programs*/
/* For the names of the global variables
   uppercase letters are used */

float DX=0.0, DY=0.0, DZ=0.0, TIME=0.0, DT=0.0, TS=0.0, PLANE_WAVE_DEPTH=0.0, PLANE_WAVE_ANGLE=0.0;
float TSNAP1=0.0, TSNAP2=0.0, TSNAPINC=0.0, *FL=NULL, TAU=0.0, FREF=0.0, REC_ARRAY_DEPTH=0.0, REC_ARRAY_DIST=0.0;
float XREC1=0.0, XREC2=0.0, YREC1=0.0, YREC2=0.0, ZREC1=0.0, ZREC2=0.0;
float REFREC[4]={0.0, 0.0, 0.0, 0.0}, DAMPING=8.0, VPPML=0.0, FPML=0.0, NPOWER=2.0, K_MAX_CPML=10.0;
float SOURCE_ALPHA=0.0, SOURCE_BETA=0.0;
float AMON=0.0, STR=0.0, DIP=0.0, RAKE=0.0;
int   SEISMO=0, NDT=1, NDTSHIFT=0, NGEOPH=0, SEIS_FORMAT[6]={0, 0, 0, 0, 0, 0}, FREE_SURF=0, READMOD=0;
int   MOD_FORMAT[6]={0, 0, 0, 0, 0, 0}, READREC=0, REC_ARRAY=0, LOG=0, FDORDER=2, FDORDER_TIME=2, FW=0, ABS_TYPE=0, BLOCK=0;
int   NX=1, NY=1, NZ=1, NT=0, SOURCE_SHAPE=0, SOURCE_TYPE=0, SNAP=0, SNAP_FORMAT=0, BOUNDARY=0, SRCREC=0;
int   CHECKPTREAD=0, CHECKPTWRITE=0, SNAP_PLANE=0;
int   NXG=1, NYG=1, NZG=1, IDX=1, IDY=1, IDZ=1, L=1, NX1=1, NX2=1, NY1=1, NY2=1, NZ1=1, NZ2=1, DRX=0, DRZ=0;
int   RUN_MULTIPLE_SHOTS=0, FDCOEFF=0, WRITE_MODELFILES=2;
int   OUTNTIMESTEPINFO=1; /*every OUTNTIMESTEPINFO th timestep, information on the time step will be given to screen/file */
int   OUTSOURCEWAVELET=0;
char  SNAP_FILE[STRING_SIZE]="", SOURCE_FILE[STRING_SIZE]="", SIGNAL_FILE[STRING_SIZE]="";
char  MFILE[STRING_SIZE]="", REC_FILE[STRING_SIZE]="", LOG_FILE[STRING_SIZE]="", CHECKPTFILE[STRING_SIZE]="";
char  SEIS_FILE[STRING_SIZE]="";
char  FILEINP[STRING_SIZE]; /* input file name (appears in SEG-Y header) */
FILE  *FP=NULL;
int   LITTLEBIG=0, ASCIIEBCDIC=0, IEEEIBM=0; /* computer charcteristics */
int   SOFI3DVERS; /* version of SOFI3D 33: current 3D isotropic elastic (SSG), version of SOFI3D 32: current 3D isotropic acoustic (SSG) */

/* Mpi-variables */
int   NP, NPSP, NPROC, NPROCX, NPROCY, NPROCZ, MYID, IENDX, IENDY, IENDZ;
int   POS[4], INDEX[7];     
const int TAG1=1,TAG2=2, TAG3=3, TAG4=4, TAG5=5,TAG6=6; 


float FC=0.0,AMP=1.0, REFSRC[3]={0.0, 0.0, 0.0}, SRC_DT, SRCTSHIFT=0.0;
int SRC_MF=0, SIGNAL_FORMAT[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
int SRCOUT_PAR[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, FSRC=1, JSRC=2147483647, LSRC=0;
char SRCOUT_FILE[STRING_SIZE]="";
