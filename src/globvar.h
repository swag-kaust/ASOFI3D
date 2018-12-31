#ifndef GLOBVAR_H
#define GLOBVAR_H
/*
 * Provide declaration of the global variables.
 * For the names of the global variables uppercase letters are used.
 * 
 * See read_par_json.c for the definition of these variables.
 */
#include <stdio.h>

#include "constants.h"

//Imaging
extern int   RTM_FLAG;


//MADAGASCAR VAR start
extern int  RSF;
extern char RSFDEN[STRING_SIZE];
//MADAGASCAR VAR end


extern float DX, DY, DZ, TIME, DT, TS, PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE;
extern float TSNAP1, TSNAP2, TSNAPINC, *FL, TAU, FREF, REC_ARRAY_DEPTH, REC_ARRAY_DIST;
extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2;
extern float REFREC[4], DAMPING, VPPML, FPML, NPOWER, K_MAX_CPML;
extern float SOURCE_ALPHA, SOURCE_BETA;
extern float AMON, STR, DIP, RAKE;

// Moment tensor components.
extern float M11, M12, M13, M22, M23, M33;

extern int SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT[6], FREE_SURF, READMOD;
extern int READREC, REC_ARRAY, LOG, FDORDER, FDORDER_TIME, FW, ABS_TYPE, BLOCK;
extern int NX, NY, NZ, NT, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, BOUNDARY, SRCREC;
extern int CHECKPTREAD, CHECKPTWRITE, SNAP_PLANE;
extern int NXG, NYG, NZG, IDX, IDY, IDZ, L, NX1, NX2, NY1, NY2, NZ1, NZ2, DRX, DRZ;
extern int RUN_MULTIPLE_SHOTS, FDCOEFF, WRITE_MODELFILES;
extern int OUTNTIMESTEPINFO; /*every OUTNTIMESTEPINFO th timestep, information on the time step will be given to screen/file */
extern int OUTSOURCEWAVELET;

extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
extern char MFILE[STRING_SIZE], REC_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
extern char SEIS_FILE[STRING_SIZE];
extern char FILEINP[STRING_SIZE]; /* input file name (appears in SEG-Y header) */
extern FILE *FP;

// Default for personal computers: Little Endian, ASCII, IEEE, which are
// encoded with "zero" values in the next line.
extern int LITTLEBIG, ASCIIEBCDIC, IEEEIBM;
extern int SOFI3DVERS; /* version of SOFI3D 33: current 3D isotropic elastic (SSG), version of SOFI3D 32: current 3D isotropic acoustic (SSG) */

// MPI variables.
extern int NP, NPSP, NPROC, NPROCX, NPROCY, NPROCZ, MYID, IENDX, IENDY, IENDZ;
extern int POS[4], INDEX[7];
extern const int TAG1, TAG2, TAG3, TAG4, TAG5, TAG6;

extern float FC, AMP, REFSRC[3], SRC_DT, SRCTSHIFT;
extern int SRC_MF, SIGNAL_FORMAT[6];
extern int SRCOUT_PAR[6], FSRC, JSRC, LSRC;
extern char SRCOUT_FILE[STRING_SIZE];

// Model parameters for model generation.
extern float VPV1, VSV1, EPSX1, EPSY1, DELX1, DELY1, DELXY1;
extern float GAMX1, GAMY1, RHO1, DH1;
extern float VPV2, VSV2, EPSX2, EPSY2, DELX2, DELY2, DELXY2;
extern float GAMX2, GAMY2, RHO2, DH2;
#endif /* ifndef GLOBVAR_H */
