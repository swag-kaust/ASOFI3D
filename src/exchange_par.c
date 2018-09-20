/*------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs                         
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){
	/* declaration of extern variables */

	extern float DX, DY, DZ, TIME, DT, *FL, TS, TAU, FREF, PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE;
	extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2;
	extern float SOURCE_ALPHA, SOURCE_BETA, VPPML;
	extern float AMON, STR, DIP, RAKE;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern float TSNAP1, TSNAP2, TSNAPINC, FW, REFREC[4], DAMPING, FPML, NPOWER, K_MAX_CPML;
	extern int SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT[6], FREE_SURF, READMOD, READREC;
	extern int BOUNDARY, REC_ARRAY, LOG, IDX, IDY, IDZ, ABS_TYPE, WRITE_MODELFILES;
	extern int   NX, NY, NZ, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, SNAP_PLANE;
	extern int DRX, DRZ, L, SRCREC, FDORDER,FDORDER_TIME;
	extern int NPROC,NPROCX,NPROCY,NPROCZ, MYID, CHECKPTREAD, CHECKPTWRITE, RUN_MULTIPLE_SHOTS, FDCOEFF;
	extern int   LITTLEBIG, ASCIIEBCDIC, IEEEIBM;
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern char  RSFDEN[STRING_SIZE]; // RSF
	extern int RSF; // RSF
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE], SEIS_FILE[STRING_SIZE];
	extern char  FILEINP[STRING_SIZE];

    // Model parameters for model generation.
    extern float VPV1, VSV1, EPSX1, EPSY1, DELX1, DELY1, DELXY1;
    extern float GAMX1, GAMY1, RHO1, DH1;
    extern float VPV2, VSV2, EPSX2, EPSY2, DELX2, DELY2, DELXY2;
    extern float GAMX2, GAMY2, RHO2, DH2;


	int idum[NPAR];
	float fdum[NPAR];

	if (MYID == 0){ 
		fdum[1]  = DX;
		fdum[2]  = DY;
		fdum[3]  = DZ;
		fdum[4]  = TIME;
		fdum[5]  = DT;
		fdum[6]  = TS;
		fdum[7]  = PLANE_WAVE_ANGLE;
		fdum[8]  = 0.0;
		fdum[9]  = 0.0;
		fdum[10]  = TAU;
		fdum[11]  = FW;
		fdum[12]  = TSNAP1;
		fdum[13]  = TSNAP2;
		fdum[14]  = TSNAPINC;
		fdum[15]  = REFREC[1];
		fdum[16]  = REFREC[2];
		fdum[17]  = REFREC[3];
		fdum[18]  = XREC1;
		fdum[19]  = YREC1;
		fdum[20]  = ZREC1;
		fdum[21]  = XREC2;
		fdum[22]  = YREC2;
		fdum[23]  = ZREC2;
		fdum[24]  = DAMPING;
		fdum[25]  = FPML;
		fdum[26]  = REC_ARRAY_DEPTH;
		fdum[27]  = REC_ARRAY_DIST;
		fdum[28]  = PLANE_WAVE_DEPTH;
		fdum[29]  = SOURCE_ALPHA;
		fdum[30]  = SOURCE_BETA;
		fdum[31]  = VPPML;
		fdum[32]  = FREF;
		fdum[33]  = NPOWER;
		fdum[34]  = K_MAX_CPML;
		
		fdum[35]  = AMON;
		fdum[36]  = STR;
		fdum[37]  = DIP;
		fdum[38]  = RAKE;

        // Model parameters, top layer.
        fdum[39]  = VPV1;
        fdum[40]  = VSV1;
        fdum[41]  = EPSX1;
        fdum[42]  = EPSY1;
        fdum[43]  = DELX1;
        fdum[44]  = DELY1;
        fdum[45]  = DELXY1;
        fdum[46]  = GAMX1;
        fdum[47]  = GAMY1;
        fdum[48]  = RHO1;
        fdum[49]  = DH1;

        // Model parameters, bottom layer.
        fdum[50]  = VPV2;
        fdum[51]  = VSV2;
        fdum[52]  = EPSX2;
        fdum[53]  = EPSY2;
        fdum[54]  = DELX2;
        fdum[55]  = DELY2;
        fdum[56]  = DELXY2;
        fdum[57]  = GAMX2;
        fdum[58]  = GAMY2;
        fdum[59]  = RHO2;
        fdum[60]  = DH2;



		idum[0]  = FDORDER;
		idum[1]  = NPROCX;
		idum[2]  = NPROCY;
		idum[3]  = NPROCZ;
		idum[4]  = NPROCX*NPROCY*NPROCZ;
		idum[5]  = NX;
		idum[6]  = NY;
		idum[7]  = NZ;
		idum[8]  = SOURCE_SHAPE;
		idum[9]  = SOURCE_TYPE;
		idum[10]  = READMOD;
		idum[11]  = L;
		idum[12]  = FREE_SURF;
		idum[13]  = SNAP;
		idum[14]  = DRX;
		idum[15]  = DRZ;
		idum[16]  = BOUNDARY;
		idum[17]  = REC_ARRAY;
		idum[18]  = SRCREC;
		idum[19]  = LOG;
		idum[20]  = IDX;
		idum[21]  = IDY;
		idum[22]  = IDZ;
		idum[23]  = SNAP_FORMAT;
		idum[24]  = SEISMO;
		idum[25]  = READREC;
		idum[26]  = NGEOPH;
		idum[27]  = NDT;
		idum[28]  = NDTSHIFT;

		idum[29]  = SEIS_FORMAT[0];
		idum[30]  = SEIS_FORMAT[1];
		idum[31]  = SEIS_FORMAT[2];
		idum[32]  = SEIS_FORMAT[3];
		idum[33]  = SEIS_FORMAT[4];
		idum[34]  = SEIS_FORMAT[5];

		idum[35]  = CHECKPTREAD;
		idum[36]  = CHECKPTWRITE;
		idum[37]  = RUN_MULTIPLE_SHOTS;
		idum[38]  = SNAP_PLANE;
		idum[39]  = ABS_TYPE;
		idum[40]  = FDCOEFF;

		idum[41]  = ASCIIEBCDIC;
		idum[42]  = LITTLEBIG;
		idum[43]  = IEEEIBM;

		idum[44]  = WRITE_MODELFILES;
        
		idum[45]  = FDORDER_TIME;
		
		// RSF
		idum[46]  = RSF;


	}

	if (MYID != 0) FL=vector(1,L);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);

	MPI_Bcast(&SOURCE_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

	// RSF
	MPI_Bcast(&RSFDEN,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	//

	MPI_Bcast(&SNAP_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&REC_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&CHECKPTFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

	MPI_Bcast(&FILEINP,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	DX=fdum[1];
	DY=fdum[2];
	DZ=fdum[3];
	TIME=fdum[4];
	DT=fdum[5];
	TS=fdum[6];
	PLANE_WAVE_ANGLE=fdum[7];

	TAU=fdum[10];
	FW=fdum[11];
	TSNAP1=fdum[12];
	TSNAP2=fdum[13];
	TSNAPINC=fdum[14];
	REFREC[1]=fdum[15];
	REFREC[2]=fdum[16];
	REFREC[3]=fdum[17];
	XREC1=fdum[18];
	YREC1=fdum[19];
	ZREC1=fdum[20];
	XREC2=fdum[21];
	YREC2=fdum[22];
	ZREC2=fdum[23];
	DAMPING=fdum[24];
	FPML=fdum[25];
	REC_ARRAY_DEPTH=fdum[26];
	REC_ARRAY_DIST=fdum[27];
	PLANE_WAVE_DEPTH=fdum[28];
	SOURCE_ALPHA = fdum[29];
	SOURCE_BETA = fdum[30];
	VPPML = fdum[31];
	FREF = fdum[32];
	NPOWER = fdum[33];
	K_MAX_CPML = fdum[34];
	AMON = fdum[35];
	STR = fdum[36];
	DIP = fdum[37];
	RAKE = fdum[38];

    // Model parameters, top layer.
    VPV1   = fdum[39];
    VSV1   = fdum[40];
    EPSX1  = fdum[41];
    EPSY1  = fdum[42];
    DELX1  = fdum[43];
    DELY1  = fdum[44];
    DELXY1 = fdum[45];
    GAMX1  = fdum[46];
    GAMY1  = fdum[47];
    RHO1   = fdum[48];
    DH1    = fdum[49];

    // Model parameters, bottom layer.
    VPV2   = fdum[50];
    VSV2   = fdum[51];
    EPSX2  = fdum[52];
    EPSY2  = fdum[53];
    DELX2  = fdum[54];
    DELY2  = fdum[55];
    DELXY2 = fdum[56];
    GAMX2  = fdum[57];
    GAMY2  = fdum[58];
    RHO2   = fdum[59];
    DH2    = fdum[60];


    // -------------------
    // Integer parameters.
	FDORDER = idum[0];
	NPROCX = idum[1];
	NPROCY = idum[2];
	NPROCZ = idum[3];
	NPROC  = idum[4];
	NX = idum[5];
	NY = idum[6];
	NZ = idum[7];
	SOURCE_SHAPE = idum[8];
	SOURCE_TYPE = idum[9];
	READMOD = idum[10];
	L = idum[11];
	FREE_SURF = idum[12];
	SNAP = idum[13];
	DRX = idum[14];
	DRZ = idum[15];
	BOUNDARY = idum[16];
	REC_ARRAY = idum[17];
	SRCREC = idum[18];
	LOG = idum[19];
	IDX = idum[20];
	IDY = idum[21];
	IDZ = idum[22];



	SNAP_FORMAT = idum[23];
	SEISMO = idum[24];
	READREC = idum[25];
	NGEOPH = idum[26];
	NDT = idum[27];
	NDTSHIFT = idum[28];

	SEIS_FORMAT[0] = idum[29];
	SEIS_FORMAT[1] = idum[30];
	SEIS_FORMAT[2] = idum[31];
	SEIS_FORMAT[3] = idum[32];
	SEIS_FORMAT[4] = idum[33];
	SEIS_FORMAT[5] = idum[34];

	CHECKPTREAD = idum[35];
	CHECKPTWRITE = idum[36];
	RUN_MULTIPLE_SHOTS= idum[37];
	SNAP_PLANE= idum[38];
	ABS_TYPE= idum[39];
	FDCOEFF= idum[40];

	ASCIIEBCDIC=idum[41];
	LITTLEBIG=idum[42];
	IEEEIBM=idum[43];

	WRITE_MODELFILES=idum[44];
    
	FDORDER_TIME=idum[45];
	
	
	RSF = idum[46];



	MPI_Bcast(&FL[1],L,MPI_FLOAT,0,MPI_COMM_WORLD);

}

