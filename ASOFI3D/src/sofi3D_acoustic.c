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
/* ----------------------------------------------------------------------
 * This is program SOFI3D. Acoustic Version
 * Parallel 3-D Viscoelastic Finite Difference Seismic Modelling
 * using the Standard Staggered Grid (SSG)
 *
 * PLEASE DO NOT DISTRIBUTE. PLEASE REFER OTHER PEOPLE TO :
 *
 * Prof. Dr. Thomas Bohlen, Karlsruhe Institute of Technology,
 * Geophysical Institute,
 * Hertzstr. 16, 76187 Karlsruhe, Germany
 * Phone/Fax: +49 (0)721 608 44416
 * mailto:thomas.bohlen@kit.edu,
 * http://www.gpi.kit.edu/
 * http://www.gpi.kit.edu/SOFI3D.php
 *
 *
 * If you want to publish synthetic data calculated with this program please
 * give a reference to the following paper:
 * Bohlen, T., 2002, Parallel 3-D viscoelastic finite-difference seismic modelling,
 * Computers @ Geopsciences, Vol. 28, No. 8, 887-889.
 *
 *  ----------------------------------------------------------------------*/



#include "fd.h"
#include "globvar.h"

int main(int argc, char **argv){
	int ns, nt, nseismograms=0, nf1, nf2, i;
	int lsnap, nsnap=0, lsamp=0, nlsamp=0, buffsize;
	int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, ishot, nshots, h;

	double 	time0=0.0, time1=0.0, time2=0.0, time3=0.03, time4=0.0;
	double * time_v_update, * time_s_update, * time_s_exchange,* time_v_exchange, * time_timestep;
	int * xa, * xb, * ya, * yb, * za, * zb;
	float  ***  sxx=NULL;
	float  ***  vx=NULL, ***  vy=NULL, ***  vz=NULL;
	float  ***  sxx1=NULL, ***  syy1=NULL, ***  szz1=NULL;
	float  ***  vx1=NULL, ***  vy1=NULL, ***  vz1=NULL;

	float  ***  sxx2=NULL, ***  syy2=NULL, ***  szz2=NULL;
	float  ***  vx2=NULL, ***  vy2=NULL, ***  vz2=NULL;

	float  ***  sxx3=NULL, ***  syy3=NULL, ***  szz3=NULL;
	float  ***  vx3=NULL, ***  vy3=NULL, ***  vz3=NULL;

	float  ***  rho, ***  pi;
	float  *** absorb_coeff=NULL, *** absorb_coeffx=NULL, *** absorb_coeffy=NULL, *** absorb_coeffz=NULL;
	float  ** srcpos=NULL, **srcpos_loc=NULL, ** srcpos1=NULL, ** signals=NULL;
	int    ** recpos=NULL, ** recpos_loc=NULL;
	float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionvz=NULL, ** sectionp=NULL,
			** sectioncurl=NULL, ** sectiondiv=NULL;
	float *** bufferlef_to_rig, *** bufferrig_to_lef;
	float *** buffertop_to_bot, *** bufferbot_to_top;
	float *** bufferfro_to_bac, *** bufferbac_to_fro;

	float *** sbufferlef_to_rig, *** sbufferrig_to_lef;
	float *** sbuffertop_to_bot, *** sbufferbot_to_top;
	float *** sbufferfro_to_bac, *** sbufferbac_to_fro;


	float ** seismo_fulldata=NULL;
	int * recswitch=NULL;

	MPI_Request *req_send, *req_rec, *sreq_send, *sreq_rec;
	MPI_Status  *send_statuses, *rec_statuses;

	float memdyn, memmodel, memseismograms, membuffer, memtotal;
	float fac1, fac2;
	char *buff_addr, ext[10];
	int * stype=NULL, * stype_loc=NULL;
	char buffer[STRING_SIZE], bufferstring[10];
	/*char cline[256];*/
	FILE * fpsrc=NULL;

	/* Initialize MPI environment */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&NP);
	MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

	setvbuf(stdout, NULL, _IONBF, 0);


	/* initialize clock for estimating runtime of program */	
	if (MYID == 0){
		time1=MPI_Wtime();
		clock();
	}


	/* print program name, version etc to stdout*/
	if (MYID == 0) info(stdout);
	SOFI3DVERS=32; /* 3D isotropic acoustic */

	/* PE 0 is reading the parameters from the input file sofi3D.inp */
	if (MYID == 0){
		if (argc>1) {
			strncpy(FILEINP,argv[1],STRING_SIZE);
			fprintf(stderr," Input parameter filename read from command line : %s. \n\n",FILEINP);
			if (strlen(FILEINP)>STRING_SIZE-1) {
				fprintf(stderr,"\n SOFI3D cannot handel pathes with more than %d characters.\n",STRING_SIZE-1);
				fprintf(stderr," Error: SOFI3D could not read input parameter file name. -> Exit. \n\n");
				return -1;
			}
		}
		else {
			strcpy(FILEINP,"sofi3D.inp");
			fprintf(stderr," Caution: input parameter filename set to default 'sofi3D.inp'. \n\n");
		}
		FP=fopen(FILEINP,"r");

		//read json formated input file
		read_par_json(stdout, FILEINP);
		fclose(FP);
	}
	/* PE 0 will broadcast the parameters to all others PEs */
	exchange_par();

	if (MYID == 0) note(stdout);

	sprintf(ext,".%i",MYID);  
	strcat(LOG_FILE,ext);

	if (MYID==0) /*  define LOG_FILE for MYID=0 */
		switch (LOG){
		case 0 : FP=fopen("/dev/null","w"); /* no logging information will be output */
		break;	
		case 1 : FP=stdout; /* logging information will be written to standard output */
		break;	
		case 2 : FP=fopen(LOG_FILE,"w"); /* logging information will be written to LOG_FILE */
		break;	
		}

	if (MYID>0) FP=fopen(LOG_FILE,"w"); /* all other nodes write logging info to LOG_FILE */

	fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);

	/* domain decomposition */
	initproc();

	/* set some time counters */
	NT=(int)ceil(TIME/DT); /* number of timesteps - replaces: NT=iround(TIME/DT); */
	TIME=(NT-1)*DT; /* TIME set to true time of the last time step */
	if (NDTSHIFT>NT) ns=0;
	else ns=(int)ceil((float)(NT-NDTSHIFT)/(float)NDT); /* number of samples per trace - replaces buggy formula: ns=iround(NT-NDTSHIFT/NDT); */
	lsnap=iround(TSNAP1/DT);  /* first snapshot at this timestep */

	/* output of parameters to stdout: */
	/*if (MYID==0) writepar(FP,recpos, ns);*/
	if (MYID==0) writepar(FP, ns);
	/* NXG, NYG NZG denote size of the entire (global) grid */
	NXG=NX;
	NYG=NY;
	NZG=NZ;

	/* In the following, NX, MY, NZ denote size of the local grid ! */
	NX = IENDX;
	NY = IENDY;
	NZ = IENDZ;

	/* compute receiver locations within each subgrid and
	   store local receiver coordinates in recpos_loc */	
	if (SEISMO){
		fprintf(FP,"\n ------------------ READING RECEIVER PARAMETERS ----------------- \n");
		recpos=receiver(FP,&ntr);
		recswitch = ivector(1,ntr);
		recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
		ntr_glob=ntr;
		ntr=ntr_loc;
	}


	/* allocate buffer for seismogram output, merged seismogram section of all PEs */
	if (SEISMO) seismo_fulldata=fmatrix(1,ntr_glob,1,ns);

	/* allocate buffer for seismogram output, seismogram section of each PE */
	/* allocation of memory for seismogramm merge */
	switch (SEISMO){
	case 1 : /* particle velocities only */
		nseismograms=3;
		break;
	case 2 : /* pressure only */
		nseismograms=1;
		break;
	case 3 : /* curl and div only */
		nseismograms=2;
		break;
	case 4 : /* everything */
		nseismograms=6;
		break;
	default : nseismograms=1;
	break;
	}

	/*allocate memory for dynamic, static and buffer arrays */
	fac1=(NZ+FDORDER)*(NY+FDORDER)*(NX+FDORDER);
	fac2=sizeof(float)*pow(2.0,-20.0);

	if(ABS_TYPE==1){ 
		memdyn=22.0*fac1*fac2;
		memmodel=5.0*fac1*fac2;
		memseismograms=nseismograms*ntr*ns*fac2;
		membuffer=2.0*5.0*((NX*NZ)+(NY*NZ)+(NX*NY))*fac2;
		buffsize=(FDORDER/2)*4.0*6.0*(max((NX*NZ),max((NY*NZ),(NX*NY))))*sizeof(MPI_FLOAT);
		memtotal=memdyn+memmodel+memseismograms+membuffer+(buffsize*pow(2.0,-20.0));


		if (MYID==0){
			fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
			fprintf(FP," Size of local grids: NX=%d \t NY=%d \t NZ=%d \n",NX,NY,NZ);
			fprintf(FP," Each process is now trying to allocate memory for:\n");
			fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
			fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
			fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
			fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
			fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
			fprintf(FP," ------------------------------------------------ \n");
			fprintf(FP," Total memory required: \t %6.2f MB.\n\n", memtotal);
		}}

	if(ABS_TYPE==2){
		memdyn=4.0*fac1*fac2;
		memmodel=3.0*fac1*fac2;
		memseismograms=nseismograms*ntr*ns*fac2;
		membuffer=2.0*5.0*((NX*NZ)+(NY*NZ)+(NX*NY))*fac2;
		buffsize=(FDORDER/2)*4.0*6.0*(max((NX*NZ),max((NY*NZ),(NX*NY))))*sizeof(MPI_FLOAT);
		memtotal=memdyn+memmodel+memseismograms+membuffer+(buffsize*pow(2.0,-20.0));


		if (MYID==0){
			fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
			fprintf(FP," Size of local grids: NX=%d \t NY=%d \t NZ=%d \n",NX,NY,NZ);
			fprintf(FP," Each process is now trying to allocate memory for:\n");
			fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
			fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
			fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
			fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
			fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
			fprintf(FP," ------------------------------------------------ \n");
			fprintf(FP," Total memory required: \t %6.2f MB.\n\n", memtotal);
		}}


	/* allocate buffer for buffering messages */
	buff_addr=malloc(buffsize);
	if (!buff_addr) err("allocation failure for buffer for MPI_Bsend !");
	MPI_Buffer_attach(buff_addr,buffsize);


	/* allocation for request and status arrays */
	req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	sreq_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	sreq_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
	rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));


	/* allocation for timing arrays used for performance analysis */
	time_v_update=dvector(1,NT);
	time_s_update=dvector(1,NT);
	time_s_exchange=dvector(1,NT);
	time_v_exchange=dvector(1,NT);
	time_timestep=dvector(1,NT);


	/* memory allocation for dynamic (wavefield) arrays */
	vx  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
	vy  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
	vz  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
	sxx =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);


	/* memory allocation for static (model) arrays */
	rho =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	pi  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	absorb_coeff=  f3tensor(1,NY,1,NX,1,NZ);

	/* PML indices*/
	xa=ivector(0,5);
	xb=ivector(0,5);
	ya=ivector(0,5);
	yb=ivector(0,5);
	za=ivector(0,5);
	zb=ivector(0,5);

	/* memory allocation for PML-splitting variables */
	if(ABS_TYPE==1){

		vx1  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		vy1  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		vz1  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		sxx1 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		syy1 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		szz1 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

		vx2  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		vy2  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		vz2  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		sxx2 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		syy2 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		szz2 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

		vx3  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		vy3  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		vz3  =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		sxx3 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		syy3 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		szz3 =  f3tensor(1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

		absorb_coeffx=  f3tensor(1,NY,1,NX,1,NZ);
		absorb_coeffy=  f3tensor(1,NY,1,NX,1,NZ);
		absorb_coeffz=  f3tensor(1,NY,1,NX,1,NZ);

	}


	/* memory allocation for buffer arrays in which the wavefield
	   information which is exchanged between neighbouring PEs is stored */


	/* number of wavefield parameters that need to be exchanged - see exchange_*_acoustic.c */
	nf1=(FDORDER/2)-1;
	nf2=FDORDER/2;

	bufferlef_to_rig = f3tensor(1,NY,1,NZ,1,nf1);
	bufferrig_to_lef = f3tensor(1,NY,1,NZ,1,nf2);
	buffertop_to_bot = f3tensor(1,NX,1,NZ,1,nf1);
	bufferbot_to_top = f3tensor(1,NX,1,NZ,1,nf2);
	bufferfro_to_bac = f3tensor(1,NY,1,NX,1,nf1);
	bufferbac_to_fro = f3tensor(1,NY,1,NX,1,nf2);


	sbufferlef_to_rig = f3tensor(1,NY,1,NZ,1,nf2);
	sbufferrig_to_lef = f3tensor(1,NY,1,NZ,1,nf1);
	sbuffertop_to_bot = f3tensor(1,NX,1,NZ,1,nf2);
	sbufferbot_to_top = f3tensor(1,NX,1,NZ,1,nf1);
	sbufferfro_to_bac = f3tensor(1,NY,1,NX,1,nf2);
	sbufferbac_to_fro = f3tensor(1,NY,1,NX,1,nf1);

	if ((ntr>0)){
		switch (SEISMO){
		case 1 : /* particle velocities only */
			sectionvx=fmatrix(1,ntr,1,ns);
			sectionvy=fmatrix(1,ntr,1,ns);
			sectionvz=fmatrix(1,ntr,1,ns);
			break;
		case 2 : /* pressure only */
			sectionp=fmatrix(1,ntr,1,ns);
			break;
		case 3 : /* curl and div only */
			sectioncurl=fmatrix(1,ntr,1,ns);
			sectiondiv=fmatrix(1,ntr,1,ns);
			break;
		case 4 : /* everything */
			sectionvx=fmatrix(1,ntr,1,ns);
			sectionvy=fmatrix(1,ntr,1,ns);
			sectionvz=fmatrix(1,ntr,1,ns);
			sectioncurl=fmatrix(1,ntr,1,ns);
			sectiondiv=fmatrix(1,ntr,1,ns);
			sectionp=fmatrix(1,ntr,1,ns);
			break;
		}
	}	

	/* memory for source position definition */
	srcpos1=fmatrix(1,6,1,1);

	if (MYID==0) 
		fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);

	/* calculate 3-D array for exponential damping of reflections
	   at the edges of the numerical mesh (PML-boundary)*/
	if(ABS_TYPE==1){
		absorb_PML(absorb_coeffx, absorb_coeffy, absorb_coeffz);
	}

	/* calculate 3-D array for exponential damping of reflections
	   at the edges of the numerical mesh (ABS-boundary) */
	if(ABS_TYPE==2){
		absorb(absorb_coeff);
	}

	/* initialisation of PML and ABS domain */
	if(ABS_TYPE==1){    
		PML_ini_acoustic(xa,xb,ya,yb,za,zb);
		printf("BLOCK = %d \n",BLOCK);
	}

	if(ABS_TYPE==2){    
		xa[0]=1;
		xb[0]=NX;
		ya[0]=1;
		yb[0]=NY;
		za[0]=1;
		zb[0]=NZ;
	}


	/* Reading source positions from SOURCE_FILE */
	fprintf(FP,"\n ------------------ READING SOURCE PARAMETERS ------------------- \n");
	if (PLANE_WAVE_DEPTH>0) {



		/*srcpos=pwsources(&nsrc);*/

		/* for unknown reasons, the pointer does not point to memory that has been allocated by a subroutine this way */
		/*stype=(int *)malloc(nsrc*sizeof(int)); 
		stype_loc=(int *)malloc(nsrc*sizeof(int));*/

		/* I replaced malloc with ivector and started with ishot=1 */
		stype=ivector(1,nsrc);
		stype_loc=ivector(1,nsrc);

		for (ishot=1;ishot<=nsrc;ishot++){
			stype[ishot]=SOURCE_TYPE;
		}	
	}
	else {
		if (MYID==0) switch (SRCREC) {
		case 1:
			fprintf(FP,"\n Reading source parameters from file: %s (SOFI3D source format)\n",SOURCE_FILE);
			if ((fpsrc=fopen(SOURCE_FILE,"r"))==NULL) err(" Source file could not be opened !");

			while(fgets(buffer, STRING_SIZE, fpsrc))
			{
				sscanf(buffer,"%s",bufferstring);
				/* checks if the line contains a '%'character which indicates a comment line,
					and if the reading of a string was successful, which is not the case for an empty line*/
				if ((strchr(buffer,'%')==0) && (sscanf(buffer,"%s",bufferstring)==1)) ++(nsrc);	
			}
			rewind(fpsrc);
			if ((nsrc)==0) fprintf(FP,"\n WARNING: Could not determine number of sources parameter sets in input file. Assuming %d.\n",(nsrc=0));
			else fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,nsrc);

			/*fgets(cline,255,fpsrc);
			if (sscanf(cline,"%d",&nsrc)==0) fprintf(FP,"\n WARNING: Could not determine number of sources parameter sets in input file. Assuming %d.\n",(nsrc=0));
			else fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,nsrc);*/
			break;
		case 2: fprintf(FP,"\n Reading source and receiver parameters from file: %s (UKOOA format)\n",SOURCE_FILE);
		err("\n  UNDER CONSTRUCTION !!! \n");
		break;

		case 3: fprintf(FP,"\n Reading source parameters from file: %s\n",SOURCE_FILE);
		err("\n  UNDER CONSTRUCTION !!! \n");
		break;		

		default: fprintf(FP,"\n WARNING: Format of source file %s unknown! Using default parameters instead of source file. \n",SOURCE_FILE);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);

		/* for unknown reasons, the pointer does not point to memory that has been allocated by a subroutine this way */
		/*stype=(int *)malloc(nsrc*sizeof(int));
	    stype_loc=(int *)malloc(nsrc*sizeof(int));*/

		/* I replaced malloc with ivector and started with ishot=1 */
		stype=ivector(1,nsrc);
		stype_loc=ivector(1,nsrc);

		srcpos=sources(fpsrc,&nsrc,stype);
	}

	if (stype==NULL) fprintf(stderr,"PE%d: Source type(s) undefined?! \n",MYID);

	if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;	

	for (ishot=1;ishot<=nshots;ishot++){

		fprintf(FP,"\n MYID=%d *****  Starting simulation for shot %d of %d  ********** \n",MYID,ishot,nshots);
		for (i=1;i<=6;i++) srcpos1[i][1]=srcpos[i][ishot];
		if (RUN_MULTIPLE_SHOTS){

			/* find this single source positions on subdomains */
			if (nsrc_loc>0) free_matrix(srcpos_loc,1,6,1,1);
			srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1, stype_loc, stype);
			if (stype_loc==NULL) stype_loc = ivector(1,nsrc);;
			srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1, stype_loc, stype);
		}


		else{
			/* Distribute multiple source positions on subdomains */
			if (stype_loc==NULL) stype_loc = ivector(1,nsrc);
			srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc, stype_loc, stype);
		}



		/* create model grids */
		if (READMOD) readmod_acoustic(rho,pi,ishot);
		else model_acoustic(rho,pi);


		/* check if the FD run will be stable and free of numerical dispersion */
		checkfd_acoustic(FP,rho,pi,srcpos,nsrc,recpos,ntr_glob);

		/* For the calculation of the material parameters beteween gridpoints
	   the have to be averaged. For this, values lying at 0 and NX+1,
		for example, are required on the local grid. These are now copied from the
		neighbouring grids */
		matcopy_acoustic(rho,pi);


		/* spatial averaging of material parameters, i.e. density */
		/*av_mat_acoustic(rho,rjp,rkp,rip);*/


		MPI_Barrier(MPI_COMM_WORLD);

		/* use of of checkpoint files is temporarily disabled
		 * there are quite some variables that are not in use in sofi3D_acoustic, e.g. syy or szz!
		if (CHECKPTREAD){
			if (MYID==0){
				time3=MPI_Wtime();
				fprintf(FP," Reading wavefield from check-point file %s \n",CHECKPTFILE);
			}

			read_checkpoint(-1, NX+2, -1, NY+2, -1, NZ+2, vx,vy,vz,sxx,syy,szz,sxx,syy,szz);
		}*/



		/* comunication initialisation for persistent communication */
		/*comm_ini(bufferlef_to_rig, bufferrig_to_lef,
	    buffertop_to_bot, bufferbot_to_top, bufferfro_to_bac,
	    bufferbac_to_fro, req_send, req_rec);

	comm_ini_acoustic(sbufferlef_to_rig, sbufferrig_to_lef,
	    sbuffertop_to_bot, sbufferbot_to_top, sbufferfro_to_bac,
	    sbufferbac_to_fro, sreq_send, sreq_rec);*/



		/* calculate wavelet for each source point */
		signals=wavelet(srcpos_loc,nsrc_loc);


		/* initialize wavefield with zero */
		zero_acoustic(1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2,vx,vy,vz,sxx);

		/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
		/* start of loop over time steps */

		if (MYID==0){
			time2=MPI_Wtime();
			fprintf(FP,"\n\n\n *********** STARTING TIME STEPPING ***************\n");
			fprintf(FP," real time before starting time loop: %4.2f s.\n",time2-time1);
		}

		lsamp=NDTSHIFT+1;
		nlsamp=1;
		for (nt=1;nt<=NT;nt++){

			time_v_update[nt]=0.0;
			time_s_update[nt]=0.0;

			/* Check if simulation is still stable */
			if (isnan(vy[NY/2][NX/2][NZ/2])) err(" Simulation is unstable !");

			if (LOG)
				if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
					fprintf(FP,"\n Computing timestep %d of %d for shot %d \n",nt,NT,ishot);
					time0=MPI_Wtime();
				}

			/* update PML boundaries */
			/* --------------------- */
			/* if BLOCK == 3 == corner, == 1 == face, == 2 == edge */


			/*if((ABS_TYPE == 1)&&((BLOCK == 3)||(BLOCK == 1)||(BLOCK == 2))){*/


			if((ABS_TYPE == 1)&&(BLOCK != 0)){

				/* update of particle velocities */

				/* update NON PML boundaries */
				time_v_update[nt]+=update_v_acoustic(xa[0],xb[0],ya[0],yb[0],za[0],zb[0],nt,vx,vy,vz,sxx,rho,srcpos_loc,signals,nsrc_loc,absorb_coeff,stype_loc);


				/* update PML boundaries */
				for(h=1;h<=BLOCK;h++){

					time_v_update[nt]+=update_v_acoustic_PML(xa[h],xb[h],ya[h],yb[h],za[h],zb[h],nt,vx,vy,vz,sxx,
							vx1,vy2,vz3,rho,srcpos_loc,signals,nsrc_loc,absorb_coeffx,
							absorb_coeffy,absorb_coeffz,stype_loc);
				}

				/* exchange values of particle velocities at grid boundaries between PEs */
				time_v_exchange[nt]=exchange_v(nt, vx, vy, vz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top,
						bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec);


				/* update of components of stress tensor */

				/* update NON PML boundaries */
				time_s_update[nt]+=update_s_acoustic(xa[0],xb[0],ya[0],yb[0],za[0],zb[0],nt,vx,vy,vz,sxx,pi);

				/* update PML boundaries */
				for(h=1;h<=BLOCK;h++){

					time_s_update[nt]+=update_s_acoustic_PML(xa[h],xb[h],ya[h],yb[h],za[h],zb[h],nt,vx,vy,vz,sxx,
							sxx1,sxx2,sxx3,pi,absorb_coeffx,absorb_coeffy,absorb_coeffz);
				}

				/* exchange values of stress at boundaries between PEs */
				time_s_exchange[nt]=exchange_s_acoustic(nt,sxx,sbufferlef_to_rig, sbufferrig_to_lef,
						sbuffertop_to_bot, sbufferbot_to_top, sbufferfro_to_bac, sbufferbac_to_fro, sreq_send, sreq_rec);

			}

			/* if BLOCK == 0 == calculation area or old ABS-Boundary*/
			/* -----------------------------------------------------*/

			if((BLOCK == 0)||(ABS_TYPE==2)){

				/* update of particle velocities */

				/* update NON PML boundaries */
				time_v_update[nt]+=update_v_acoustic(xa[0],xb[0],ya[0],yb[0],za[0],zb[0],nt,vx,vy,vz,sxx,
						rho,srcpos_loc,signals,nsrc_loc,absorb_coeff, stype_loc);


				/* exchange values of particle velocities at grid boundaries between PEs */
				time_v_exchange[nt]=exchange_v(nt, vx, vy, vz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top,
						bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec);

				/* update of components of stress tensor */

				/* update NON PML boundaries */
				time_s_update[nt]+=update_s_acoustic(xa[0],xb[0],ya[0],yb[0],za[0],zb[0],nt,vx,vy,vz,sxx,pi);


				/* exchange values of stress at boundaries between PEs */
				time_s_exchange[nt]=exchange_s_acoustic(nt,sxx,sbufferlef_to_rig, sbufferrig_to_lef,
						sbuffertop_to_bot, sbufferbot_to_top, sbufferfro_to_bac, sbufferbac_to_fro, sreq_send, sreq_rec);

			}

			/* explosive source */
			if (!CHECKPTREAD)
				psource_acoustic(nt,sxx,srcpos_loc,signals,nsrc_loc,stype_loc);


			/* stress free surface ? */
			if ((FREE_SURF) && (POS[2]==0))
				surface_acoustic(1,pi,sxx,vx,vy,vz);


			/* store amplitudes at receivers in sectionvx-sectionvz */
			if ((SEISMO) && (ntr>0) && (nt==lsamp)){
				seismo_acoustic(nlsamp,ntr,recpos_loc,sectionvx,sectionvy,sectionvz,
						sectiondiv,sectioncurl,sectionp,vx,vy,vz,sxx,pi);
				nlsamp++;
				lsamp+=NDT;
			}


			/* save snapshot in file */
			if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
				snap_acoustic(FP,nt,++nsnap,SNAP_FORMAT,SNAP,vx,vy,vz,sxx,pi,IDX,IDY,IDZ,1,1,1,NX,NY,NZ);
				lsnap=lsnap+iround(TSNAPINC/DT);
			}

			if (LOG)
				if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
					time3=MPI_Wtime();
					time_timestep[nt]=(time3-time0);
					fprintf(FP," total real time for timestep %d : \t\t %4.2f s.\n",nt,time3-time0);
				}




		} /* end of loop over timesteps */
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
		fprintf(FP, "\n\n *********** Finish TIME STEPPING ****************\n");
		fprintf(FP, " **************************************************\n\n");


		/* write seismograms to file(s) */
		if (SEISMO){

			/* saves seismograms portion of each PE individually to file */
			//if (ntr>0) saveseis(FP,sectionvx,sectionvy,sectionvz,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns);

			/* merge of seismogram data from all PE and output data collectively */
			switch (SEISMO){
			case 1 : /* particle velocities only */
				catseis(sectionvx, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,1);
				catseis(sectionvy, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,2);
				catseis(sectionvz, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,3);

				break;
			case 2 : /* pressure only */
				catseis(sectionp, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,4);

				break;
			case 3 : /* curl and div only */
				catseis(sectiondiv, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,5);
				catseis(sectioncurl, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,6);

				break;
			case 4 : /* everything */
				/*fprintf(FP," start merging, ntr= %d : \n",ntr_glob);
						fprintf(stdout,"Message from PE %d\n",MYID);*/
				catseis(sectionvx, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,1);
				catseis(sectionvy, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,2);
				catseis(sectionvz, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,3);
				catseis(sectionp, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,4);
				catseis(sectiondiv, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,5);
				catseis(sectioncurl, seismo_fulldata, recswitch, ntr_glob,ns);
				if (MYID==0) saveseis_glob(FP,seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,6);

				break;
			default :	break;

			}
			fprintf(FP, "\n\n");

		}


		/* output timing information (real times for update and exchange) */
		if (LOG)
			if (MYID==0) timing(time_v_update,time_s_update, time_s_exchange,time_v_exchange,time_timestep, ishot);



	} /* end of loop over shots */

	/* output of checkpoint files is temporarily disabled
	 * there are quite some variables that are not in use in sofi3D_acoustic, e.g. syy or szz!
	if (CHECKPTWRITE){
		if (MYID==0){
			time3=MPI_Wtime();
 			fprintf(FP," Saving wavefield to check-point file %s \n",CHECKPTFILE);	
		}

		save_checkpoint(-1, NX+2, -1, NY+2, -1, NZ+2, vx,vy,vz,sxx,syy,szz,sxx,syy,szz);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0){
			time4=MPI_Wtime();
      			fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
		}
	}*/

	/*de-allocation of memory */
	free_f3tensor(vx,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
	free_f3tensor(vy,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
	free_f3tensor(vz,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
	free_f3tensor(sxx,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

	free_f3tensor(rho,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(pi,0,NY+1,0,NX+1,0,NZ+1);

	if(ABS_TYPE==1){
		free_f3tensor(vx1,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(vy1,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(vz1,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(sxx1,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(syy1,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(szz1,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

		free_f3tensor(vx2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(vy2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(vz2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(sxx2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(syy2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(szz2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

		free_f3tensor(vx3,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(vy3,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(vz3,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(sxx3,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(syy3,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);
		free_f3tensor(szz3,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2);

		free_f3tensor(absorb_coeffx,1,NY,1,NX,1,NZ);
		free_f3tensor(absorb_coeffy,1,NY,1,NX,1,NZ);
		free_f3tensor(absorb_coeffz,1,NY,1,NX,1,NZ);
	}


	if(ABS_TYPE==2){
		free_f3tensor(absorb_coeff,1,NY,1,NX,1,NZ);
	}

	if (nsrc_loc>0){	
		free_matrix(signals,1,nsrc_loc,1,NT);
		free_matrix(srcpos_loc,1,6,1,nsrc_loc);
	}


	free_f3tensor(bufferlef_to_rig,1,NY,1,NZ,1,nf1);
	free_f3tensor(bufferrig_to_lef,1,NY,1,NZ,1,nf2);
	free_f3tensor(buffertop_to_bot,1,NX,1,NZ,1,nf1);
	free_f3tensor(bufferbot_to_top,1,NX,1,NZ,1,nf2);
	free_f3tensor(bufferfro_to_bac,1,NY,1,NX,1,nf1);
	free_f3tensor(bufferbac_to_fro,1,NY,1,NX,1,nf2);

	free_f3tensor(sbufferlef_to_rig,1,NY,1,NZ,1,nf2);
	free_f3tensor(sbufferrig_to_lef,1,NY,1,NZ,1,nf1);
	free_f3tensor(sbuffertop_to_bot,1,NX,1,NZ,1,nf2);
	free_f3tensor(sbufferbot_to_top,1,NX,1,NZ,1,nf1);
	free_f3tensor(sbufferfro_to_bac,1,NY,1,NX,1,nf2);
	free_f3tensor(sbufferbac_to_fro,1,NY,1,NX,1,nf1);


	/* free memory for global source positions */
	free_imatrix(recpos,1,3,1,ntr_glob);

	/* free memory for global source positions */
	free_matrix(srcpos,1,6,1,nsrc);


	if ((ntr>0) && (SEISMO)){	
		free_matrix(seismo_fulldata,1,ntr_glob,1,ns);
		free_imatrix(recpos_loc,1,3,1,ntr);

		switch (SEISMO){
		case 1 : /* particle velocities only */
			free_matrix(sectionvx,1,ntr,1,ns);
			free_matrix(sectionvy,1,ntr,1,ns);		
			free_matrix(sectionvz,1,ntr,1,ns);		
			break;	
		case 2 : /* pressure only */
			free_matrix(sectionp,1,ntr,1,ns);
			break;	
		case 3 : /* curl and div only */
			free_matrix(sectioncurl,1,ntr,1,ns);
			free_matrix(sectiondiv,1,ntr,1,ns);
			break;	
		case 4 : /* everything */
			free_matrix(sectionvx,1,ntr,1,ns);
			free_matrix(sectionvy,1,ntr,1,ns);
			free_matrix(sectionvz,1,ntr,1,ns);
			free_matrix(sectionp,1,ntr,1,ns);
			free_matrix(sectioncurl,1,ntr,1,ns);
			free_matrix(sectiondiv,1,ntr,1,ns);		
			break;
		}	

	}	


	/* free memory for source position definition */
	free_matrix(srcpos1,1,6,1,1);

	/* de-allocate buffer for messages */
	MPI_Buffer_detach(buff_addr,&buffsize);



	MPI_Barrier(MPI_COMM_WORLD);

	/* merge snapshot files created by the PEs into one file */
	/* if ((SNAP) && (MYID==0)) snapmerge(nsnap);*/


	/* free PML indices */
	free_ivector(xa,0,5);
	free_ivector(xb,0,5);
	free_ivector(ya,0,5);
	free_ivector(yb,0,5);
	free_ivector(za,0,5);
	free_ivector(zb,0,5);
	free_ivector(stype,1,nsrc);
	free_ivector(stype_loc,1,nsrc);

	/* free timing arrays */
	free_dvector(time_v_update,1,NT);
	free_dvector(time_s_update,1,NT);
	free_dvector(time_s_exchange,1,NT);
	free_dvector(time_v_exchange,1,NT);
	free_dvector(time_timestep,1,NT);


	if (MYID==0){
		fprintf(FP,"\n **Info from main (written by PE %d): \n",MYID);
		time4=MPI_Wtime();
		fprintf(FP," Total real time of program: %4.2f seconds.\n\n",time4-time1);
		fprintf(FP," ******************************************************\n");
		fprintf(FP," *********** SOFI3D_acoustic has finished *************\n");
		fprintf(FP," ******************************************************\n\n");
	}


	fclose(FP);

	MPI_Finalize();

	return 0;

}
