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
 * This is program SOFI3D.
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
	int ns, nt, nseismograms=0, nf1, nf2;
	int lsnap, nsnap=0, lsamp=0, nlsamp=0, buffsize;
	int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, ishot, nshots; /* removed variable "h", not in use*/

	double 	time1=0.0, time2=0.0, time3=0.0, time4=0.0;
	double * time_v_update, * time_s_update, * time_s_exchange,* time_v_exchange, * time_timestep;	
	int * xb, * yb, * zb, l;

	float  *** absorb_coeff=NULL;
	float  ***  sxy=NULL, ***  syz=NULL, ***  sxz=NULL;
	float  ***  sxx=NULL, ***  syy=NULL, ***  szz=NULL;
	float  ***  rxy=NULL, ***  ryz=NULL, ***  rxz=NULL;
	float  ***  rxx=NULL, ***  ryy=NULL, ***  rzz=NULL;
	float  ***  vx=NULL, ***  vy=NULL, ***  vz=NULL;
    
    /* Save old spatial derivations of velocity for Adam Bashforth */
    float *** vxyyx=NULL, *** vyzzy=NULL, *** vxzzx=NULL, *** vxxyyzz=NULL, *** vyyzz=NULL, *** vxxzz=NULL, *** vxxyy=NULL;
    float *** vxyyx_2=NULL, *** vyzzy_2=NULL, *** vxzzx_2=NULL, *** vxxyyzz_2=NULL, *** vyyzz_2=NULL, *** vxxzz_2=NULL, *** vxxyy_2=NULL;
    float *** vxyyx_3=NULL, *** vyzzy_3=NULL, *** vxzzx_3=NULL, *** vxxyyzz_3=NULL, *** vyyzz_3=NULL, *** vxxzz_3=NULL, *** vxxyy_3=NULL;
    float *** vxyyx_4=NULL, *** vyzzy_4=NULL, *** vxzzx_4=NULL, *** vxxyyzz_4=NULL, *** vyyzz_4=NULL, *** vxxzz_4=NULL, *** vxxyy_4=NULL;
    
    /* Save old derivation of the stress for Adam Bashforth */
    float *** svx=NULL,  *** svy=NULL, *** svz=NULL;
    float *** svx_2=NULL, *** svy_2=NULL, *** svz_2=NULL, *** svx_3=NULL, *** svy_3=NULL, *** svz_3=NULL;
    float *** svx_4=NULL, *** svy_4=NULL, *** svz_4=NULL;

    /* We need some pointers for the time shift for Adam Bashforth*/
    float *** shift_s1=NULL,*** shift_s2=NULL,*** shift_s3=NULL;
    float *** shift_v1=NULL,*** shift_v2=NULL,*** shift_v3=NULL,*** shift_v4=NULL,*** shift_v5=NULL,*** shift_v6=NULL,*** shift_v7=NULL;
    float *** shift_r1=NULL,*** shift_r2=NULL,*** shift_r3=NULL,*** shift_r4=NULL,*** shift_r5=NULL,*** shift_r6=NULL;
    
    /* We need some pointes for the memory variables for Adams Bashforth */
    float  ***  rxy_2=NULL, ***  ryz_2=NULL, ***  rxz_2=NULL;
    float  ***  rxx_2=NULL, ***  ryy_2=NULL, ***  rzz_2=NULL;
    float  ***  rxy_3=NULL, ***  ryz_3=NULL, ***  rxz_3=NULL;
    float  ***  rxx_3=NULL, ***  ryy_3=NULL, ***  rzz_3=NULL;
    float  ***  rxy_4=NULL, ***  ryz_4=NULL, ***  rxz_4=NULL;
    float  ***  rxx_4=NULL, ***  ryy_4=NULL, ***  rzz_4=NULL;
    
	float  ***  rho, ***  pi, ***  u;
	float  ***  taus=NULL, ***  taup=NULL, *eta=NULL;
	float  *** uipjp, *** ujpkp, *** uipkp, *** tausipjp=NULL, *** tausjpkp=NULL, *** tausipkp=NULL,*** rjp, *** rkp, *** rip;

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

	float * K_x=NULL, * alpha_prime_x=NULL, * a_x=NULL, * b_x=NULL, * K_x_half=NULL, * alpha_prime_x_half=NULL, * a_x_half=NULL, * b_x_half=NULL, * K_y=NULL, * alpha_prime_y=NULL, * a_y=NULL, * b_y=NULL, * K_y_half=NULL, * alpha_prime_y_half=NULL, * a_y_half=NULL, * b_y_half=NULL, * K_z=NULL, * alpha_prime_z=NULL, * a_z=NULL, * b_z=NULL, * K_z_half=NULL, * alpha_prime_z_half=NULL, * a_z_half=NULL, * b_z_half=NULL;
	float *** psi_sxx_x=NULL, *** psi_syy_y=NULL, *** psi_szz_z=NULL, *** psi_sxy_y=NULL, *** psi_sxy_x=NULL, *** psi_sxz_x=NULL, *** psi_sxz_z=NULL, *** psi_syz_y=NULL, *** psi_syz_z=NULL, *** psi_vxx=NULL, *** psi_vyy=NULL, *** psi_vzz=NULL, *** psi_vxy=NULL, *** psi_vxz=NULL, *** psi_vyx=NULL, *** psi_vyz=NULL, *** psi_vzx=NULL, *** psi_vzy=NULL;

	MPI_Request *req_send, *req_rec, *sreq_send, *sreq_rec;

	float memdyn, memmodel, memseismograms, membuffer, memcpml=0.0, memtotal;
	float amon=0.0, str=0.0, dip=0.0, rake=0.0;
	float fac1, fac2;
	char *buff_addr, ext[10];
	int * stype=NULL, *stype_loc=NULL;
	char buffer[STRING_SIZE], bufferstring[10];
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

	/* print program name, version, author etc to stdout*/
	if (MYID == 0) info(stdout);
	SOFI3DVERS=33; /* 3D isotropic elastic */

	/* PE 0 is reading the parameters from the input file, default is sofi3D.json */
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
			strcpy(FILEINP,"sofi3D.json");
			fprintf(stderr," Caution: input parameter filename set to default 'sofi3D.json'. \n\n");
		}
		FP=fopen(FILEINP,"r");

		//read json formated input file
		read_par_json(stdout, FILEINP);
		fclose(FP);
	} 
	/* PE 0 will broadcast the parameters to all others PEs */
	exchange_par(); 
    
	/* Print info on log-files to stdout */
	if (MYID == 0) note(stdout);

	/* open log-file (each PE is using different file) */
	/*	fp=stdout; */
	sprintf(ext,".%i",MYID);  
	strcat(LOG_FILE,ext);	

	/* nodes MYIDo writes logging info to LOG_FILE or stdout */	
	if (MYID==0) 
		switch (LOG){
		case 0 : FP=fopen("/dev/null","w"); /* no logging information will be output */
		break;
		case 1 : FP=stdout; /* logging information will be written to standard output */
		break;
		case 2 : if ((FP=fopen(LOG_FILE,"w"))==NULL) err(" Opening log-file failed.");
		/* logging information will be written to LOG_FILE */
		break;
		}

	/* all other nodes write logging info to LOG_FILE */		
	if (MYID>0) 
		if ((FP=fopen(LOG_FILE,"w"))==NULL) err(" Opening log-file failed.");

	fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);

	/* domain decomposition */
	initproc();

	/* set some time counters */
	NT=(int)ceil(TIME/DT); /* number of timesteps - replaces: NT=iround(TIME/DT); */
	TIME=(NT-1)*DT; /* TIME set to true time of the last time step */

	if (NDTSHIFT>NT) {
		ns=0;
	}
	else {
		ns=(int)ceil((float)(NT-NDTSHIFT)/(float)NDT); /* number of samples per trace - replaces buggy formula: ns=iround(NT-NDTSHIFT/NDT); */
	}
	lsnap=iround(TSNAP1/DT); /* first snapshot at this timestep */
	if (((ns<0) && (MYID==0)) && (SEISMO>0)) {
		/*fprintf(FP," \nSampling interval for seismogram output (DT) : %f \n\n",DT);
		fprintf(FP," \nSampling interval shift (NDTSHIFT) : %i , TIME : %f NT : %i NDT : %i \n\n",NDTSHIFT,TIME, NT,NDT);*/
		fprintf(FP," \nSampling rate for seismogram output (NDT) is out of limit : %i \n\n",ns);
		err(" Check Sampling rate for seismogram output (NDT)!");
	}

	/* output of parameters to stdout: */
	if (MYID==0) writepar(FP,ns);

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

	/* number of seismogram sections which have to be stored in core memory*/
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
	default : nseismograms=0;
	break;
	}

	/*allocate memory for dynamic, static and buffer arrays */
	fac1=(NZ+FDORDER)*(NY+FDORDER)*(NX+FDORDER);
    fac2=sizeof(float)*pow(2.0,-20.0);
    
    if (L>0){ /*viscoelastic case*/
        if(FDORDER_TIME==2){
            memdyn=15.0*fac1*fac2;
            memmodel=16.0*fac1*fac2;
        } else {
            if(FDORDER_TIME==3) {
                memdyn=57*fac1*fac2;
                memmodel=16.0*fac1*fac2;
            } else {
                memdyn=78*fac1*fac2;
                memmodel=16.0*fac1*fac2;
            }
        }
    } else{   /* elastic case*/
        if(FDORDER_TIME==2){
            memdyn=9.0*fac1*fac2;
            memmodel=10.0*fac1*fac2;
        } else {
            if(FDORDER_TIME==3) {
                memdyn=39*fac1*fac2;
                memmodel=10.0*fac1*fac2;
            } else {
                memdyn=49*fac1*fac2;
                memmodel=10.0*fac1*fac2;
            }
        }
    }

	memseismograms=nseismograms*ntr_glob*ns*fac2 + nseismograms*ntr*ns*fac2;
	membuffer=(2.0*(3.0*FDORDER/2-1)*(NY*NZ+NX*NZ+NY*NX)+2.0*(3.0*FDORDER/2-2)*(NY*NZ+NX*NZ+NY*NX))*fac2;
	membuffer=4.0*6.0*((NX*NZ)+(NY*NZ)+(NX*NY))*fac2;
	if (ABS_TYPE==1) memcpml=2.0*FW*6.0*(NY*NZ+NX*NZ+NY*NX)*fac2+24.0*2.0*FW*fac2;
	buffsize=(FDORDER)*4.0*6.0*(max((NX*NZ),max((NY*NZ),(NX*NY))))*sizeof(MPI_FLOAT);
	memtotal=memdyn+memmodel+memseismograms+membuffer+memcpml+(buffsize*pow(2.0,-20.0));


	if (MYID==0){
		fprintf(FP,"\n ----------------------------------------------------------------");
		fprintf(FP,"\n ------------------ MEMORY ALLOCATION --------------------------- \n");
		fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
		fprintf(FP," Size of local grids: NX=%d \t NY=%d \t NZ=%d \n",NX,NY,NZ);
		fprintf(FP," Each process is now trying to allocate memory for:\n");
		fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
		fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
		fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
		if (ABS_TYPE==1) fprintf(FP," CPML memory variables : \t %6.2f MB\n", memcpml);
		fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
		fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n\n", buffsize*pow(2.0,-20.0));
		fprintf(FP," Total memory required: \t %6.2f MB.\n", memtotal);
		fprintf(FP," Please note that the memory consumption is only a good estimate! ");
		fprintf(FP,"\n ---------------------------------------------------------------- \n\n");
	}

	/* allocate buffer for buffering messages */
	buff_addr=malloc(buffsize);
	if (!buff_addr) err("allocation failure for buffer for MPI_Bsend !");
	MPI_Buffer_attach(buff_addr,buffsize);


	/* allocation for request and status arrays */
	req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	sreq_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	sreq_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));

	/* allocation for timing arrays used for performance analysis */
	time_v_update=dvector(1,NT);
	time_s_update=dvector(1,NT);
	time_s_exchange=dvector(1,NT);
	time_v_exchange=dvector(1,NT);
	time_timestep=dvector(1,NT);

	l=1;
	if(ABS_TYPE==1 && FDORDER==2){l=2;}
    
	/* memory allocation for dynamic (wavefield) arrays */
	if(POS[2]==0){
		vx  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		vy  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		vz  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        
        if(FDORDER_TIME != 2){ /* Allocate memory for Adams Bashforth */
            vxyyx  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyzzy  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxzzx  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyyzz  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyyzz  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxzz  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyy  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            vxyyx_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyzzy_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxzzx_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyyzz_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyyzz_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxzz_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyy_2  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            vxyyx_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyzzy_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxzzx_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyyzz_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyyzz_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxzz_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyy_3  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            svx=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svy=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svz=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            svx_2=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svy_2=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svz_2=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            svx_3=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svy_3=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svz_3=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            if(FDORDER_TIME==4){
                svx_4=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                svy_4=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                svz_4=f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                
                vxyyx_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vyzzy_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxzzx_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxxyyzz_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vyyzz_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxxzz_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxxyy_4  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            }
        }
        
		sxy =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		syz =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	}

    if(POS[2]>0){
        vx  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        vy  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        vz  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        
        if(FDORDER_TIME != 2){ /* Allocate memory for Adams Bashforth */
            vxyyx  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyzzy  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxzzx  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyyzz  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyyzz  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxzz  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyy  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            vxyyx_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyzzy_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxzzx_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyyzz_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyyzz_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxzz_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyy_2  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            vxyyx_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyzzy_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxzzx_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyyzz_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vyyzz_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxzz_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            vxxyy_3  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            svx=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svy=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svz=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            svx_2=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svy_2=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svz_2=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            svx_3=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svy_3=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            svz_3=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            if(FDORDER_TIME==4){
                svx_4=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                svy_4=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                svz_4=f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                
                vxyyx_4  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vyzzy_4  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxzzx_4  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxxyyzz_4  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vyyzz_4  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxxzz_4  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                vxxyy_4  = f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            }
        }
        
        sxy =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        syz =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
    }

	sxz =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	sxx =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	syy =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	szz =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);

	xb=ivector(0,1);
	yb=ivector(0,1);
	zb=ivector(0,1);

	if (L){ /* no allocation in case of purely elastic simulation */
		/* memory allocation for dynamic (model) arrays */
		rxy =  f3tensor(1,NY,1,NX,1,NZ);
		ryz =  f3tensor(1,NY,1,NX,1,NZ);
		rxz =  f3tensor(1,NY,1,NX,1,NZ);
		rxx =  f3tensor(1,NY,1,NX,1,NZ);
		ryy =  f3tensor(1,NY,1,NX,1,NZ);
		rzz =  f3tensor(1,NY,1,NX,1,NZ);
        
        if(FDORDER_TIME != 2){ /* Allocate memory for Adams Bashforth */
            rxy_2 =  f3tensor(1,NY,1,NX,1,NZ);
            ryz_2 =  f3tensor(1,NY,1,NX,1,NZ);
            rxz_2 =  f3tensor(1,NY,1,NX,1,NZ);
            rxx_2 =  f3tensor(1,NY,1,NX,1,NZ);
            ryy_2 =  f3tensor(1,NY,1,NX,1,NZ);
            rzz_2 =  f3tensor(1,NY,1,NX,1,NZ);
            
            rxy_3 =  f3tensor(1,NY,1,NX,1,NZ);
            ryz_3 =  f3tensor(1,NY,1,NX,1,NZ);
            rxz_3 =  f3tensor(1,NY,1,NX,1,NZ);
            rxx_3 =  f3tensor(1,NY,1,NX,1,NZ);
            ryy_3 =  f3tensor(1,NY,1,NX,1,NZ);
            rzz_3 =  f3tensor(1,NY,1,NX,1,NZ);
            if(FDORDER_TIME==4){
                rxy_4 =  f3tensor(1,NY,1,NX,1,NZ);
                ryz_4 =  f3tensor(1,NY,1,NX,1,NZ);
                rxz_4 =  f3tensor(1,NY,1,NX,1,NZ);
                rxx_4 =  f3tensor(1,NY,1,NX,1,NZ);
                ryy_4 =  f3tensor(1,NY,1,NX,1,NZ);
                rzz_4 =  f3tensor(1,NY,1,NX,1,NZ);
            }
        }
        
		/* memory allocation for static (model) arrays */
		taus=  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		taup=  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		tausipjp=f3tensor(1,NY,1,NX,1,NZ);
		tausjpkp=f3tensor(1,NY,1,NX,1,NZ);
		tausipkp=f3tensor(1,NY,1,NX,1,NZ);
		eta =  vector(1,L);
	}

	/* memory allocation for static (model) arrays */
	rho =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	pi  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	u   =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);

	absorb_coeff=  f3tensor(1,NY,1,NX,1,NZ);

	/* averaged material parameters */
	uipjp=f3tensor(1,NY,1,NX,1,NZ);
	ujpkp=f3tensor(1,NY,1,NX,1,NZ);
	uipkp=f3tensor(1,NY,1,NX,1,NZ);
	rjp=f3tensor(1,NY,1,NX,1,NZ);
	rkp=f3tensor(1,NY,1,NX,1,NZ);
	rip=f3tensor(1,NY,1,NX,1,NZ);


	/* memory allocation for CPML variables*/
	if(ABS_TYPE==1){

		K_x = vector(1,2*FW);
		alpha_prime_x = vector(1,2*FW);
		a_x = vector(1,2*FW);
		b_x = vector(1,2*FW);
		K_x_half = vector(1,2*FW);
		alpha_prime_x_half = vector(1,2*FW);
		a_x_half = vector(1,2*FW);
		b_x_half = vector(1,2*FW);

		K_y = vector(1,2*FW);
		alpha_prime_y = vector(1,2*FW);
		a_y = vector(1,2*FW);
		b_y = vector(1,2*FW);
		K_y_half = vector(1,2*FW);
		alpha_prime_y_half = vector(1,2*FW);
		a_y_half = vector(1,2*FW);
		b_y_half = vector(1,2*FW);

		K_z = vector(1,2*FW);
		alpha_prime_z = vector(1,2*FW);
		a_z = vector(1,2*FW);
		b_z = vector(1,2*FW);
		K_z_half = vector(1,2*FW);
		alpha_prime_z_half = vector(1,2*FW);
		a_z_half = vector(1,2*FW);
		b_z_half = vector(1,2*FW);


		psi_sxx_x =  f3tensor(1,NY,1,2*FW,1,NZ);
		psi_sxy_x =  f3tensor(1,NY,1,2*FW,1,NZ);
		psi_sxz_x =  f3tensor(1,NY,1,2*FW,1,NZ);
		psi_syy_y =  f3tensor(1,2*FW,1,NX,1,NZ);
		psi_sxy_y =  f3tensor(1,2*FW,1,NX,1,NZ);
		psi_syz_y =  f3tensor(1,2*FW,1,NX,1,NZ);
		psi_szz_z =  f3tensor(1,NY,1,NX,1,2*FW);
		psi_sxz_z =  f3tensor(1,NY,1,NX,1,2*FW);
		psi_syz_z =  f3tensor(1,NY,1,NX,1,2*FW);


		psi_vxx   =  f3tensor(1,NY,1,2*FW,1,NZ);
		psi_vyy   =  f3tensor(1,2*FW,1,NX,1,NZ);
		psi_vzz   =  f3tensor(1,NY,1,NX,1,2*FW);
		psi_vxy   =  f3tensor(1,2*FW,1,NX,1,NZ);
		psi_vxz   =  f3tensor(1,NY,1,NX,1,2*FW);
		psi_vyx   =  f3tensor(1,NY,1,2*FW,1,NZ);
		psi_vyz   =  f3tensor(1,NY,1,NX,1,2*FW);
		psi_vzx   =  f3tensor(1,NY,1,2*FW,1,NZ);
		psi_vzy   =  f3tensor(1,2*FW,1,NX,1,NZ);

	}


	/* memory allocation for buffer arrays in which the wavefield information which is exchanged between neighboring PEs is stored */

	/* number of wavefield parameters that need to be exchanged - see exchange_v.c */
	nf1=(3*FDORDER/2)-1;
	nf2=nf1-1;

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

	/* allocate buffer for seismogram output, merged seismogram section of all PEs */
	if (SEISMO) seismo_fulldata=fmatrix(1,ntr_glob,1,ns);

	/* allocate buffer for seismogram output, seismogram section of each PE */
	if (ntr>0){
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


	if (MYID==0) 
		fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);


	/* memory for source position definition */
	srcpos1=fmatrix(1,6,1,1);

	/* Reading source positions from SOURCE_FILE */ 	
	fprintf(FP,"\n ------------------ READING SOURCE PARAMETERS ------------------- \n");
	switch (SRCREC) {
	case 0:
		if (MYID==0) err("SRCREC parameter is invalid (SRCREC!=1)! No source parameters specified!");
		break;
	case 1:
		if (MYID==0) {
			fprintf(FP,"\n Reading source parameters from file: %s (SOFI3D source format)\n",SOURCE_FILE);

			if ((fpsrc=fopen(SOURCE_FILE,"r"))==NULL) err(" Source file could not be opened !");
			while(fgets(buffer, STRING_SIZE, fpsrc))
			{
				sscanf(buffer,"%s",bufferstring);
				/* checks if the line contains a '%'character which indicates a comment line,
					and if the reading of a string was successful, which is not the case for an empty line*/
				if ((strchr(buffer,'#')==0) && (sscanf(buffer,"%s",bufferstring)==1)) ++(nsrc);
			}
			rewind(fpsrc);
			if ((nsrc)==0) fprintf(FP,"\n WARNING: Could not determine number of sources parameter sets in input file. Assuming %d.\n",(nsrc=0));
			else fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,nsrc);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);

		stype = ivector(1,nsrc);
		srcpos=sources(fpsrc,&nsrc,stype);

		/*originally, SOURCE_TYPE=stype is defined in the source file, if not, SOURCE_TYPE is taken from the input file */
		/*if (stype==NULL) printf("PE%d: Source type(s) undefined?! \n",MYID);*/

		break;
	case 2:
		if ((PLANE_WAVE_DEPTH>0)) {
			if (MYID==0) {
				/*stype=(int *)malloc(nsrc*sizeof(int));*/ /* for unknown reasons, the pointer does not point to memory that has been allocated by a subroutine this way */

				/*determining the number of sources in the specified plane normal/tilted to the surface/upper model boundary*/
				nsrc=(NXG-2*FW+1)*(NZG-2*FW+1);
				/*fprintf(FP,"\n nsrc= %i with NGX=%i, NYG=%i and FW=%i. \n",nsrc,NXG,NYG,FW);*/
			}

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);

			stype = ivector(1,nsrc);
			srcpos=pwsources(&nsrc,stype);

		}
		else {
			err("SRCREC parameter specifies PLANE_WAVE excitation, but PLANE_WAVE_DEPTH<=0!");
		}
		break;
		/*more source file formats or source file options can be implemented here*/
	default: err("SRCREC parameter is invalid (SRCREC!=1 or SRCREC!=2)! No source parameters specified!");
	break;
	}

	/* create model grids */
	fprintf(FP,"\n ------------------ MODEL CREATION AND OUTPUT-------------------- \n");
	if (READMOD) readmod(rho,pi,u,taus,taup,eta);
	else {
		if (L==0) model_elastic(rho,pi,u,taus,taup,eta); /* elastic modeling, L is specified in input file*/
        else {
            model_visco(rho,pi,u,taus,taup,eta); /* viscoelastic modeling, L is specified in input file*/
        }
	}
	
	
	if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;		
	/*printf("\n ------------------ checkfd by MYID %i -------------------- \n", MYID);*/
	/* check if the FD run will be stable and free of numerical dispersion */
	checkfd(FP,rho,pi,u,taus,taup,eta,srcpos,nsrc,recpos,ntr_glob);

	/* calculate 3-D array for exponential damping of reflections
           at the edges of the numerical mesh (PML-boundary)*/
	/*if(ABS_TYPE==1){   
	  absorb_PML(absorb_coeffx, absorb_coeffy, absorb_coeffz);
        }*/

	/* calculate damping coefficients for CPML boundary*/
	if(ABS_TYPE==1){
		CPML_coeff(K_x,alpha_prime_x,a_x,b_x,K_x_half,alpha_prime_x_half,a_x_half,b_x_half,K_y,alpha_prime_y,a_y,b_y,K_y_half,alpha_prime_y_half,a_y_half,b_y_half,K_z,alpha_prime_z,a_z,b_z,K_z_half,alpha_prime_z_half,a_z_half,b_z_half);
	}

	/* calculate 3-D array for exponential damping of reflections
	   at the edges of the numerical mesh */
	if(ABS_TYPE==2){   
		absorb(absorb_coeff);
	}

	/* For the calculation of the material parameters beteween gridpoints
	   the have to be averaged. For this, values lying at 0 and NX+1,
	for example, are required on the local grid. These are now copied from the
	neighbouring grids */
	matcopy(rho,pi,u,taus,taup);

	/* spatial averaging of material parameters, i.e. Tau for S-waves, shear modulus, and density */
	av_mat(rho,pi,u,taus,taup,uipjp,ujpkp,uipkp,tausipjp,tausjpkp,tausipkp,rjp,rkp,rip);
    
    
	MPI_Barrier(MPI_COMM_WORLD);

	if (CHECKPTREAD){
		if (MYID==0){
			time3=MPI_Wtime();
			fprintf(FP," Reading wavefield from check-point file %s \n",CHECKPTFILE);
		}

		read_checkpoint(-1, NX+2, -1, NY+2, -1, NZ+2, vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,
			psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z,
			psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0){
			time4=MPI_Wtime();
			fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
		}
	}



	/* comunication initialisation for persistent communication */
	/*comm_ini(bufferlef_to_rig, bufferrig_to_lef,
	    buffertop_to_bot, bufferbot_to_top, bufferfro_to_bac,
	    bufferbac_to_fro, req_send, req_rec);

	comm_ini_s(sbufferlef_to_rig, sbufferrig_to_lef,
	    sbuffertop_to_bot, sbufferbot_to_top, sbufferfro_to_bac,
	    sbufferbac_to_fro, sreq_send, sreq_rec);*/

	/* initialisation of PML and ABS domain */
	if(ABS_TYPE==1){    
		CPML_ini_elastic(xb,yb,zb);
	}

	if(ABS_TYPE==2){    
		xb[0]=1; xb[1]=NX;
		yb[0]=1; yb[1]=NY;
		zb[0]=1; zb[1]=NZ;
	}

	if (MYID==0){
		time2=MPI_Wtime();
		fprintf(FP,"\n\n\n **************************************************\n");
		fprintf(FP," *********** STARTING TIME STEPPING ***************\n");
		fprintf(FP," **************************************************\n");
		fprintf(FP," real time before starting time loop: %4.2f s.\n",time2-time1);
	}


	for (ishot=1;ishot<=nshots;ishot++){

		fprintf(FP,"\n MYID=%d *****  Starting simulation for shot %d of %d  ********** \n",MYID,ishot,nshots);
		for (nt=1;nt<=6;nt++) srcpos1[nt][1]=srcpos[nt][ishot];
		if (RUN_MULTIPLE_SHOTS){
			/* find this single source positions on subdomains */
			if (nsrc_loc>0) free_matrix(srcpos_loc,1,6,1,1);
			if (stype_loc==NULL) stype_loc = ivector(1,nsrc);;
			srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1, stype_loc, stype);
		}
		else {
			/* Distribute multiple source positions on subdomains */
			if (stype_loc==NULL) stype_loc = ivector(1,nsrc);
			srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc, stype_loc, stype);
		}

		/* calculate wavelet for each source point */
		signals=wavelet(srcpos_loc,nsrc_loc);

		/* output of calculated wavelet for each source point */

		if ((OUTSOURCEWAVELET !=0 ) && (nsrc_loc>0)) {
			char  source_signal_file[STRING_SIZE];
			sprintf(source_signal_file,"source_signal.MYID%d.shot%d.su",MYID,ishot);
			fprintf(FP,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
			output_source_signal(fopen(source_signal_file,"w"), signals, NT, 1);
		}
        
        /* initialize wavefield with zero */
        if ((L==1) && (ABS_TYPE==2) && (CHECKPTREAD==0)){
            zero(1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,
                 vxyyx,vyzzy,vxzzx,vxxyyzz,vyyzz,vxxzz,vxxyy,vxyyx_2,vyzzy_2,vxzzx_2,vxxyyzz_2,vyyzz_2,vxxzz_2,vxxyy_2,vxyyx_3,
                 vyzzy_3,vxzzx_3,vxxyyzz_3,vyyzz_3,vxxzz_3,vxxyy_3,vxyyx_4,vyzzy_4,vxzzx_4,
                 vxxyyzz_4,vyyzz_4,vxxzz_4,vxxyy_4,svx,svy,svz,svx_2,svy_2,svz_2,svx_3,svy_3,svz_3,svx_4,svy_4,svz_4,rxx,ryy,rzz,rxy,ryz,rxz,rxx_2,ryy_2,rzz_2,rxy_2,ryz_2,rxz_2,rxx_3,ryy_3,rzz_3,rxy_3,ryz_3,rxz_3,rxx_4,ryy_4,rzz_4,rxy_4,ryz_4,rxz_4);
        }
        
		if((L==0) && (ABS_TYPE==2) && (CHECKPTREAD==0)){
            zero_elastic(1-FDORDER/2,NX+FDORDER/2,1-FDORDER/2,NY+FDORDER/2,1-FDORDER/2,NZ+FDORDER/2,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,
                         vxyyx,vyzzy,vxzzx,vxxyyzz,vyyzz,vxxzz,vxxyy,vxyyx_2,vyzzy_2,vxzzx_2,vxxyyzz_2,vyyzz_2,vxxzz_2,vxxyy_2,vxyyx_3,
                         vyzzy_3,vxzzx_3,vxxyyzz_3,vyyzz_3,vxxzz_3,vxxyy_3,vxyyx_4,vyzzy_4,vxzzx_4,
                         vxxyyzz_4,vyyzz_4,vxxzz_4,vxxyy_4,svx,svy,svz,svx_2,svy_2,svz_2,svx_3,svy_3,svz_3,svx_4,svy_4,svz_4);
		}
		if((ABS_TYPE==1) && (CHECKPTREAD==0)){
			zero_elastic_CPML(NX,NY,NZ,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz,rxx_2,ryy_2,rzz_2,rxy_2,ryz_2,rxz_2,rxx_3,ryy_3,rzz_3,rxy_3,ryz_3,rxz_3,rxx_4,ryy_4,rzz_4,rxy_4,ryz_4,rxz_4);

		}
        
		/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
		/* start of loop over time steps */

		lsamp=NDTSHIFT+1;
		nlsamp=1;
		for (nt=1;nt<=NT;nt++){

			time_v_update[nt]=0.0;
			time_s_update[nt]=0.0;
			time_v_exchange[nt]=0.0;
			time_s_exchange[nt]	=0.0;

			/* Check if simulation is still stable */
			if (isnan(vy[NY/2][NX/2][NZ/2])) err(" Simulation is unstable !"); /* maybe just breaking the loop would be better */

			if (LOG)
				if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
					fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
					time2=MPI_Wtime();
				}

			/* update of particle velocities */
			time_v_update[nt]=update_v(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rho,rjp, rkp,rip,srcpos_loc,signals,nsrc_loc,absorb_coeff,stype_loc,svx,svy,svz,svx_2,svy_2,svz_2,svx_3,svy_3,svz_3,svx_4,svy_4,svz_4);

			if(ABS_TYPE==1){
				update_v_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rho,rjp,rkp,rip,srcpos_loc,signals,nsrc_loc,absorb_coeff,stype_loc,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half, K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z);}
			;
            
            /* Shift spartial derivations of the stress one time step back */
            if (FDORDER_TIME==4){
                shift_s1=svx_4;svx_4=svx_3;svx_3=svx_2;svx_2=svx;svx=shift_s1;
                shift_s2=svy_4;svy_4=svy_3;svy_3=svy_2;svy_2=svy;svy=shift_s2;
                shift_s3=svz_4;svz_4=svz_3;svz_3=svz_2;svz_2=svz;svz=shift_s3;
            }
            if (FDORDER_TIME==3){
                shift_s1=svx_3; svx_3=svx_2;svx_2=svx;svx=shift_s1;
                shift_s2=svy_3;svy_3=svy_2;svy_2=svy;svy=shift_s2;
                shift_s3=svz_3;svz_3=svz_2; svz_2=svz;svz=shift_s3;
            }

			/* exchange values of particle velocities at grid boundaries between PEs */

			time_v_exchange[nt]=exchange_v(nt,vx,vy,vz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top,
					bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec);

			/* update of components of stress tensor */

            
            /* update NON PML boundaries */
            if (L>0){
                
                time_s_update[nt]=update_s(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,
                                           pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta,vxyyx,vyzzy,vxzzx,vxxyyzz,vyyzz,vxxzz,vxxyy,vxyyx_2,vyzzy_2,
                                           vxzzx_2,vxxyyzz_2,vyyzz_2,vxxzz_2,vxxyy_2,vxyyx_3,vyzzy_3,vxzzx_3,vxxyyzz_3,vyyzz_3,vxxzz_3,vxxyy_3,vxyyx_4,vyzzy_4,vxzzx_4,
                                           vxxyyzz_4,vyyzz_4,vxxzz_4,vxxyy_4,rxx_2,ryy_2,rzz_2,rxy_2,ryz_2,rxz_2,rxx_3,ryy_3,rzz_3,rxy_3,ryz_3,rxz_3,rxx_4,ryy_4,rzz_4,rxy_4,ryz_4,rxz_4);
                if(ABS_TYPE==1) update_s_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,
                                              rxy,ryz,rxz,pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta,K_x,a_x,b_x,K_x_half,a_x_half,
                                              b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,
                                              psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
            } else{
                time_s_update[nt]=update_s_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,
                                                   pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta,vxyyx,vyzzy,vxzzx,vxxyyzz,vyyzz,vxxzz,vxxyy,vxyyx_2,vyzzy_2,
                                                   vxzzx_2,vxxyyzz_2,vyyzz_2,vxxzz_2,vxxyy_2,vxyyx_3,vyzzy_3,vxzzx_3,vxxyyzz_3,vyyzz_3,vxxzz_3,vxxyy_3,vxyyx_4,vyzzy_4,vxzzx_4,
                                                   vxxyyzz_4,vyyzz_4,vxxzz_4,vxxyy_4);
                if(ABS_TYPE==1) update_s_CPML_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,
                                                      pi,u,uipjp,ujpkp,uipkp,K_x,a_x,b_x,K_x_half,a_x_half,
                                                      b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,
                                                      psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
            }
           
            /* Shift spartial derivations from the velocity one time step back */
            if (FDORDER_TIME==4){
                shift_v1=vxyyx_4;vxyyx_4=vxyyx_3;vxyyx_3=vxyyx_2;vxyyx_2=vxyyx;vxyyx=shift_v1;
                shift_v2=vyzzy_4;vyzzy_4=vyzzy_3;vyzzy_3=vyzzy_2;vyzzy_2=vyzzy;vyzzy=shift_v2;
                shift_v3=vxzzx_4;vxzzx_4=vxzzx_3;vxzzx_3=vxzzx_2;vxzzx_2=vxzzx;vxzzx=shift_v3;
                shift_v4=vxxyyzz_4;vxxyyzz_4=vxxyyzz_3;vxxyyzz_3=vxxyyzz_2;vxxyyzz_2=vxxyyzz;vxxyyzz=shift_v4;
                shift_v5=vyyzz_4;vyyzz_4=vyyzz_3;vyyzz_3=vyyzz_2;vyyzz_2=vyyzz;vyyzz=shift_v5;
                shift_v6=vxxzz_4;vxxzz_4=vxxzz_3;vxxzz_3=vxxzz_2; vxxzz_2=vxxzz;vxxzz=shift_v6;
                shift_v7=vxxyy_4;vxxyy_4=vxxyy_3;vxxyy_3=vxxyy_2;vxxyy_2=vxxyy;vxxyy=shift_v7;
                if (L==1) {
                    shift_r1=rxy_4;rxy_4=rxy_3; rxy_3=rxy_2;rxy_2=rxy;rxy=shift_r1;
                    shift_r2=ryz_4;ryz_4=ryz_3;ryz_3=ryz_2;ryz_2=ryz;ryz=shift_r2;
                    shift_r3=rxz_4;rxz_4=rxz_3;rxz_3=rxz_2;rxz_2=rxz;rxz=shift_r3;
                    shift_r4=rxx_4;rxx_4=rxx_3;rxx_3=rxx_2;rxx_2=rxx;rxx=shift_r4;
                    shift_r5=ryy_4;ryy_4=ryy_3;ryy_3=ryy_2;ryy_2=ryy;ryy=shift_r5;
                    shift_r6=rzz_4;rzz_4=rzz_3;rzz_3=rzz_2;rzz_2=rzz;rzz=shift_r6;
                }
                
            }
            if (FDORDER_TIME==3){
                shift_v1=vxyyx_3;vxyyx_3=vxyyx_2;vxyyx_2=vxyyx;vxyyx=shift_v1;
                shift_v2=vyzzy_3;vyzzy_3=vyzzy_2;vyzzy_2=vyzzy;vyzzy=shift_v2;
                shift_v3=vxzzx_3;vxzzx_3=vxzzx_2;vxzzx_2=vxzzx;vxzzx=shift_v3;
                shift_v4=vxxyyzz_3;vxxyyzz_3=vxxyyzz_2;vxxyyzz_2=vxxyyzz;vxxyyzz=shift_v4;
                shift_v5=vyyzz_3;vyyzz_3=vyyzz_2;vyyzz_2=vyyzz;vyyzz=shift_v5;
                shift_v6=vxxzz_3;vxxzz_3=vxxzz_2;vxxzz_2=vxxzz;vxxzz=shift_v6;
                shift_v7=vxxyy_3;vxxyy_3=vxxyy_2;vxxyy_2=vxxyy;vxxyy=shift_v7;
                if (L==1) {
                    shift_r1=rxy_3;rxy_3=rxy_2;rxy_2=rxy;rxy=shift_r1;
                    shift_r2=ryz_3;ryz_3=ryz_2;ryz_2=ryz;ryz=shift_r2;
                    shift_r3=rxz_3;rxz_3=rxz_2;rxz_2=rxz;rxz=shift_r3;
                    shift_r4=rxx_3;rxx_3=rxx_2;rxx_2=rxx;rxx=shift_r4;
                    shift_r5=ryy_3;ryy_3=ryy_2;ryy_2=ryy;ryy=shift_r5;
                    shift_r6=rzz_3;rzz_3=rzz_2;rzz_2=rzz;rzz=shift_r6;
                }
            }

			/* explosive source */
			if (!CHECKPTREAD){
				psource(nt,sxx,syy,szz,srcpos_loc,signals,nsrc_loc,stype_loc);
				/*
				 * eqsource is a implementation of moment tensor points sources. */
				 eqsource(nt,sxx,syy,szz,sxy, syz, sxz, srcpos_loc,signals,nsrc_loc,stype_loc, amon, str, dip, rake);
			}
			
			/* stress free surface ? */
			if ((FREE_SURF) && (POS[2]==0)){
				if (L) surface(1,u,pi,taus,taup,eta,sxx,syy,szz,sxy,syz,rxx,ryy,rzz,vx,vy,vz,K_x,a_x,b_x,
                                                      K_z,a_z,b_z,psi_vxx, psi_vzz);
				else surface_elastic(1,u,pi,sxx,syy,szz,sxy,syz,vx,vy,vz,K_x,a_x,b_x,
                                                      K_z,a_z,b_z,psi_vxx, psi_vzz);
			}

			/* exchange values of stress at boundaries between PEs */
			time_s_exchange[nt]=exchange_s(nt,sxx,syy,szz,sxy,syz,sxz,sbufferlef_to_rig,sbufferrig_to_lef,sbuffertop_to_bot,sbufferbot_to_top,sbufferfro_to_bac,sbufferbac_to_fro, sreq_send, sreq_rec);


			/* store amplitudes at receivers in e.g. sectionvx, sectionvz, sectiondiv, ...*/
			if ((SEISMO) && (ntr>0) && (nt==lsamp)){
				seismo(nlsamp,ntr,recpos_loc,sectionvx,sectionvy,sectionvz,
						sectiondiv,sectioncurl,sectionp,vx,vy,vz,sxx,syy,szz,pi,u);
				nlsamp++;
				lsamp+=NDT;
			}


			/* save snapshot in file */
			if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
				snap(FP,nt,++nsnap,SNAP_FORMAT,SNAP,vx,vy,vz,sxx,syy,szz,u,pi,
						IDX,IDY,IDZ,1,1,1,NX,NY,NZ);
				lsnap=lsnap+iround(TSNAPINC/DT);
			}

			if (LOG)
				if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
					time3=MPI_Wtime();
					time_timestep[nt]=(time3-time2);
					fprintf(FP," total real time for timestep %d : \t\t %4.2f s.\n",nt,time3-time2);
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

	if (CHECKPTWRITE){
		if (MYID==0){
			time3=MPI_Wtime();
			fprintf(FP," Saving wavefield to check-point file %s \n",CHECKPTFILE);
		}

		save_checkpoint(-1, NX+2, -1, NY+2, -1, NZ+2, vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,
				psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z,
				psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0){
			time4=MPI_Wtime();
			fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
		}
	}

	l=1;
	if(ABS_TYPE==1 && FDORDER==2){l=2;}

	/*de-allocation of memory */
	if(POS[2]==0){
		free_f3tensor(vx,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(vy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(vz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(sxy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(syz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        
        if(FDORDER_TIME != 2){
            free_f3tensor(vxyyx,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyzzy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxzzx,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyyzz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyyzz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxzz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);

            free_f3tensor(vxyyx_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyzzy_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxzzx_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyyzz_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyyzz_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxzz_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyy_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(vxyyx_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyzzy_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxzzx_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyyzz_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyyzz_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxzz_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyy_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(svx,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(svx_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svy_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svz_2,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(svx_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svy_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svz_3,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            
            if(FDORDER_TIME==4){
                free_f3tensor(vxyyx_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vyzzy_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxzzx_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxxyyzz_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vyyzz_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxxzz_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxxyy_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                
                free_f3tensor(svx_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(svy_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(svz_4,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            }

        }

	}

	if(POS[2]>0){
		free_f3tensor(vx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(vy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(vz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(sxy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f3tensor(syz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        
        if(FDORDER_TIME != 2){
            free_f3tensor(vxyyx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyzzy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxzzx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyyzz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyyzz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxzz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(vxyyx_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyzzy_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxzzx_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyyzz_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyyzz_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxzz_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyy_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(vxyyx_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyzzy_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxzzx_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyyzz_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vyyzz_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxzz_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(vxxyy_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(svx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(svx_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svy_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svz_2,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            free_f3tensor(svx_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svy_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            free_f3tensor(svz_3,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            
            
            if(FDORDER_TIME==4){
                free_f3tensor(vxyyx_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vyzzy_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxzzx_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxxyyzz_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vyyzz_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxxzz_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(vxxyy_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                
                free_f3tensor(svx_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(svy_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
                free_f3tensor(svz_4,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
            }
            
        }
        
        
	}

	free_f3tensor(sxz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(sxx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(syy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(szz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);

	if(ABS_TYPE==1){

		free_vector(K_x,1,2*FW);
		free_vector(alpha_prime_x,1,2*FW);
		free_vector(a_x,1,2*FW);
		free_vector(b_x,1,2*FW);
		free_vector(K_x_half,1,2*FW);
		free_vector(alpha_prime_x_half,1,2*FW);
		free_vector(a_x_half,1,2*FW);
		free_vector(b_x_half,1,2*FW);

		free_vector(K_y,1,2*FW);
		free_vector(alpha_prime_y,1,2*FW);
		free_vector(a_y,1,2*FW);
		free_vector(b_y,1,2*FW);
		free_vector(K_y_half,1,2*FW);
		free_vector(alpha_prime_y_half,1,2*FW);
		free_vector(a_y_half,1,2*FW);
		free_vector(b_y_half,1,2*FW);

		free_vector(K_z,1,2*FW);
		free_vector(alpha_prime_z,1,2*FW);
		free_vector(a_z,1,2*FW);
		free_vector(b_z,1,2*FW);
		free_vector(K_z_half,1,2*FW);
		free_vector(alpha_prime_z_half,1,2*FW);
		free_vector(a_z_half,1,2*FW);
		free_vector(b_z_half,1,2*FW);

		free_f3tensor(psi_sxx_x,1,NY,1,2*FW,1,NZ);
		free_f3tensor(psi_syy_y,1,2*FW,1,NX,1,NZ);
		free_f3tensor(psi_szz_z,1,NY,1,NX,1,2*FW);
		free_f3tensor(psi_sxy_x,1,NY,1,2*FW,1,NZ);
		free_f3tensor(psi_sxy_y,1,2*FW,1,NX,1,NZ);
		free_f3tensor(psi_sxz_x,1,NY,1,2*FW,1,NZ);
		free_f3tensor(psi_sxz_z,1,NY,1,NX,1,2*FW);
		free_f3tensor(psi_syz_y,1,2*FW,1,NX,1,NZ);
		free_f3tensor(psi_syz_z,1,NY,1,NX,1,2*FW);

		free_f3tensor(psi_vxx,1,NY,1,2*FW,1,NZ);
		free_f3tensor(psi_vyy,1,2*FW,1,NX,1,NZ);
		free_f3tensor(psi_vzz,1,NY,1,NX,1,2*FW);
		free_f3tensor(psi_vxy,1,2*FW,1,NX,1,NZ);
		free_f3tensor(psi_vyx,1,NY,1,2*FW,1,NZ);
		free_f3tensor(psi_vxz,1,NY,1,NX,1,2*FW);
		free_f3tensor(psi_vzx,1,NY,1,2*FW,1,NZ);
		free_f3tensor(psi_vyz,1,NY,1,NX,1,2*FW);
		free_f3tensor(psi_vzy,1,2*FW,1,NX,1,NZ);
	}

	if (L) {
        free_f3tensor(rxx,1,NY,1,NX,1,NZ);
        free_f3tensor(ryy,1,NY,1,NX,1,NZ);
        free_f3tensor(rzz,1,NY,1,NX,1,NZ);
        free_f3tensor(rxy,1,NY,1,NX,1,NZ);
        free_f3tensor(ryz,1,NY,1,NX,1,NZ);
        free_f3tensor(rxz,1,NY,1,NX,1,NZ);
        if (FDORDER_TIME>2) {
            free_f3tensor(rxx_2,1,NY,1,NX,1,NZ);
            free_f3tensor(ryy_2,1,NY,1,NX,1,NZ);
            free_f3tensor(rzz_2,1,NY,1,NX,1,NZ);
            free_f3tensor(rxy_2,1,NY,1,NX,1,NZ);
            free_f3tensor(ryz_2,1,NY,1,NX,1,NZ);
            free_f3tensor(rxz_2,1,NY,1,NX,1,NZ);
            free_f3tensor(rxx_3,1,NY,1,NX,1,NZ);
            free_f3tensor(ryy_3,1,NY,1,NX,1,NZ);
            free_f3tensor(rzz_3,1,NY,1,NX,1,NZ);
            free_f3tensor(rxy_3,1,NY,1,NX,1,NZ);
            free_f3tensor(ryz_3,1,NY,1,NX,1,NZ);
            free_f3tensor(rxz_3,1,NY,1,NX,1,NZ);
            if (FDORDER_TIME==4) {
                free_f3tensor(rxx_4,1,NY,1,NX,1,NZ);
                free_f3tensor(ryy_4,1,NY,1,NX,1,NZ);
                free_f3tensor(rzz_4,1,NY,1,NX,1,NZ);
                free_f3tensor(rxy_4,1,NY,1,NX,1,NZ);
                free_f3tensor(ryz_4,1,NY,1,NX,1,NZ);
                free_f3tensor(rxz_4,1,NY,1,NX,1,NZ);
            }
        }
		free_f3tensor(taus,0,NY+1,0,NX+1,0,NZ+1);
		free_f3tensor(taup,0,NY+1,0,NX+1,0,NZ+1);
		free_vector(eta,1,L);
		free_f3tensor(tausipjp,1,NY,1,NX,1,NZ);
		free_f3tensor(tausjpkp,1,NY,1,NX,1,NZ);
		free_f3tensor(tausipkp,1,NY,1,NX,1,NZ);

	}

	free_f3tensor(rho,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(pi,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(u,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(absorb_coeff,1,NY,1,NX,1,NZ);

	/* averaged material parameters */
	free_f3tensor(uipjp,1,NY,1,NX,1,NZ);
	free_f3tensor(ujpkp,1,NY,1,NX,1,NZ);
	free_f3tensor(uipkp,1,NY,1,NX,1,NZ);


	free_f3tensor(rjp,1,NY,1,NX,1,NZ);
	free_f3tensor(rkp,1,NY,1,NX,1,NZ);
	free_f3tensor(rip,1,NY,1,NX,1,NZ);

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

	/* free memory for global receiver source positions */
	if (SEISMO>0) {
		free_imatrix(recpos,1,3,1,ntr_glob);
	}

	/* free memory for global source positions */
	free_matrix(srcpos,1,6,1,nsrc);


	if ((ntr>0) && (SEISMO>0)){

		free_imatrix(recpos_loc,1,3,1,ntr);
		free_matrix(seismo_fulldata,1,ntr_glob,1,ns);

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
	free_ivector(stype,1,nsrc);
	free_ivector(stype_loc,1,nsrc);

	/* de-allocate buffer for messages */
	MPI_Buffer_detach(buff_addr,&buffsize);


	MPI_Barrier(MPI_COMM_WORLD);

	/* merge snapshot files created by the PEs into one file */
	/* if ((SNAP) && (MYID==0)) snapmerge(nsnap);*/

	free_ivector(xb,0,1);
	free_ivector(yb,0,1);
	free_ivector(zb,0,1);

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
		fprintf(FP," **************** SOFI3D has finished *****************\n");
		fprintf(FP," ******************************************************\n\n");
	}


	fclose(FP);

	MPI_Finalize();

	return 0;

}
