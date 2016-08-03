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
/*-------------------------------------------------------------
 * - program SOFI3D, Check SSG O(2,4) for stability and grid dispersion.
 *   If the stability criterion is not fulfilled the program will
 *   terminate.                   
 * - Check data output directories and files for accessibility.
 +   If any is not accessible from any PE the program will terminate.
 *  ----------------------------------------------------------*/

#include <libgen.h>
#include <unistd.h>
#include "fd.h"

void checkfd(FILE *fp, float *** prho, float *** ppi, float *** pu,
		float *** ptaus, float *** ptaup, float *peta, float **srcpos, int nsrc, int **recpos, int ntr){

	/* external variables */
	extern float DX, DY, DZ, DT, TS, TIME, TSNAP2;
	//extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2;
	extern int NX, NY, NZ, L, MYID, IDX, IDY, IDZ, FW, POS[4], NT, NDT, NDTSHIFT;
	extern int FDCOEFF, ABS_TYPE;
	//extern int READREC;
	extern int NPROCX, NPROCY,NPROCZ, FW, SRCREC, FREE_SURF;
	extern int SNAP, SEISMO, CHECKPTREAD, CHECKPTWRITE, SEIS_FORMAT[6], SNAP_FORMAT;
	extern int FDORDER, FDORDER_TIME;
	/*extern int RUN_MULTIPLE_SHOTS; no determination is done for the output check whether the simulation runs with one or multiple shot
		-> directorys specified in input file should work in both cases */
	extern char SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE], SNAP_FILE[STRING_SIZE];
	extern char SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];

	extern int BOUNDARY;

	/* local variables */
	float  c=0.0, cmax_p=0.0, cmin_p=1e9, cmax_s=0.0, cmin_s=1.0e9, fmax, cwater=1.0e-1;
	float  q=0.0, qmax_p=0.0, qmin_p=1e9, qmax_s=0.0, qmin_s=1.0e9;
	float snapoutx=0.0, snapouty=0.0, snapoutz=0.0, dhmax, dhmin;
    float  cmax=0.0, cmin=1e9, sum, ts, cmax_r, cmin_r, qmax_r, g=0.0;
    // float gamma=0.0; // isnt needed anymore
	float srec_minx=DX*NX*NPROCX+1, srec_miny=DY*NY*NPROCY+1, srec_minz=DZ*NZ*NPROCZ+1;
	float srec_maxx=-1.0, srec_maxy=-1.0, srec_maxz=-1.0;
	const float w=2.0*PI/TS; /*center frequency of source*/
    
    /* Variables to check stability */
    float CFL;
    float (*CFL_STAB)[3];
    float CFL_STAB_Taylor[6][3]= {
        {0.577, 0.494, 0.384},
        {0.494, 0.424, 0.329},
        {0.464, 0.398, 0.309},
        {0.448, 0.484, 0.299},
        {0.438, 0.375, 0.292},
        {0.431, 0.369, 0.287}
    };
    float CFL_STAB_Holberg[6][3]= {
        {0.577, 0.494, 0.384},
        {0.487, 0.418, 0.325},
        {0.45, 0.386, 0.3},
        {0.429, 0.486, 0.286},
        {0.416, 0.357, 0.277},
        {0.407, 0.349, 0.287}
    };
    
    
	int nfw=FW;
	int i, j, k, l, ny1=1, nx, ny, nz;

	char xfile[STRING_SIZE], errormessage[STRING_SIZE], xmod[4], file_ext[8];
	FILE *fpcheck;


	nx=NX; ny=NY; nz=NZ;

	/* low Q frame not yet applied as a absorbing boundary */
	/* if (!FREE_SURF) ny1=1+nfw;*/
	nfw=0;


	/*printf(" checkfd: TS= %f \n",TS);*/
	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ************ CHECKS OF INPUT FILE PARAMETERS  ************");
	fprintf(fp,"\n **********************************************************\n\n");
	fprintf(fp,"\n **Message from checkfd (printed by PE %d):\n",MYID);
	fprintf(fp,"\n\n ------------------ CHECK OUTPUT FILES --------------------------\n");	

	/* The original checks might delete files accidentally that would not be overwritten anyway. 
	   and did not test accessibility from all CPUs which may be vary, especially in distributed clusters */

	/*Checking SNAP Output */
	/*-------------------- */

	if ((SNAP>0) && (MYID==0)) {

		switch(SNAP_FORMAT){
		case 1:
			sprintf(file_ext,".su");
			strcpy(xmod,"ab");
			break;
		case 2:
			sprintf(file_ext,".asc");
			strcpy(xmod,"a");
			break;
		case 3:
			sprintf(file_ext,".bin");
			strcpy(xmod,"ab");
			break;
		default: err(" Sorry. Snapshot format (SNAP_FORMAT) unknown. \n");
		break;
		}		


		fprintf(fp," Check accessibility for snapshot files ... \n");

		switch(SNAP){
		case 1 :
			sprintf(xfile,"%s%s.vx.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0  cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			sprintf(xfile,"%s%s.vy.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			sprintf(xfile,"%s%s.vz.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			break;
		case 2 :
			sprintf(xfile,"%s%s.p.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			break;
		case 4 :
			sprintf(xfile,"%s%s.vx.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			sprintf(xfile,"%s%s.vy.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			sprintf(xfile,"%s%s.vz.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			sprintf(xfile,"%s%s.p.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
		case 3 :
			sprintf(xfile,"%s%s.div.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			sprintf(xfile,"%s%s.curl.%i%i%i",SNAP_FILE,file_ext,POS[1],POS[2],POS[3]);
			fprintf(fp,"    Check accessibility for snapshot file %s... \n",xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);
			else fclose(fpcheck);
			break;
		}
	}


	/*Checking SEISMOGRAM Output Particle velocities */
	/*-------------------------------------- */
	if ((SEISMO>0) && (MYID==0)) {
		switch (SEIS_FORMAT[0]){
		case 0: sprintf(file_ext,"sgy"); break;
		case 1: sprintf(file_ext,"su");  break;
		case 2: sprintf(file_ext,"txt"); break;
		case 3: sprintf(file_ext,"bin"); break;
		case 4: sprintf(file_ext,"sgy"); break;
		case 5: sprintf(file_ext,"sgy"); break;
		}

		fprintf(fp," Check accessibility for seismogram files of each PE ... \n");
		fprintf(fp," However, the list below refers only to PE0  ... \n");

		if (SEIS_FORMAT[0]==2) strcpy(xmod,"a");
		else strcpy(xmod,"w+b");

		/* MYID=0 is checking if seismogram files can be written to file
		 * since seismogram pieces from all PEs are exchanged and merged before the output
		 * only MYID0 has to check this	 */

		switch (SEISMO){
		case 1: /* particle velocities only */
			sprintf(xfile,"%s_vx.%s.%d",SEIS_FILE,file_ext,MYID);
			/*in case of number of PE's=500, there will be 500 messages, too many to be displayed! */
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			//if (access(xfile,W_OK|X_OK)==-1) err(errormessage);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			sprintf(xfile,"%s_vy.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			sprintf(xfile,"%s_vz.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			break;

		case 2 : /* pressure only */
			sprintf(xfile,"%s_p.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			break;

		case 4 : /* everything */
			sprintf(xfile,"%s_vx.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			sprintf(xfile,"%s_vy.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			sprintf(xfile,"%s_vz.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			sprintf(xfile,"%s_p.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

		case 3 : /* curl and div only */
			sprintf(xfile,"%s_div.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			sprintf(xfile,"%s_curl.%s.%d",SEIS_FILE,file_ext,MYID);
			fprintf(fp,"    Check accessibility for seismogram file %s... \n",xfile);
			sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);
			if ((fpcheck=fopen(xfile,xmod))==NULL) err(errormessage);
			else fclose(fpcheck);
			remove(xfile);

			break;
		}
	}

	/*Checking CHECKPOINT Output */
	/*-------------------------- */
	if (CHECKPTREAD>0) {
		strcpy(xmod,"rb");
		sprintf(xfile,"%s.%d",CHECKPTFILE,MYID);
		fprintf(fp," Check readability for checkpoint files %s... \n",xfile);
		if (((fpcheck=fopen(xfile,xmod))==NULL) && (MYID==0)) err(" PE 0 cannot read checkpoints!");
		else fclose(fpcheck);
	}
	if ((CHECKPTWRITE>0)){ 
		strcpy(xmod,"ab");
		sprintf(xfile,"%s.%d",CHECKPTFILE,MYID);
		fprintf(fp," Check writability for checkpoint files %s... \n",xfile);
		if (((fpcheck=fopen(xfile,xmod))==NULL) && (MYID==0)) err(" PE 0 cannot write checkpoints!");
		else fclose(fpcheck); /* Is there any reason to remove it? */	
		if(MYID==0){
			if (FDORDER_TIME > 2){warning(" Checkpoint writing for FDORDER_TIME > 2 is not implemented,\n wavefields calculated from the checkpoint will slightly differ.\n");}
		}	
		
	}

	fprintf(fp," Accessibility of output files from PE %d has been checked successfully.\n", MYID);

	fprintf(fp,"\n\n --------- DETERMININATION OF MIN AND MAX VELOCITIES -----------\n");	
	fprintf(fp," Minimum and maximum P-wave and S-wave velocities within subvolume: \n ");
	/* find maximum model phase velocity of shear waves at infinite
	      frequency within the whole model */
	if (L){ /*viscoelastic simulation */
		for (k=1+nfw;k<=(nz-nfw);k++){
			for (i=1+nfw;i<=(nx-nfw);i++){
				for (j=ny1;j<=(ny-nfw);j++){
					if (prho[j][i][k]==0.0)
						printf("rho=0 detected from MYID=%d at i=%d,j=%d,k=%d\n",MYID,i,j,k);
					else
						c=sqrt(pu[j][i][k]*(1.0+L*ptaus[j][i][k])/prho[j][i][k]);
					if (cmax_s<c) cmax_s=c;
					q=2/ptaus[j][i][k];
					if (qmax_s<q) qmax_s=q;
					/* find minimum model phase velocity of shear waves at center
					       frequency of the source */
					sum=0.0;
					for (l=1;l<=L;l++){
						ts=DT/peta[l];
						sum=sum+((w*w*ptaus[j][i][k]*ts*ts)/(1.0+w*w*ts*ts));
					}
					if (prho[j][i][k]==0.0)
						printf("rho=0 detected from MYID=%d at i=%d,j=%d,k=%d\n",MYID,i,j,k);
					else
						c=sqrt(pu[j][i][k]/prho[j][i][k]);
					if ((c>cwater)&&(cmin_s>c)) cmin_s=c;
					if (qmin_s>q) qmin_s=q;
				}
			}
		}
	}
	else { /*elastic simulation */
		for (k=1+nfw;k<=(nz-nfw);k++){
			for (i=1+nfw;i<=(nx-nfw);i++){
				for (j=ny1;j<=(ny-nfw);j++){
					if (prho[j][i][k]==0.0)
						printf("rho=0 detected from MYID=%d at i=%d,j=%d,k=%d\n",MYID,i,j,k);
					else
						c=sqrt(pu[j][i][k]/prho[j][i][k]);
					if ((c>cwater) && (cmax_s<c)) cmax_s=c;
					if ((c>cwater) && (cmin_s>c)) cmin_s=c;
				}
			}
		}
	}


	/* find maximum model phase velocity of P-waves at infinite
		 frequency within the whole model */
	if (L){ /*viscoelastic simulation */
		for (k=1+nfw;k<=(nz-nfw);k++){
			for (i=1+nfw;i<=(nx-nfw);i++){
				for (j=ny1;j<=(ny-nfw);j++){
					if (prho[j][i][k]==0.0)
						printf("prho=0 detected from MYID=%d at i=%d,j=%d,k=%d\n",MYID,i,j,k);
					else
						c=sqrt(ppi[j][i][k]*(1.0+L*ptaup[j][i][k])/prho[j][i][k]);
					q=2/ptaus[j][i][k];
					if (qmax_p<q) qmax_p=q;
					if ((c>cwater)&&(cmax_p<c)) cmax_p=c;
					/* find minimum model phase velocity of P-waves at center
					       frequency of the source */
					sum=0.0;
					for (l=1;l<=L;l++){
						ts=DT/peta[l];
						sum=sum+((w*w*ptaup[j][i][k]*ts*ts)/(1.0+w*w*ts*ts));
					}
					if (prho[j][i][k]==0.0)
						printf("prho=0 detected from MYID=%d at i=%d,j=%d,k=%d\n",MYID,i,j,k);
					else
						c=sqrt(ppi[j][i][k]/prho[j][i][k]);
					if ((c>cwater)&&(cmin_p>c)) cmin_p=c;
					if (qmin_p<q) qmin_p=q;
				}
			}
		}
	}
	else { /*elastic simulation */
		for (k=1+nfw;k<=(nz-nfw);k++){
			for (i=1+nfw;i<=(nx-nfw);i++){
				for (j=ny1;j<=(ny-nfw);j++){
					if (prho[j][i][k]==0.0)
						printf("rho=0 detected from MYID=%d at i=%d,j=%d,k=%d\n",MYID,i,j,k);
					else
						c=sqrt(ppi[j][i][k]/prho[j][i][k]);
					if ((c>cwater) && (cmax_p<c)) cmax_p=c;
					if ((c>cwater) && (cmin_p>c)) cmin_p=c;					}
			}
		}
	}


	fprintf(fp," MYID\t Vp_min(f=fc) \t Vp_max(f=inf) \t Vs_min(f=fc) \t Vs_max(f=inf) \n");
	fprintf(fp," %d \t %8.2f \t %8.2f \t %8.2f \t %8.2f \n\n\n", MYID, cmin_p, cmax_p, cmin_s, cmax_s);

	fprintf(fp," Note : if any P- or S-wave velocity is set below 1.0 m/s to simulate water or air,\n");
	fprintf(fp," this minimum velocity will be ignored for determining stable DH and DT.\n\n");

	if (cmax_s>cmax_p) cmax=cmax_s; 
	else cmax=cmax_p;
	if (cmin_s<cmin_p) cmin=cmin_s; 
	else cmin=cmin_p;

	/* find global maximum for Vp and global minimum for Vs*/
	MPI_Allreduce(&cmax,&cmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&cmin,&cmin_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	cmax=cmax_r;
	cmin=cmin_r;
	/* find global maximum for qp and qs and global minimum for qp and qs*/
	MPI_Allreduce(&qmax_p,&cmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	qmax_p=cmax_r;
	MPI_Allreduce(&qmax_s,&qmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	qmax_s=cmax_r;
	MPI_Allreduce(&qmin_p,&qmax_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	qmin_p=cmax_r;
	MPI_Allreduce(&qmin_s,&qmax_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	qmin_s=cmax_r;

	/* if (MYID==0){		checkfd is performed by MID=0 only */

	fprintf(fp," Global values for entire model: \n");
	fprintf(fp," Vp_max= %8.2f m/s \t Vs_min=%8.2f m/s \n", cmax,cmin);
	if (L) {
		fprintf(fp," \n Please note that for viscoelastic modeling V_max and V_min can slightly differ \n");
		fprintf(fp," from your actual velocity extrema due to velocity dispersion by viscoelastic relaxation!\n");
		fprintf(fp,"\n Minimum and maximum Qp and Qs damping within subvolume: \n ");
		fprintf(fp," MYID\t Qp_min \t Qp_max \t Qs_min \t Qs_max \n");
		fprintf(fp," %d \t %3.5f \t %3.5f \t %3.5f \t %3.5f \n\n", MYID, qmin_p, qmax_p, qmin_s, qmax_s);
	}


	/* estimate number of grid points per minimum wavelength to avoid numerical dispersion */
	/* Taylor */
	if (FDCOEFF==1) {
		switch (FDORDER){
		case 2: g=12.0;
		break;
		case 4: g=8.0;
		break;
		case 6: g=7.0;
		break;
		case 8: g=6.0;
		break;
		case 10: g=5.0;
		break;
		case 12: g=4.0;
		break;
		default: err(" invalid FDORDER in checkfd_ssg() !");
		break;
		}
	}

	/* Holberg */
	if (FDCOEFF==2) {
		switch (FDORDER){
		case 2: g=12.0;
		break;
		case 4: g=8.32;
		break;
		case 6: g=4.77;
		break;
		case 8: g=3.69;
		break;
		case 10: g=3.19;	
		break;
		case 12: g=2.91;
		break;
		default: err(" invalid FDORDER in checkfd_ssg() !");
		break;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/******************************************************************************************/

	/* calculate maximum grid point distance */
	dhmax=DX;
	if(DY>dhmax){dhmax=DY;}
	if(DZ>dhmax){dhmax=DZ;}	

	fprintf(fp,"\n ------------------ CHECK FOR GRID DISPERSION --------------------\n");
	fprintf(fp," To satisfactorily limit grid dispersion the number of gridpoints \n");
	fprintf(fp,"    per minimum wavelength (of S-waves) should be 6 (better more).\n");
	fprintf(fp," Here the minimum wavelength is assumed to be minimum model phase velocity \n");
	fprintf(fp,"    (of S-waves) at maximum frequency of the source\n");
	fprintf(fp,"    devided by maximum frequency of the source.\n");
	fmax=2.0/TS;
	fprintf(fp," Maximum frequency of the source is approximately %8.2f Hz\n",fmax);
	fprintf(fp,"    (= two times FC specified in the sources file)\n");
	fprintf(fp," The minimum wavelength (of P- or S-waves) in the following simulation will\n");
	fprintf(fp,"    be %4.4f meter.\n", cmin/fmax);
	fprintf(fp," Thus, at the chosen FD_ORDER of %i the recommended value for DH is %4.4f meter \n    (%2.2f gridpoints per minimum wavelength).\n", FDORDER,cmin/fmax/g,g);
	fprintf(fp," You have specified DH= %4.4f meter.\n\n", dhmax);
	if ((dhmax>(cmin/fmax/g))&&(MYID==0))
		warning(" Grid dispersion will influence wave propagation, choose smaller grid spacing (DH).");

	/* calculate minimum grid point distance */
	dhmin=DX;
	if(DY<dhmin){dhmin=DY;}
	if(DZ<dhmin){dhmin=DZ;}
    
    /* calculate the courant friedrich lewy number */
    CFL=cmax*DT/dhmax;
    /* Stability limits for the simulations in terms of Courant-Friedrich-Lewy: Switch Coefficients */
    switch (FDCOEFF) {
        case 1: /* Taylor */
            CFL_STAB=CFL_STAB_Taylor;
            break;
            
        case 2: /* Holberg */
            CFL_STAB=CFL_STAB_Holberg;
            break;
            
        default:
            CFL_STAB=CFL_STAB_Taylor;
    }
    
	fprintf(fp," \n\n ----------------------- CHECK FOR STABILITY ---------------------\n");
	fprintf(fp," The following simulation is stable provided that\n\n");
	fprintf(fp," \t p=cmax*DT/DH <= %f \n\n",CFL_STAB[FDORDER/2-1][FDORDER_TIME-2]);
	fprintf(fp," where cmax is the maximum phase velocity at infinite frequency,\n");
	fprintf(fp," In the current simulation cmax is %8.2f m/s .\n",cmax);
	fprintf(fp," DT is the timestep and DH is the grid size.\n\n");
    fprintf(fp," In this simulation the maximum Courant-Friedrich-Lewey (CFL) number p is: %f \n",CFL_STAB[FDORDER/2-1][FDORDER_TIME-2]);
    fprintf(fp," The CFL-number in this simulation will be p= %f \n",CFL);
	fprintf(fp," In this simulation the stability limit for timestep DT is %e seconds .\n",CFL_STAB[FDORDER/2-1][FDORDER_TIME-2]*DX/cmax);
	fprintf(fp," You have specified DT= %e s.\n", DT);

    if (CFL>CFL_STAB[FDORDER/2-1][FDORDER_TIME-2])
        err(" The simulation will get unstable, choose smaller DT. ");
    else fprintf(fp," The simulation will be stable.\n");
    

	fprintf(fp," \n\n --------------------- CHECK FOR INPUT ERRORS ---------------------\n");
    
//    /* Higher temporal orders only work in the elastic case */
//    if(L==1 && FDORDER_TIME>2){
//        err("Higher temporal orders (FDORDER_TIME>2) only work in the elastic case.");
//    }
    
	if (SNAP){
		fprintf(fp," Checking the snapshot parameters. \n");
		if (TSNAP2>TIME) {
			sprintf(errormessage,"\nTSNAP2 = %e (last snapshot) > Time of wave propagation %e. TSNAP2 was changed to be equal to TIME.\n",TSNAP2, TIME); 
			TSNAP2=TIME; 
		if  (MYID==0)
			warning(errormessage); /* if TSNAP2>simulation TIME, snapmerge will generate "additional" snapshots out of nowhere, thus, snapshot files size blow up */
		}

		snapoutx=NX/(float)IDX;
		snapouty=NY/(float)IDY;
		snapoutz=NZ/(float)IDZ;
		fprintf(fp,"    Output of snapshot gridpoints per node (NX/NPROCX/IDX) %8.2f .\n", snapoutx);
		fprintf(fp,"    Output of snapshot gridpoints per node (NY/NPROCY/IDY) %8.2f .\n", snapouty);
		fprintf(fp,"    Output of snapshot gridpoints per node (NZ/NPROCZ/IDZ) %8.2f .\n", snapoutz);

		if (snapoutx-(int)snapoutx>0)
			err("\n\n Ratio NX-NPROCX-IDX must be whole-numbered \n\n");
		if (snapouty-(int)snapouty>0)
			err("\n\n Ratio NY-NPROCY-IDY must be whole-numbered \n\n");
		if (snapoutz-(int)snapoutz>0)
			err("\n\n Ratio NZ-NPROCZ-IDZ must be whole-numbered \n\n");
	}

	if ((SEISMO)&& (MYID==0)){
		fprintf(fp,"\n Checking the number of seismogram samples. \n");
		fprintf(fp,"    Number of timesteps %d.\n", NT);
		fprintf(fp,"    Seismogram sampling rate in timesteps %d.\n", NDT);
		fprintf(fp,"    Number of seismogram output samples %d.\n", NT/NDT-NDTSHIFT);
	} 

	if ((SEISMO>0) && (MYID==0)) {
		srec_minx=DX*NX*NPROCX+1, srec_miny=DY*NY*NPROCY+1, srec_minz=DZ*NZ*NPROCZ+1;
		srec_maxx=-1.0, srec_maxy=-1.0, srec_maxz=-1.0;
		fprintf(fp,"\n Checking for receiver position(s) specified in input file.\n");
		fprintf(fp,"    Global grid size in m : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) : %5.2f-%5.2f (z in m).\n",DX,NX*DX*NPROCX,DY,NY*DY*NPROCY,DZ,NZ*DZ*NPROCZ);
		if (FREE_SURF==0) fprintf(fp,"    Global grid size in m (-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) : %5.2f-%5.2f (z in m).\n",FW*DX,NX*DX*NPROCX-FW*DX,FW*DY,NY*DY*NPROCY-FW*DY,FW*DZ,NZ*DZ*NPROCZ-FW*DZ);
		if (FREE_SURF==1) fprintf(fp,"    Global grid size in m (-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) : %5.2f-%5.2f (z in m).\n",FW*DX,NX*DX*NPROCX-FW*DX,DY,NY*DY*NPROCY-FW*DY,FW*DZ,NZ*DZ*NPROCZ-FW*DZ);

		/* find maximum and minimum source positions coordinate ---- from input file*/
		/*note that "y" is used for the vertical coordinate */

		/* in case of mod(XREC1,DX)!=0 or any other Receiver coordinate - spatial discretization pair
		 * specified in the input file, it might be useful to actually determine the min and max values
		 * by looking in the rec_pos variable */

		/*if (READREC==0) {
			if (XREC1>XREC2) {
				srec_maxx=XREC1;
				srec_minx=XREC2;
			}
			else {
				srec_maxx=XREC2;
				srec_minx=XREC1;
			}
			if (YREC1>YREC2) {
				srec_maxy=YREC1;
				srec_miny=YREC2;
			}
			else {
				srec_maxy=YREC2;
				srec_miny=YREC1;
			}
			if (ZREC1>ZREC2) {
				srec_maxz=ZREC1;
				srec_minz=ZREC2;
			}
			else {
				srec_maxz=ZREC2;
				srec_minz=ZREC1;
			}
		}*/
		//if (READREC==1) {

		/* find maximum and minimum source positions coordinate ---- from matrix rec_pos*/
		for (k=1;k<=ntr;k++){
			/* find maximum source positions coordinate*/
			if ((recpos[1][k]*DX)>srec_maxx) srec_maxx=(float)recpos[1][k]*DX;
			if ((recpos[2][k]*DY)>srec_maxy) srec_maxy=(float)recpos[2][k]*DY;
			if ((recpos[3][k]*DZ)>srec_maxz) srec_maxz=(float)recpos[3][k]*DZ;
			/* find minimum source positions coordinate*/
			if ((recpos[1][k]*DX)<srec_minx) srec_minx=(float)recpos[1][k]*DX;
			if ((recpos[2][k]*DY)<srec_miny) srec_miny=(float)recpos[2][k]*DY;
			if ((recpos[3][k]*DZ)<srec_minz) srec_minz=(float)recpos[3][k]*DZ;
		}
		fprintf(fp,"    Number of receiver positions in receiver file %s : %i \n", REC_FILE, ntr);


		fprintf(fp,"    Minimum receiver position coordinates : %5.2f (x) : %5.2f (y) : %5.2f (z).\n",srec_minx,srec_miny, srec_minz);
		fprintf(fp,"    Maximum receiver position coordinates : %5.2f (x) : %5.2f (y) : %5.2f (z).\n",srec_maxx,srec_maxy, srec_maxz);
		/* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
		if ((((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_maxz<0.0))) || (((srec_minx<0.0) || (srec_miny<0.0)) || ((srec_minz<0.0)))) {
			err("\n\n Coordinate of at least one receiver location is outside the global grid. \n\n");
		}
		if (((srec_maxx>NX*DX*NPROCX) || (srec_maxz>NZ*DZ*NPROCZ)) || ((srec_maxy>NY*DY*NPROCY))) {
			err("\n\n Coordinate of at least one receiver location is outside the global grid. \n\n");
		}
		/* checking if receiver coordinate of first receiver in line specified in input-file is inside the Absorbing Boundary  */
		if ((((srec_maxx<(FW*DX)) || (srec_maxz<(FW*DZ))) || ((srec_maxy<(FW*DY))&& !(FREE_SURF==1))) || (((srec_minx<(FW*DX)) || (srec_minz<(FW*DZ)))) || ((srec_miny<(FW*DY)) && !(FREE_SURF==1))) {
			/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
			warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 1). \n\n");
		}
		if (((srec_maxx>(NX*DX*NPROCX-FW*DX)) || (srec_maxy>(NY*DY*NPROCY-FW*DY)) || ((srec_maxz>(NZ*DZ*NPROCZ-FW*DZ))))) {
			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
			warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 2). \n\n");
		}

		fprintf(fp," ... complete, receiver position specified in input file are located within the global grid.\n\n");

	}

	if ((SRCREC==1)&& (MYID==0)){
		srec_minx=DX*NX*NPROCX+1, srec_miny=DY*NY*NPROCY+1, srec_minz=DZ*NZ*NPROCZ+1;
		srec_maxx=-1.0, srec_maxy=-1.0, srec_maxz=-1.0;
		fprintf(fp,"\n Checking for source position(s) specified in source file. \n");
		fprintf(fp,"    Global grid size in m: %5.2f (x) : %5.2f (y) : %5.2f (z) :.\n",NX*DX*NPROCX,NY*DY*NPROCY,NZ*DZ*NPROCZ);
		if (FREE_SURF==0) fprintf(fp,"    Global grid size in m (-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) : %5.2f-%5.2f (z in	m).\n",FW*DX,NX*DX*NPROCX-FW*DX,FW*DZ,NY*DY*NPROCY-FW*DY,FW*DZ,NZ*DZ*NPROCZ-FW*DZ);
		if (FREE_SURF==1) fprintf(fp,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) : %5.2f-%5.2f (z in m).\n",FW*DX,NX*DX*NPROCX-FW*DX,FW*DY,NZ*DY*NPROCY-FW*DY,DZ,NZ*DZ*NPROCZ-FW*DZ);

		/*only for testing
		fprintf(fp," initial min : %5.2f (x) %5.2f (y) %5.2f (z). \n",srec_minx,srec_miny,srec_minz);
		fprintf(fp," initial max : %5.2f (x) %5.2f (y) %5.2f (z). \n",srec_maxx,srec_maxy,srec_maxz);

		fprintf(fp," initial model size  : %i (Nx) %i (Ny) %i (Nz). \n",NX,NY,NZ);
		fprintf(fp," initial model size  : %5.2f (dx) %5.2f (dy) %5.2f (dz). \n",DX,DY,DZ);*/
		for (k=1;k<=nsrc;k++){
			/* find maximum source positions coordinate*/
			if (srcpos[1][k]>srec_maxx) srec_maxx=srcpos[1][k];
			if (srcpos[2][k]>srec_maxy) srec_maxy=srcpos[2][k];
			if (srcpos[3][k]>srec_maxz) srec_maxz=srcpos[3][k];
			/* find minimum source positions coordinate*/
			if (srcpos[1][k]<srec_minx) srec_minx=srcpos[1][k];
			if (srcpos[2][k]<srec_miny) srec_miny=srcpos[2][k];
			if (srcpos[3][k]<srec_minz) srec_minz=srcpos[3][k];
		}
		/*only for testing
		fprintf(fp," new min : %5.2f (x) %5.2f (y) %5.2f (z). \n",srec_minx,srec_miny,srec_minz);
		fprintf(fp," new max : %5.2f (x) %5.2f (y) %5.2f (z). \n",srec_maxx,srec_maxy,srec_maxz);*/

		fprintf(fp,"    Number of source positions in source file %s. : %i.\n", SOURCE_FILE, nsrc);
		fprintf(fp,"    Minimum source position coordinates : %5.2f (x) : %5.2f (y) : %5.2f (z).\n",srec_minx,srec_miny, srec_minz);
		fprintf(fp,"    Maximum source position coordinates : %5.2f (x) : %5.2f (y) : %5.2f (z).\n",srec_maxx,srec_maxy, srec_maxz);
		/* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
		if ((((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_maxz<0.0))) || (((srec_minx<0.0) || (srec_miny<0.0)) || ((srec_minz<0.0)))) {
			err("\n\n Coordinate of at least one source location is outside the global grid. \n\n");
		}
		if (((srec_maxx>NX*DX*NPROCX) || (srec_maxy>NY*DY*NPROCY)) || ((srec_maxz>NZ*DZ*NPROCZ))) {
			err("\n\n Coordinate of at least one source location is outside the global grid. \n\n");
		}
		/* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
		if ((((srec_maxx<(FW*DX)) || (srec_maxz<(FW*DZ))) || ((srec_maxy<(FW*DY))&& !(FREE_SURF==1))) || (((srec_minx<(FW*DX)) || (srec_minz<(FW*DZ)))) || ((srec_miny<(FW*DY)) && !(FREE_SURF==1))) {
			/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
			warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 1). \n\n");
		}
		if (((srec_maxx>(NX*DX*NPROCX-FW*DX)) || (srec_maxy>(NY*DY*NPROCY-FW*DY)) || ((srec_maxz>(NZ*DZ*NPROCZ-FW*DZ))))) {
			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
			warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 2). \n\n");
		}


		fprintf(fp," ...complete, all source position(s) specified in source file are located within the global grid.\n");


	}

	/* 		if ((NT/NDT-NDTSHIFT)>32000)
	err("\n\n Seismogram samples must be less than 32.0000 \n\n"); */ /*moved here from inside the previous statement for reason below */
	/* checked in writepar.c now. SU and SEG-Y allow 32767 samples, furthermore
	the exist programs allowing 65535 samples and pseudo-SEG-Y formats allowing
	almost arbitrarily long traces. And for binary and textual output the limit is
	arbitrary. */


	fprintf(fp,"\n\n ----------------------- ABSORBING BOUNDARY ------------------------\n");

	if (ABS_TYPE==2){

		fprintf(fp," Width (FW) of absorbing frame should be at least 30 gridpoints.\n");
		fprintf(fp," You have specified a width of %i gridpoints.\n",FW);
		if ((FW<30)&&(MYID==0))
			warning(" Be aware of strong artificial reflections from grid boundaries ! \n");
	}

	if (ABS_TYPE==1){
		fprintf(fp," Width (FW) of CPML-frame should be at least 20 gridpoints.\n");
		fprintf(fp," You have specified a width of %i gridpoints.\n",FW);
		if ((FW<20)&&(MYID==0))
			warning(" Be aware of strong artificial reflections from grid boundaries ! \n");
		if(FDORDER>4)
			warning(" CPML Update of FDORDER > 4 not implemented yet, coefficients of FDORDER==4 are used instead!");
	}
	if (((((NX)<FW) || ((NY)<FW)) || ((NZ)<FW)) && (MYID==0))	{
		fprintf(fp," \n Width of boundary area (FW=%i gridpoints) is larger than at least one subdomain dimension: \n",FW);
		fprintf(fp," \t NX/NPROCX = %i, NY/NPROCY = %i, NZ/NPROCZ = %i gridpoints.\n",NX,NY,NZ);
		err(" Choose smaller width of absorbing frame (ABS) or increase subdomain dimensions");
	}

	// Check if there is a conflict between BOUNDARY and ABS_TYPE
	if(BOUNDARY)
		if(ABS_TYPE==1 || ABS_TYPE==2)
		{
			warning(" You have activated a periodic boundary and set an absorbing boundary at the same time! \n");
		}

}

