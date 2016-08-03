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
 *   Write 3D snapshot for current timestep  to disk                                   
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void snap_acoustic(FILE *fp, int nt, int nsnap, int format, int type, 
float ***vx, float ***vy, float ***vz, float ***sxx, float ***pi,
int idx, int idy, int idz, int nx1, int ny1, int nz1, int nx2, 
int ny2, int nz2){

	/* 
	different data formats of output available:
	format=1  :  SU (IEEE)
	format=2  :  ASCII
	format=3  :  BINARY (IEEE)
	
	different types:
	type=1 : values in vx, vy, and vz
	type=2 : -(sxx+syy+szz) (pressure field)
	type=3 : divergence of vx, vy and vz (energy of compressional waves)
	         and curl of vx, vy and vz (energy of shear waves)
	type=4 : both particle velocities (type=1) and energy (type=3)
	*/


	
	char xfile[STRING_SIZE], yfile[STRING_SIZE], zfile[STRING_SIZE];
	char rotfile[STRING_SIZE], ext[8], wm[2];
	char  divfile[STRING_SIZE], pfile[STRING_SIZE];
	FILE *fpx1, *fpy1, *fpz1, *fpp /*, *fpx2, *fpy2*/;
	int i,j,k;
	float /*a=0.0,*/ amp; /* dh24x, dh24y, dh24z, vyx, vxy, vxx, vyy, vzx, vyz, vxz, vzy, vzz not used here*/


	extern float /*DX, DY, DZ,*/ DT;
	extern char SNAP_FILE[STRING_SIZE];
	extern int MYID, POS[4], LOG; /* SNAP_PLANE not used here*/



	switch(format){
	case 1: 
		sprintf(ext,".su");
		break;
	case 2: 
		sprintf(ext,".asc");
		break;
	case 3: 
		sprintf(ext,".bin");
		break;
	}


	sprintf(xfile,"%s%s.x.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(yfile,"%s%s.z.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(zfile,"%s%s.y.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(divfile,"%s%s.div.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(rotfile,"%s%s.rot.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(pfile,"%s%s.p.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);

        if (LOG){
	fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);}


	if (nsnap==1) 
		sprintf(wm,"w");
	else 
		sprintf(wm,"a");

	switch(type){
	case 1 :
		fprintf(fp,"\t%s\n", xfile);
		fprintf(fp,"\t%s\n", yfile);
		fprintf(fp,"\t%s\n\n", zfile);
		fpx1=fopen(xfile,wm);
		fpy1=fopen(yfile,wm);
		fpz1=fopen(zfile,wm);
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
				
			
					writedsk(fpx1,vx[j][i][k],format);
					writedsk(fpy1,vy[j][i][k],format);
					writedsk(fpz1,vz[j][i][k],format);
					
					
				}
		fclose(fpx1);
		fclose(fpy1);
		fclose(fpz1);
		break;
	case 2 :
		fprintf(fp,"\t%s\n\n", pfile);
		fpp=fopen(pfile,wm);
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
					amp=-3.0*sxx[j][i][k];
					
				
					writedsk(fpp,amp,format);

				}
		fclose(fpp);
		break;
	case 4 :
		fprintf(fp,"\t%s\n", xfile);
		fprintf(fp,"\t%s\n", yfile);
		fprintf(fp,"\t%s\n\n", zfile);
		fprintf(fp,"\t%s\n\n", pfile);
		fpx1=fopen(xfile,wm);
		fpy1=fopen(yfile,wm);
		fpz1=fopen(zfile,wm);
		fpp=fopen(pfile,wm);
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
					amp=-3.0*sxx[j][i][k];

					writedsk(fpx1,vx[j][i][k],format);
					writedsk(fpy1,vy[j][i][k],format);
					writedsk(fpz1,vz[j][i][k],format);
					
					writedsk(fpp,amp,format);
				}
		fclose(fpx1);
		fclose(fpy1);
		fclose(fpz1);
		fclose(fpp);
		break;
		
	}
}


