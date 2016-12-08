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
/* ------------------------------------------------------------------------
 *   Create dissipative boundarie around the model grid
 *   The dissipative coefficients are stored in the 3-D array
 *   absorb_coeff. The interior of the model is weighted by the
 *   coefficient 1. In the absorbing frame the coefficients 
 *   are less than one. Coefficients are computed using 
 *   exponential damping (see Cerjan et al., 1985, 
 *   Geophysics, 50, 705-708)
 *
 ------------------------------------------------------------------------*/

#include "fd.h"

void absorb(float *** absorb_coeff){

	/* extern variables */

	extern float DAMPING;
	extern int FREE_SURF, NX, NY, NZ, BOUNDARY, FW;
	extern int NPROCX, NPROCY, NPROCZ, MYID, POS[4];
	
	/* local variables */
	int i, j, k, ifw, ii, jj, kk, xb, yb, zb, xe, ye, ze;
	float amp, a, *coeff;
	char modfile[STRING_SIZE];
	extern FILE *FP;
	
	if (MYID==0){
		fprintf(FP,"\n **Message from absorb (printed by PE %d):\n",MYID);
		fprintf(FP," Coefficients for absorbing frame are now calculated.\n");
		fprintf(FP," Width of dissipative frame (grid points)= %i\n",FW);
		fprintf(FP," Percentage of exponential damping = %5.2f\n",DAMPING);
	}

	amp=1.0-DAMPING/100.0;   /* amplitude at the edge of the numerical grid */
	ifw=FW;  /* frame width in gridpoints */
	coeff=vector(1,ifw);
	a=sqrt(-log(amp)/((ifw-1)*(ifw-1)));
	
	for (i=1;i<=ifw;i++)
		coeff[i]=exp(-(a*a*(ifw-i)*(ifw-i)));
	
	if (MYID==0){
		fprintf(FP," Table of coefficients \n # \t coeff \n");
		/*fprintf(FP," ifw=%d \t a=%f amp=%f \n", ifw,a,amp); */
		for (i=1;i<=ifw;i++)
			fprintf(FP," %d \t %5.3f \n", i, coeff[i]);
	}	
	

	/* initialize array of coefficients with one */
	for (k=1;k<=NZ;k++)
	for (j=1;j<=NY;j++)
	for (i=1;i<=NX;i++) absorb_coeff[j][i][k]=1.0;


	/* compute coefficients for left and right grid boundaries (x-direction) */
	if ((!BOUNDARY) && (POS[1]==0)){
		zb=1; yb=1; ye=NY; ze=NZ;
		for (i=1;i<=ifw;i++){
			if (POS[3]==0) zb=i;
			if (POS[3]==NPROCZ-1) ze=NZ-i+1;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i;
			if (POS[2]==NPROCY-1) ye=NY-i+1;
			for (k=zb;k<=ze;k++)
				for (j=yb;j<=ye;j++)
					absorb_coeff[j][i][k]=coeff[i];
		}
	}
			
	if ((!BOUNDARY) && (POS[1]==NPROCX-1)){
		zb=1; yb=1; ye=NY; ze=NZ;
		for (i=1;i<=ifw;i++){
			ii=NX-i+1;
			if (POS[3]==0) zb=i;
			if (POS[3]==NPROCZ-1) ze=NZ-i+1;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i;
			if (POS[2]==NPROCY-1) ye=NY-i+1;
			for (k=zb;k<=ze;k++)
				for (j=yb;j<=ye;j++)
					absorb_coeff[j][ii][k]=coeff[i];
		}
	}
	
	/* compute coefficients for front and backward grid boundaries (z-direction) */
			
	if ((!BOUNDARY) && (POS[3]==0)){
		xb=1; yb=1; ye=NY; xe=NX;
		for (k=1;k<=ifw;k++){
			if (POS[1]==0) xb=k;
			if (POS[1]==NPROCX-1) xe=NX-k+1;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=k;
			if (POS[2]==NPROCY-1) ye=NY-k+1;
			for (i=xb;i<=xe;i++)
				for (j=yb;j<=ye;j++)
					absorb_coeff[j][i][k]=coeff[k];
		}
	}

			
	if ((!BOUNDARY) && (POS[3]==NPROCZ-1)){
		xb=1; yb=1; ye=NY; xe=NX;
		for (k=1;k<=ifw;k++){
			kk=NZ-k+1;
			if (POS[1]==0) xb=k;
			if (POS[1]==NPROCX-1) xe=NX-k+1;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=k;
			if (POS[2]==NPROCY-1) ye=NY-k+1;
			for (i=xb;i<=xe;i++)
				for (j=yb;j<=ye;j++)
					absorb_coeff[j][i][kk]=coeff[k];
		}
	}

	/* compute coefficients for top and bottom grid boundaries (y-direction) */

	if ((POS[2]==0) && (!(FREE_SURF))){
		xb=1; zb=1; ze=NZ; xe=NX;
		for (j=1;j<=ifw;j++){
			if ((!BOUNDARY) && (POS[1]==0)) xb=j;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1;
			if ((!BOUNDARY) && (POS[3]==0)) zb=j;
			if ((!BOUNDARY) && (POS[3]==NPROCZ-1)) ze=NZ-j+1;
			for (i=xb;i<=xe;i++)
				for (k=zb;k<=ze;k++)
					absorb_coeff[j][i][k]=coeff[j];
		}
	}

	if (POS[2]==NPROCY-1){
		xb=1; zb=1; ze=NZ; xe=NX;
		for (j=1;j<=ifw;j++){
			jj=NY-j+1;
			if ((!BOUNDARY) && (POS[1]==0)) xb=j;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1;
			if ((!BOUNDARY) && (POS[3]==0)) zb=j;
		   if ((!BOUNDARY) && (POS[3]==NPROCZ-1)) ze=NZ-j+1;
			for (i=xb;i<=xe;i++)
				for (k=zb;k<=ze;k++)
					absorb_coeff[jj][i][k]=coeff[j];
		}
	}


	sprintf(modfile,"model/absorb.bin");

	/*writemod(modfile,absorb_coeff,3); 

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3); */


	free_vector(coeff,1,ifw);
}



