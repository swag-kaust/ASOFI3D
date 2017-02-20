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

void absorb_PML(float *** absorb_coeffx, float *** absorb_coeffy, float *** absorb_coeffz){

	/* extern variables */

	extern float DAMPING, DX, VPPML;
	extern int FREE_SURF, NX, NY, NZ, BOUNDARY, FW;
	extern int NPROCX, NPROCY, NPROCZ, MYID, POS[4];
	
	/* local variables */
	int i, j, k, ifw, ii, jj, kk, xb, yb, zb, xe, ye, ze;
	float a, b, *coeff, c, dump, alpha, ffw; /*variable "DT1, c1" removed, not in use */
	//float amp, a0, pis, f0; /* variable "R" removed, not in use*/
	char modfile[STRING_SIZE];
	extern FILE *FP;
	
	if (MYID==0){
		fprintf(FP,"\n **Message from absorb (printed by PE %d):\n",MYID);
		fprintf(FP," Coefficcients for absorbing frame are now calculated.\n");
		fprintf(FP," Width of dissipative frame (meter)= %i\n",FW);
		fprintf(FP," Percentage of exponential damping = %5.2f\n",DAMPING);
	}

	ifw=FW;  /* frame width in gridpoints */
	ffw=(float)(FW*DX); /* frame width in m*/
	coeff=vector(1,ifw);
	//amp=1.0-DAMPING/100.0;   /* amplitude at the edge of the numerical grid */
	//pis=4.0*atan(1.0);
	//f0=467.0;
	//a0=1.79;

	alpha=1e-4;
	a=0.25;
	b=0.75;
	
	/*c=log(1.0/DAMPING)*(3.0*Vp/(2.0*FW));*/ /* Collino */
	c=-VPPML*log(alpha)/ffw;  /* Wang 2003 */
	
	for (i=1;i<=ifw;i++){
		   coeff[i]=c*((a*(((i-1)*DX)/ffw))+b*(((i-1)*DX)/ffw)*(((i-1)*DX)/ffw));  /* Wang 2003 */
		   	 
	}
	
	for (i=1;i<=(ifw/2);i++){
		dump=coeff[i];
		coeff[i]=coeff[ifw-i+1];
		coeff[ifw-i+1]=dump;	
	}
	
	if (MYID==0){
		fprintf(FP," Table of coefficients \n # \t coeff \n");
		/*fprintf(FP," ifw=%d \t a=%f amp=%f \n", ifw,a,amp); */
		for (i=1;i<=ifw;i++)
			fprintf(FP," %d \t %5.3f \n", i, coeff[i]);
	}	
	

	/* initialize array of coefficients with one */
	for (k=1;k<=NZ;k++)
	for (j=1;j<=NY;j++)
	for (i=1;i<=NX;i++) absorb_coeffx[j][i][k]=absorb_coeffy[j][i][k]=absorb_coeffz[j][i][k]=0.0;


	/* compute coefficients for left and right grid boundaries (x-direction) */
	if ((!BOUNDARY) && (POS[1]==0)){
		zb=1; yb=1; ye=NY; ze=NZ;
		for (i=1;i<=ifw;i++){
			if (POS[3]==0) zb=i; zb=1;
			if (POS[3]==NPROCZ-1) ze=NZ-i+1; ze=NZ;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i; yb=1;
			if (POS[2]==NPROCY-1) ye=NY-i+1; ye=NY;
			for (k=zb;k<=ze;k++)
				for (j=yb;j<=ye;j++)
					absorb_coeffx[j][i][k]=coeff[i];
		}
	}
			
	if ((!BOUNDARY) && (POS[1]==NPROCX-1)){
		zb=1; yb=1; ye=NY; ze=NZ;
		for (i=1;i<=ifw;i++){
			ii=NX-i+1;
			if (POS[3]==0) zb=i; zb=1;
			if (POS[3]==NPROCZ-1) ze=NZ-i+1; ze=NZ;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=i; yb=1;
			if (POS[2]==NPROCY-1) ye=NY-i+1; ye=NY;
			for (k=zb;k<=ze;k++)
				for (j=yb;j<=ye;j++)
					absorb_coeffx[j][ii][k]=coeff[i];
		}
	}
	
	/* compute coefficients for front and backward grid boundaries (z-direction) */
			
	if ((!BOUNDARY) && (POS[3]==0)){
		xb=1; yb=1; ye=NY; xe=NX;
		for (k=1;k<=ifw;k++){
			if (POS[1]==0) xb=k; xb=1;
			if (POS[1]==NPROCX-1) xe=NX-k+1; xe=NX;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=k; yb=1;
			if (POS[2]==NPROCY-1) ye=NY-k+1; ye=NY;
			for (i=xb;i<=xe;i++)
				for (j=yb;j<=ye;j++)
					absorb_coeffz[j][i][k]=coeff[k];
		}
	}

			
	if ((!BOUNDARY) && (POS[3]==NPROCZ-1)){
		xb=1; yb=1; ye=NY; xe=NX;
		for (k=1;k<=ifw;k++){
			kk=NZ-k+1;
			if (POS[1]==0) xb=k; xb=1;
			if (POS[1]==NPROCX-1) xe=NX-k+1; xe=NX;
			if ((POS[2]==0) && (!(FREE_SURF))) yb=k; yb=1;
			if (POS[2]==NPROCY-1) ye=NY-k+1; ye=NY;
			for (i=xb;i<=xe;i++)
				for (j=yb;j<=ye;j++)
					absorb_coeffz[j][i][kk]=coeff[k];
		}
	}

	/* compute coefficients for top and bottom grid boundaries (y-direction) */

	if ((POS[2]==0) && (!(FREE_SURF))){
		xb=1; zb=1; ze=NZ; xe=NX;
		for (j=1;j<=ifw;j++){
			if ((!BOUNDARY) && (POS[1]==0)) xb=j; xb=1;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1; xe=NX;
			if ((!BOUNDARY) && (POS[3]==0)) zb=j; zb=1;
			if ((!BOUNDARY) && (POS[3]==NPROCZ-1)) ze=NZ-j+1; ze=NZ;
			for (i=xb;i<=xe;i++)
				for (k=zb;k<=ze;k++)
					absorb_coeffy[j][i][k]=coeff[j];
		}
	}

	if (POS[2]==NPROCY-1){
		xb=1; zb=1; ze=NZ; xe=NX;
		for (j=1;j<=ifw;j++){
			jj=NY-j+1;
			if ((!BOUNDARY) && (POS[1]==0)) xb=j; xb=1;
			if ((!BOUNDARY) && (POS[1]==NPROCX-1)) xe=NX-j+1; xe=NX;
			if ((!BOUNDARY) && (POS[3]==0)) zb=j; zb=1;
		   if ((!BOUNDARY) && (POS[3]==NPROCZ-1)) ze=NZ-j+1; ze=NZ;
			for (i=xb;i<=xe;i++)
				for (k=zb;k<=ze;k++)
					absorb_coeffy[jj][i][k]=coeff[j];
		}
	}


	sprintf(modfile,"model/absorb.bin");

	writemod(modfile,absorb_coeffx,3); 

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3); 


	free_vector(coeff,1,ifw);
}



