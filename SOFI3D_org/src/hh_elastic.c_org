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
 *   Generates elastic model properties (vp,vs,density) on the fly
 *
 *   depending on model dimension in vertical direction and local variable "h"
 *   this function can generate a
 *   	-> homogeneneous full space
 *   	-> layer over half space
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void model_elastic(float  ***  rho, float ***  pi, float ***  u,
		float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float DY;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	extern int WRITE_MODELFILES;

	/* local variables */
	float muv, piv;
	float Vp, Vs, Rho;
	float *** pwavemod=NULL, *** swavemod=NULL;
	float y;
	int i, j, k, ii, jj, kk;
	char modfile[STRING_SIZE];

	/*-----------------material property definition -------------------------*/

	/* parameters for layer 1 */
	const float vp1=3500.0, vs1=2000.0, rho1=2000.0, h=100000.0;

	/* parameters for layer 2 */
	const float vp2=5700.0, vs2=3400.0, rho2=2500.0;


	if (WRITE_MODELFILES==1) {
		pwavemod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		swavemod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	}

	/*elastic simulation */
	if (L==0) {
		/* loop over global grid */
		for (k=1;k<=NZG;k++){
			for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){

					/*=========================================================
					 * modify below this point for ELASTIC model definition
					 *=========================================================
					 */
					/*note that "y" is used for the vertical coordinate*/
					/* calculate vertical coordinate in m */

					y=(float)j*DY;

					/* two layer case */
					if (y<=h){
						Vp=vp1; Vs=vs1; Rho=rho1; }


					else{
						Vp=vp2; Vs=vs2; Rho=rho2;}

					/*=========================================================
					 * modify up to this point for ELASTIC model definition
					 *=========================================================
					 */

					muv=Vs*Vs*Rho;
					piv=Vp*Vp*Rho;

					/* only the PE which belongs to the current global gridpoint
							is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) &&
							(POS[2]==((j-1)/NY)) &&
							(POS[3]==((k-1)/NZ))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						kk=k-POS[3]*NZ;

						u[jj][ii][kk]=muv;
						rho[jj][ii][kk]=Rho;
						pi[jj][ii][kk]=piv;

						if (WRITE_MODELFILES==1) {
							pwavemod[jj][ii][kk]=Vp;
							swavemod[jj][ii][kk]=Vs;
						}
					}
				}
			}
		}
	}


	/* each PE writes his model to disk */

	/* all models are written to file */
	if (WRITE_MODELFILES==1) {
		sprintf(modfile,"%s.SOFI3D.pi",MFILE);
		writemod(modfile,pi,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.u",MFILE);
		writemod(modfile,u,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.vp",MFILE);
		writemod(modfile,pwavemod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.vs",MFILE);
		writemod(modfile,swavemod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	/* only density is written to file */
	if (WRITE_MODELFILES==2) {
		sprintf(modfile,"%s.SOFI3D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}


	if (WRITE_MODELFILES==1) {
		free_f3tensor(pwavemod,0,NY+1,0,NX+1,0,NZ+1);
		free_f3tensor(swavemod,0,NY+1,0,NX+1,0,NZ+1);
	}

}
