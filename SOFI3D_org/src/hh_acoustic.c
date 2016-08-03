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
/*----------------------------------------------------------------------
 *   homogeneous acoustic model
 ---------------------------------------------------------------------- */

#include "fd.h"

void model_acoustic(float  ***  rho, float ***  pi){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], MYID;
	extern char  MFILE[STRING_SIZE];
	extern int WRITE_MODELFILES;

	/* local variables */
	float piv;
	float Vp, Rhov;
	int i, j, k, ii, jj, kk;
	char modfile[STRING_SIZE];

	/*-----------------------------------------------------------------------*/




	for (j=1;j<=NYG;j++){
		for (i=1;i<=NXG;i++){
			for (k=1;k<=NZG;k++){

				/*=========================================================
				 * modify below this point for ACOUSTIC model definition
				 *=========================================================
				 */


				Vp=3500.0;  Rhov=2000.0;

				/*=========================================================
				 * modify up to this this point for ACOUSTIC model definition
				 *=========================================================
				 */

				piv=Vp*Vp*Rhov;

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
						(POS[2]==((j-1)/NY)) &&
						(POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					rho[jj][ii][kk]=Rhov;
					pi[jj][ii][kk]=piv;
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
}



