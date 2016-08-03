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
 * Averaging of material parameters
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void av_mat(float *** rho, float *** pi, float *** u,
		float *** taus, float *** taup,
		float  *** uipjp, float *** ujpkp, float *** uipkp, float *** tausipjp,
		float  *** tausjpkp, float  *** tausipkp, float  *** rjp, float  *** rkp, float  *** rip ){


	extern int NX, NY, NZ, MYID, L;
	extern FILE *FP;
	double time1=0.0, time2=0.0;
	int i, j, k;

	if (MYID==0){
		fprintf(FP,"\n\n **Message from av_mat (written by PE %d):",MYID);
		fprintf(FP,"\n Averaging of material parameters ... \n");
		time1=MPI_Wtime();
	}	

	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){


				/* harmonic averaging of shear modulus */
				uipjp[j][i][k]=4.0/((1.0/u[j][i][k])+(1.0/u[j][i+1][k])+(1.0/u[j+1][i+1][k])+(1.0/u[j+1][i][k]));
				ujpkp[j][i][k]=4.0/((1.0/u[j][i][k])+(1.0/u[j][i][k+1])+(1.0/u[j+1][i][k+1])+(1.0/u[j+1][i][k]));
				uipkp[j][i][k]=4.0/((1.0/u[j][i][k])+(1.0/u[j][i][k+1])+(1.0/u[j][i+1][k+1])+(1.0/u[j][i+1][k]));

				/* arithmetic averaging of TAU for S-waves and density */

				if (L){
					tausipjp[j][i][k]=0.25*(taus[j][i][k]+taus[j][i+1][k]+taus[j+1][i+1][k]+taus[j+1][i][k]);
					tausjpkp[j][i][k]=0.25*(taus[j][i][k]+taus[j+1][i][k]+taus[j+1][i][k+1]+taus[j][i][k+1]);
					tausipkp[j][i][k]=0.25*(taus[j][i][k]+taus[j][i+1][k]+taus[j][i+1][k+1]+taus[j][i][k+1]);
				}

				rjp[j][i][k]=0.5*(rho[j][i][k]+rho[j+1][i][k]);
				rkp[j][i][k]=0.5*(rho[j][i][k]+rho[j][i][k+1]);
				rip[j][i][k]=0.5*(rho[j][i][k]+rho[j][i+1][k]);


			}
		}
	}



	if (MYID==0){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}



}




















