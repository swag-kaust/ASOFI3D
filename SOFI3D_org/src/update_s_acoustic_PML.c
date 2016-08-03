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
 *   updating stress at gridpoints [nx1...nx2][ny1...ny2][nz1...nz2]
 *   by a staggered grid finite difference scheme of 4th order accuracy in space
 *   and second order accuracy in time
 *   acoustic version
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s_acoustic_PML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
		float *** vx, float *** vy, float *** vz, float *** sxx,
		float *** sxx1, float *** sxx2, float *** sxx3, float ***  pi, float *** absorb_coeffx,
		float *** absorb_coeffy, float *** absorb_coeffz){


	extern float DT, DX, DY, DZ;
	extern int MYID, FDORDER, FDCOEFF;
	extern FILE *FP;
	/*float DH;*/
	extern int OUTNTIMESTEPINFO;

	register int i, j, k;
	double time=0.0, time1=0.0, time2=0.0;
	register float vxx, vyy, vzz, PML1, PML2;
	register float b1, b2, b3, b4, b5, b6, dx, dy, dz;

	dx=1.0/DX;
	dy=1.0/DY;
	dz=1.0/DZ;

	if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) time1=MPI_Wtime();

	switch (FDORDER){
	case 2 :


		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){				

					vxx=(vx[j][i][k]-vx[j][i-1][k]);
					vyy=(vy[j][i][k]-vy[j-1][i][k]);
					vzz=(vz[j][i][k]-vz[j][i][k-1]);

					PML1=((((1.0/DT)-(absorb_coeffx[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));

					sxx1[j][i][k] = (PML1*sxx1[j][i][k]) + PML2*(dx*pi[j][i][k]*vxx);


					PML1=((((1.0/DT)-(absorb_coeffy[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));

					sxx2[j][i][k] = (PML1*sxx2[j][i][k]) + PML2*(dy*pi[j][i][k]*vyy);


					PML1=((((1.0/DT)-(absorb_coeffz[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));

					sxx3[j][i][k] = (PML1*sxx3[j][i][k]) + PML2*(dz*pi[j][i][k]*vzz);


					sxx[j][i][k] = sxx1[j][i][k] + sxx2[j][i][k] + sxx3[j][i][k];

				}
			}
		}
		break;



	case 4 :

		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */ 
		if(FDCOEFF==2){
			b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/


			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					for (k=nz1;k<=nz2;k++){


						vxx=b1*(vx[j][i][k]-vx[j][i-1][k])+
								b2*(vx[j][i+1][k]-vx[j][i-2][k]);
						vyy=b1*(vy[j][i][k]-vy[j-1][i][k])+
								b2*(vy[j+1][i][k]-vy[j-2][i][k]);
						vzz=b1*(vz[j][i][k]-vz[j][i][k-1])+
								b2*(vz[j][i][k+1]-vz[j][i][k-2]);

						PML1=((((1.0/DT)-(absorb_coeffx[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));
						PML2=(1.0/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));

						sxx1[j][i][k] = (PML1*sxx1[j][i][k]) + PML2*(dx*pi[j][i][k]*vxx);


						PML1=((((1.0/DT)-(absorb_coeffy[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));
						PML2=(1.0/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));

						sxx2[j][i][k] = (PML1*sxx2[j][i][k]) + PML2*(dy*pi[j][i][k]*vyy);


						PML1=((((1.0/DT)-(absorb_coeffz[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));
						PML2=(1.0/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));

						sxx3[j][i][k] = (PML1*sxx3[j][i][k]) + PML2*(dz*pi[j][i][k]*vzz);


						sxx[j][i][k] = sxx1[j][i][k] + sxx2[j][i][k] + sxx3[j][i][k];

					}
				}
			}
			break;

	case 6 :
		b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/


		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){


					vxx=b1*(vx[j][i][k]-vx[j][i-1][k])+
							b2*(vx[j][i+1][k]-vx[j][i-2][k])+
							b3*(vx[j][i+2][k]-vx[j][i-3][k]);
					vyy=b1*(vy[j][i][k]-vy[j-1][i][k])+
							b2*(vy[j+1][i][k]-vy[j-2][i][k])+
							b3*(vy[j+2][i][k]-vy[j-3][i][k]);
					vzz=b1*(vz[j][i][k]-vz[j][i][k-1])+
							b2*(vz[j][i][k+1]-vz[j][i][k-2])+
							b3*(vz[j][i][k+2]-vz[j][i][k-3]);

					PML1=((((1.0/DT)-(absorb_coeffx[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));

					sxx1[j][i][k] = (PML1*sxx1[j][i][k]) + PML2*(dx*pi[j][i][k]*vxx);


					PML1=((((1.0/DT)-(absorb_coeffy[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));

					sxx2[j][i][k] = (PML1*sxx2[j][i][k]) + PML2*(dy*pi[j][i][k]*vyy);


					PML1=((((1.0/DT)-(absorb_coeffz[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));

					sxx3[j][i][k] = (PML1*sxx3[j][i][k]) + PML2*(dz*pi[j][i][k]*vzz);


					sxx[j][i][k] = sxx1[j][i][k] + sxx2[j][i][k] + sxx3[j][i][k];


				}
			}
		}
		break;
	case 8 :
		b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/


		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){


					vxx=b1*(vx[j][i][k]-vx[j][i-1][k])+
							b2*(vx[j][i+1][k]-vx[j][i-2][k])+
							b3*(vx[j][i+2][k]-vx[j][i-3][k])+
							b4*(vx[j][i+3][k]-vx[j][i-4][k]);
					vyy=b1*(vy[j][i][k]-vy[j-1][i][k])+
							b2*(vy[j+1][i][k]-vy[j-2][i][k])+
							b3*(vy[j+2][i][k]-vy[j-3][i][k])+
							b4*(vy[j+3][i][k]-vy[j-4][i][k]);
					vzz=b1*(vz[j][i][k]-vz[j][i][k-1])+
							b2*(vz[j][i][k+1]-vz[j][i][k-2])+
							b3*(vz[j][i][k+2]-vz[j][i][k-3])+
							b4*(vz[j][i][k+3]-vz[j][i][k-4]);

					PML1=((((1.0/DT)-(absorb_coeffx[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));

					sxx1[j][i][k] = (PML1*sxx1[j][i][k]) + PML2*(dx*pi[j][i][k]*vxx);


					PML1=((((1.0/DT)-(absorb_coeffy[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));

					sxx2[j][i][k] = (PML1*sxx2[j][i][k]) + PML2*(dy*pi[j][i][k]*vyy);


					PML1=((((1.0/DT)-(absorb_coeffz[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));

					sxx3[j][i][k] = (PML1*sxx3[j][i][k]) + PML2*(dz*pi[j][i][k]*vzz);


					sxx[j][i][k] = sxx1[j][i][k] + sxx2[j][i][k] + sxx3[j][i][k];


				}
			}
		}
		break;

	case 10 :
		b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/


		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){


					vxx=b1*(vx[j][i][k]-vx[j][i-1][k])+
							b2*(vx[j][i+1][k]-vx[j][i-2][k])+
							b3*(vx[j][i+2][k]-vx[j][i-3][k])+
							b4*(vx[j][i+3][k]-vx[j][i-4][k])+
							b5*(vx[j][i+4][k]-vx[j][i-5][k]);
					vyy=b1*(vy[j][i][k]-vy[j-1][i][k])+
							b2*(vy[j+1][i][k]-vy[j-2][i][k])+
							b3*(vy[j+2][i][k]-vy[j-3][i][k])+
							b4*(vy[j+3][i][k]-vy[j-4][i][k])+
							b5*(vy[j+4][i][k]-vy[j-5][i][k]);
					vzz=b1*(vz[j][i][k]-vz[j][i][k-1])+
							b2*(vz[j][i][k+1]-vz[j][i][k-2])+
							b3*(vz[j][i][k+2]-vz[j][i][k-3])+
							b4*(vz[j][i][k+3]-vz[j][i][k-4])+
							b5*(vz[j][i][k+4]-vz[j][i][k-5]);

					PML1=((((1.0/DT)-(absorb_coeffx[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));

					sxx1[j][i][k] = (PML1*sxx1[j][i][k]) + PML2*(dx*pi[j][i][k]*vxx);


					PML1=((((1.0/DT)-(absorb_coeffy[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));

					sxx2[j][i][k] = (PML1*sxx2[j][i][k]) + PML2*(dy*pi[j][i][k]*vyy);


					PML1=((((1.0/DT)-(absorb_coeffz[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));

					sxx3[j][i][k] = (PML1*sxx3[j][i][k]) + PML2*(dz*pi[j][i][k]*vzz);


					sxx[j][i][k] = sxx1[j][i][k] + sxx2[j][i][k] + sxx3[j][i][k];


				}
			}
		}
		break;

	case 12 :

		/* Taylor coefficients */
		b1=160083.0/131072.0; b2=-12705.0/131072.0; b3=22869.0/1310720.0; 
		b4=-5445.0/1835008.0; b5=847.0/2359296.0; b6=-63.0/2883584;

		/* Holberg coefficients E=0.1 %*/
		if(FDCOEFF==2){
			b1=1.2508; b2=-0.12034; b3=0.032131; b4=-0.010142; b5=0.0029857; b6=-0.00066667;}

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){


					vxx=b1*(vx[j][i][k]-vx[j][i-1][k])+
							b2*(vx[j][i+1][k]-vx[j][i-2][k])+
							b3*(vx[j][i+2][k]-vx[j][i-3][k])+
							b4*(vx[j][i+3][k]-vx[j][i-4][k])+
							b5*(vx[j][i+4][k]-vx[j][i-5][k])+
							b6*(vx[j][i+5][k]-vx[j][i-6][k]);
					vyy=b1*(vy[j][i][k]-vy[j-1][i][k])+
							b2*(vy[j+1][i][k]-vy[j-2][i][k])+
							b3*(vy[j+2][i][k]-vy[j-3][i][k])+
							b4*(vy[j+3][i][k]-vy[j-4][i][k])+
							b5*(vy[j+4][i][k]-vy[j-5][i][k])+
							b6*(vy[j+5][i][k]-vy[j-6][i][k]);
					vzz=b1*(vz[j][i][k]-vz[j][i][k-1])+
							b2*(vz[j][i][k+1]-vz[j][i][k-2])+
							b3*(vz[j][i][k+2]-vz[j][i][k-3])+
							b4*(vz[j][i][k+3]-vz[j][i][k-4])+
							b5*(vz[j][i][k+4]-vz[j][i][k-5])+
							b6*(vz[j][i][k+5]-vz[j][i][k-6]);

					PML1=((((1.0/DT)-(absorb_coeffx[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffx[j][i][k]/2.0))));

					sxx1[j][i][k] = (PML1*sxx1[j][i][k]) + PML2*(dx*pi[j][i][k]*vxx);


					PML1=((((1.0/DT)-(absorb_coeffy[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffy[j][i][k]/2.0))));

					sxx2[j][i][k] = (PML1*sxx2[j][i][k]) + PML2*(dy*pi[j][i][k]*vyy);


					PML1=((((1.0/DT)-(absorb_coeffz[j][i][k]/2.0)))/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));
					PML2=(1.0/(((1.0/DT)+(absorb_coeffz[j][i][k]/2.0))));

					sxx3[j][i][k] = (PML1*sxx3[j][i][k]) + PML2*(dz*pi[j][i][k]*vzz);


					sxx[j][i][k] = sxx1[j][i][k] + sxx2[j][i][k] + sxx3[j][i][k];


				}
			}
		}
		break;


	default:
		err(" error in particle velocity update: wrong FDORDER. ");
		break;
	} /* end of switch (FDORDER) */


	if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0){
		time2=MPI_Wtime();
		time=time2-time1;
		fprintf(FP," Real time for stress tensor update: \t\t %4.2f s.\n",time);
	}
	return time;

}
