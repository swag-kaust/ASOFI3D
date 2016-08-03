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
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2][nz1...nz2]
 *   by a staggered grid finite difference scheme of FDORDER order accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_v_acoustic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz, float *** sxx, float  ***  rho, float **  srcpos_loc, float ** signals, int nsrc,
		float *** absorb_coeff, int * stype){


	extern float DT, DX, DY, DZ, SOURCE_ALPHA, SOURCE_BETA;
	extern int MYID, FDORDER, ABS_TYPE, LOG, FDCOEFF;
	extern int OUTNTIMESTEPINFO;

	extern FILE *FP;

	register int i, j, k, l;
	float  amp, alpha_rad, beta_rad, rjp, rkp, rip;
	double time=0.0, time1=0.0, time2=0.0;
	register float b1, b2, b3, b4, b5, b6, dx, dy, dz;
	register float sxx_x, syy_y, szz_z;

	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) time1=MPI_Wtime();

	switch (FDORDER){
	case 2 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){

					rjp=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
					rkp=0.5*(rho[j][i][k]+rho[j][i][k+1]);
					rip=0.5*(rho[j][i][k]+rho[j][i+1][k]);

					vx[j][i][k]+= dx*(sxx[j][i+1][k]-sxx[j][i][k])/rip;
					vy[j][i][k]+= dy*(sxx[j+1][i][k]-sxx[j][i][k])/rjp;
					vz[j][i][k]+= dz*(sxx[j][i][k+1]-sxx[j][i][k])/rkp;

				}
			}
		}
		break;

	case 4 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;

		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){

					rjp=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
					rkp=0.5*(rho[j][i][k]+rho[j][i][k+1]);
					rip=0.5*(rho[j][i][k]+rho[j][i+1][k]);

					sxx_x = b1*(sxx[j][i+1][k]-sxx[j][i][k])+
							b2*(sxx[j][i+2][k]-sxx[j][i-1][k]);
					syy_y =	b1*(sxx[j+1][i][k]-sxx[j][i][k])+
							b2*(sxx[j+2][i][k]-sxx[j-1][i][k]);
					szz_z = b1*(sxx[j][i][k+1]-sxx[j][i][k])+
							b2*(sxx[j][i][k+2]-sxx[j][i][k-1]);

					vx[j][i][k]+=(dx*sxx_x/rip);
					vy[j][i][k]+=(dy*syy_y/rjp);
					vz[j][i][k]+=(dz*szz_z/rkp);

				}
			}
		}
		break;

	case 6 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;

		b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){

					rjp=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
					rkp=0.5*(rho[j][i][k]+rho[j][i][k+1]);
					rip=0.5*(rho[j][i][k]+rho[j][i+1][k]);

					sxx_x = b1*(sxx[j][i+1][k]-sxx[j][i][k])+
							b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
							b3*(sxx[j][i+3][k]-sxx[j][i-2][k]);
					syy_y =	b1*(sxx[j+1][i][k]-sxx[j][i][k])+
							b2*(sxx[j+2][i][k]-sxx[j-1][i][k])+
							b3*(sxx[j+3][i][k]-sxx[j-2][i][k]);
					szz_z = b1*(sxx[j][i][k+1]-sxx[j][i][k])+
							b2*(sxx[j][i][k+2]-sxx[j][i][k-1])+
							b3*(sxx[j][i][k+3]-sxx[j][i][k-2]);

					vx[j][i][k]+=(dx*sxx_x/rip);
					vy[j][i][k]+=(dy*syy_y/rjp);
					vz[j][i][k]+=(dz*szz_z/rkp);

				}
			}
		}
		break;

	case 8 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;

		b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){

					rjp=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
					rkp=0.5*(rho[j][i][k]+rho[j][i][k+1]);
					rip=0.5*(rho[j][i][k]+rho[j][i+1][k]);

					sxx_x = b1*(sxx[j][i+1][k]-sxx[j][i][k])+
							b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
							b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
							b4*(sxx[j][i+4][k]-sxx[j][i-3][k]);
					syy_y =	b1*(sxx[j+1][i][k]-sxx[j][i][k])+
							b2*(sxx[j+2][i][k]-sxx[j-1][i][k])+
							b3*(sxx[j+3][i][k]-sxx[j-2][i][k])+
							b4*(sxx[j+4][i][k]-sxx[j-3][i][k]);
					szz_z = b1*(sxx[j][i][k+1]-sxx[j][i][k])+
							b2*(sxx[j][i][k+2]-sxx[j][i][k-1])+
							b3*(sxx[j][i][k+3]-sxx[j][i][k-2])+
							b4*(sxx[j][i][k+4]-sxx[j][i][k-3]);

					vx[j][i][k]+=(dx*sxx_x/rip);
					vy[j][i][k]+=(dy*syy_y/rjp);
					vz[j][i][k]+=(dz*szz_z/rkp);

				}
			}
		}
		break;

	case 10 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;

		b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){

					rjp=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
					rkp=0.5*(rho[j][i][k]+rho[j][i][k+1]);
					rip=0.5*(rho[j][i][k]+rho[j][i+1][k]);

					sxx_x = b1*(sxx[j][i+1][k]-sxx[j][i][k])+
							b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
							b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
							b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
							b5*(sxx[j][i+5][k]-sxx[j][i-4][k]);
					syy_y =	b1*(sxx[j+1][i][k]-sxx[j][i][k])+
							b2*(sxx[j+2][i][k]-sxx[j-1][i][k])+
							b3*(sxx[j+3][i][k]-sxx[j-2][i][k])+
							b4*(sxx[j+4][i][k]-sxx[j-3][i][k])+
							b5*(sxx[j+5][i][k]-sxx[j-4][i][k]);
					szz_z = b1*(sxx[j][i][k+1]-sxx[j][i][k])+
							b2*(sxx[j][i][k+2]-sxx[j][i][k-1])+
							b3*(sxx[j][i][k+3]-sxx[j][i][k-2])+
							b4*(sxx[j][i][k+4]-sxx[j][i][k-3])+
							b5*(sxx[j][i][k+5]-sxx[j][i][k-4]);

					vx[j][i][k]+=(dx*sxx_x/rip);
					vy[j][i][k]+=(dy*syy_y/rjp);
					vz[j][i][k]+=(dz*szz_z/rkp);

				}
			}
		}
		break;

	case 12 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;


		/* Taylor coefficients */
		b1=160083.0/131072.0; b2=-12705.0/131072.0; b3=22869.0/1310720.0; 
		b4=-5445.0/1835008.0; b5=847.0/2359296.0; b6=-63.0/2883584;

		/* Holberg coefficients E=0.1 %*/
		if(FDCOEFF==2){
			b1=1.2508; b2=-0.12034; b3=0.032131; b4=-0.010142; b5=0.0029857; b6=-0.00066667;}

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){

					rjp=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
					rkp=0.5*(rho[j][i][k]+rho[j][i][k+1]);
					rip=0.5*(rho[j][i][k]+rho[j][i+1][k]);

					sxx_x = b1*(sxx[j][i+1][k]-sxx[j][i][k])+
							b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
							b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
							b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
							b5*(sxx[j][i+5][k]-sxx[j][i-4][k])+
							b6*(sxx[j][i+6][k]-sxx[j][i-5][k]);
					syy_y =	b1*(sxx[j+1][i][k]-sxx[j][i][k])+
							b2*(sxx[j+2][i][k]-sxx[j-1][i][k])+
							b3*(sxx[j+3][i][k]-sxx[j-2][i][k])+
							b4*(sxx[j+4][i][k]-sxx[j-3][i][k])+
							b5*(sxx[j+5][i][k]-sxx[j-4][i][k])+
							b6*(sxx[j+6][i][k]-sxx[j-5][i][k]);
					szz_z = b1*(sxx[j][i][k+1]-sxx[j][i][k])+
							b2*(sxx[j][i][k+2]-sxx[j][i][k-1])+
							b3*(sxx[j][i][k+3]-sxx[j][i][k-2])+
							b4*(sxx[j][i][k+4]-sxx[j][i][k-3])+
							b5*(sxx[j][i][k+5]-sxx[j][i][k-4])+
							b6*(sxx[j][i][k+6]-sxx[j][i][k-5]);

					vx[j][i][k]+=(dx*sxx_x/rip);
					vy[j][i][k]+=(dy*syy_y/rjp);
					vz[j][i][k]+=(dz*szz_z/rkp);

				}
			}
		}
		break;


	default:
		err(" error in particle velocity update: wrong FDORDER. ");
		break;
	} /* end of switch (FDORDER) */


	/* Adding body force components to corresponding particle velocities */
	for (l=1;l<=nsrc;l++) {
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		k=(int)srcpos_loc[3][l];

		//amp=signals[l][nt]; // unscaled force amplitude
		amp=(DT*rho[j][i][k]*signals[l][nt])/(DX*DY*DZ);// scaled force amplitude with F= 1N


		switch (stype[l]){
		case 2 : 
			vx[j][i][k]+=amp;  /* single force in x */
			break;
		case 3 : 
			vy[j][i][k]+=amp;  /* single force in y */
			/*note that "y" is used for the vertical coordinate */
			break;
		case 4 : 
			vz[j][i][k]+=amp;  /* single force in z */
			break;
		case 5 : 

			alpha_rad=SOURCE_ALPHA*PI/180; /* custom force */
			beta_rad=SOURCE_BETA*PI/180;
			vx[j][i][k]+=cos(alpha_rad)*sin(beta_rad)*amp;
			vy[j][i][k]+=cos(beta_rad)*amp; /*vertical component*/
			vz[j][i][k]+=sin(alpha_rad)*sin(beta_rad)*amp;
			break;
		}
	}






	/* absorbing boundary condition (exponential damping) */

	if (ABS_TYPE==2){
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
					vx[j][i][k]*=absorb_coeff[j][i][k];
					vy[j][i][k]*=absorb_coeff[j][i][k];
					vz[j][i][k]*=absorb_coeff[j][i][k];
					sxx[j][i][k]*=absorb_coeff[j][i][k];

				}
			}
		}
	}

	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {

			time2=MPI_Wtime();
			time=time2-time1;
			fprintf(FP," Real time for particle velocity update: \t %4.2f s.\n",time);
		}
	return time;

}


