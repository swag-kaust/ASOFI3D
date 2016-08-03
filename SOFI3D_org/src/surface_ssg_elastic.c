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
 *   stress free surface condition, elastic case
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_elastic(int ndepth, float *** u, float *** pi,
		float *** sxx, float ***syy, float ***szz, float *** sxy,float *** syz,
		float *** vx, float *** vy, float *** vz,
		float * K_x, float * a_x, float * b_x, float * K_z, float * a_z, float * b_z, 
		float *** psi_vxx, float *** psi_vzz ){

	int i, k ,j, fdoh,m,h1;
	float  vxx, vyy, vzz;
	float f, g, h; /* variables "dthalbe, dh24y, dh24z" removed, not in use */


	extern int NX, NZ, FDORDER, FDCOEFF, FW, NPROCX, NPROCZ, POS[4];
	register float b1, b2, b3, b4, b5, b6;
	extern float DT, DX, DY, DZ;
	extern int ABS_TYPE;
	j=ndepth;     /* The free surface is located exactly in y=(ndepth-1/2)*dh meter!! */

	//dthalbe=DT/2.0;
	fdoh=FDORDER/2;


	switch (FDORDER){
	case 2 :

		for (k=1;k<=NZ;k++){
			for (i=1;i<=NX;i++){


				/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
				syy[j][i][k]=0.0;

				syy[j-1][i][k]=-syy[j+1][i][k];
				sxy[j-1][i][k]=-sxy[j][i][k];
				sxy[j-2][i][k]=-sxy[j+1][i][k];
				syz[j-1][i][k]=-syz[j][i][k];
				syz[j-2][i][k]=-syz[j+1][i][k];

				vxx = (vx[j][i][k]-vx[j][i-1][k])/DX;
				vyy = (vy[j][i][k]-vy[j-1][i][k])/DY;
				vzz = (vz[j][i][k]-vz[j][i][k-1])/DZ;
				
				if (ABS_TYPE==1){
					if((POS[1]==0) && (i<=FW)){
						psi_vxx[j][i][k] = b_x[i] * psi_vxx[j][i][k] + a_x[i] * vxx;
						vxx = vxx / K_x[i] + psi_vxx[j][i][k];
					}
					if((POS[1]==NPROCX-1) && (i>=NX-FW+1)){
						h1 = i-NX+2*FW;

						psi_vxx[j][h1][k] = b_x[h1] * psi_vxx[j][h1][k] + a_x[h1] * vxx;
						vxx = vxx /K_x[h1] + psi_vxx[j][h1][k];
					}
					if((POS[3]==0) && (k<=FW)){
						psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
						vzz = vzz / K_z[k] + psi_vzz[j][i][k];					  
					}
					if((POS[3]==NPROCZ-1) && (k>=NZ-FW+1)){

						h1 = (k-NZ+2*FW);
						psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
						vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];
					}
					
				}

				f=u[j][i][k]*2.0;
				g=pi[j][i][k];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				sxx[j][i][k]+=h;
				szz[j][i][k]+=h;


			}
		}
		break;

	case 4 :

		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */

		if(FDCOEFF==2){
			b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/

		for (k=1;k<=NZ;k++){
			for (i=1;i<=NX;i++){


				/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
				syy[j][i][k]=0.0;

				for (m=1; m<=fdoh; m++) {
					syy[j-m][i][k]=-syy[j+m][i][k];
					sxy[j-m][i][k]=-sxy[j+m-1][i][k];
					syz[j-m][i][k]=-syz[j+m-1][i][k];
				}
				vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
				vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
				vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

				if (ABS_TYPE==1){
					if((POS[1]==0) && (i<=FW)){
						psi_vxx[j][i][k] = b_x[i] * psi_vxx[j][i][k] + a_x[i] * vxx;
						vxx = vxx / K_x[i] + psi_vxx[j][i][k];
					}
					if((POS[1]==NPROCX-1) && (i>=NX-FW+1)){
						h1 = i-NX+2*FW;

						psi_vxx[j][h1][k] = b_x[h1] * psi_vxx[j][h1][k] + a_x[h1] * vxx;
						vxx = vxx /K_x[h1] + psi_vxx[j][h1][k];
					}
					if((POS[3]==0) && (k<=FW)){
						psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
						vzz = vzz / K_z[k] + psi_vzz[j][i][k];					  
					}
					if((POS[3]==NPROCZ-1) && (k>=NZ-FW+1)){

						h1 = (k-NZ+2*FW);
						psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
						vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];
					}
					
				}
				
				f=u[j][i][k]*2.0;
				g=pi[j][i][k];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				sxx[j][i][k]+=h;
				szz[j][i][k]+=h;


			}
		}
		break;

	case 6 :

		b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/

		for (k=1;k<=NZ;k++){
			for (i=1;i<=NX;i++){


				/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
				syy[j][i][k]=0.0;

				for (m=1; m<=fdoh; m++) {
					syy[j-m][i][k]=-syy[j+m][i][k];
					sxy[j-m][i][k]=-sxy[j+m-1][i][k];
					syz[j-m][i][k]=-syz[j+m-1][i][k];
				}


				vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
						b2*(vx[j][i+1][k]-vx[j][i-2][k])+
						b3*(vx[j][i+2][k]-vx[j][i-3][k]))/DX;

				vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
						b2*(vy[j+1][i][k]-vy[j-2][i][k])+
						b3*(vy[j+2][i][k]-vy[j-3][i][k]))/DY;

				vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
						b2*(vz[j][i][k+1]-vz[j][i][k-2])+
						b3*(vz[j][i][k+2]-vz[j][i][k-3]))/DZ;


				f=u[j][i][k]*2.0;
				g=pi[j][i][k];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				sxx[j][i][k]+=h;
				szz[j][i][k]+=h;
			}
		}
		break;

	case 8 :

		b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients */

		if(FDCOEFF==2){
			b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/

		for (k=1;k<=NZ;k++){
			for (i=1;i<=NX;i++){


				/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
				syy[j][i][k]=0.0;

				for (m=1; m<=fdoh; m++) {
					syy[j-m][i][k]=-syy[j+m][i][k];
					sxy[j-m][i][k]=-sxy[j+m-1][i][k];
					syz[j-m][i][k]=-syz[j+m-1][i][k];
				}


				vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
						b2*(vx[j][i+1][k]-vx[j][i-2][k])+
						b3*(vx[j][i+2][k]-vx[j][i-3][k])+
						b4*(vx[j][i+3][k]-vx[j][i-4][k]))/DX;

				vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
						b2*(vy[j+1][i][k]-vy[j-2][i][k])+
						b3*(vy[j+2][i][k]-vy[j-3][i][k])+
						b4*(vy[j+3][i][k]-vy[j-4][i][k]))/DY;

				vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
						b2*(vz[j][i][k+1]-vz[j][i][k-2])+
						b3*(vz[j][i][k+2]-vz[j][i][k-3])+
						b4*(vz[j][i][k+3]-vz[j][i][k-4]))/DZ;


				f=u[j][i][k]*2.0;
				g=pi[j][i][k];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				sxx[j][i][k]+=h;
				szz[j][i][k]+=h;


			}
		}
		break;

	case 10 :

		b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor coefficients */

		if(FDCOEFF==2){
			b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/

		for (k=1;k<=NZ;k++){
			for (i=1;i<=NX;i++){


				/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
				syy[j][i][k]=0.0;

				for (m=1; m<=fdoh; m++) {
					syy[j-m][i][k]=-syy[j+m][i][k];
					sxy[j-m][i][k]=-sxy[j+m-1][i][k];
					syz[j-m][i][k]=-syz[j+m-1][i][k];
				}


				vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
						b2*(vx[j][i+1][k]-vx[j][i-2][k])+
						b3*(vx[j][i+2][k]-vx[j][i-3][k])+
						b4*(vx[j][i+3][k]-vx[j][i-4][k])+
						b5*(vx[j][i+4][k]-vx[j][i-5][k]))/DX;

				vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
						b2*(vy[j+1][i][k]-vy[j-2][i][k])+
						b3*(vy[j+2][i][k]-vy[j-3][i][k])+
						b4*(vy[j+3][i][k]-vy[j-4][i][k])+
						b5*(vy[j+4][i][k]-vy[j-5][i][k]))/DY;

				vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
						b2*(vz[j][i][k+1]-vz[j][i][k-2])+
						b3*(vz[j][i][k+2]-vz[j][i][k-3])+
						b4*(vz[j][i][k+3]-vz[j][i][k-4])+
						b5*(vz[j][i][k+4]-vz[j][i][k-5]))/DZ;

				f=u[j][i][k]*2.0;
				g=pi[j][i][k];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				sxx[j][i][k]+=h;
				szz[j][i][k]+=h;
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

		for (k=1;k<=NZ;k++){
			for (i=1;i<=NX;i++){


				/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
				syy[j][i][k]=0.0;

				for (m=1; m<=fdoh; m++) {
					syy[j-m][i][k]=-syy[j+m][i][k];
					sxy[j-m][i][k]=-sxy[j+m-1][i][k];
					syz[j-m][i][k]=-syz[j+m-1][i][k];
				}


				vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
						b2*(vx[j][i+1][k]-vx[j][i-2][k])+
						b3*(vx[j][i+2][k]-vx[j][i-3][k])+
						b4*(vx[j][i+3][k]-vx[j][i-4][k])+
						b5*(vx[j][i+4][k]-vx[j][i-5][k])+
						b6*(vx[j][i+5][k]-vx[j][i-6][k]))/DX;

				vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
						b2*(vy[j+1][i][k]-vy[j-2][i][k])+
						b3*(vy[j+2][i][k]-vy[j-3][i][k])+
						b4*(vy[j+3][i][k]-vy[j-4][i][k])+
						b5*(vy[j+4][i][k]-vy[j-5][i][k])+
						b6*(vy[j+5][i][k]-vy[j-6][i][k]))/DY;

				vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
						b2*(vz[j][i][k+1]-vz[j][i][k-2])+
						b3*(vz[j][i][k+2]-vz[j][i][k-3])+
						b4*(vz[j][i][k+3]-vz[j][i][k-4])+
						b5*(vz[j][i][k+4]-vz[j][i][k-5])+
						b6*(vz[j][i][k+5]-vz[j][i][k-6]))/DZ;

				f=u[j][i][k]*2.0;
				g=pi[j][i][k];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				sxx[j][i][k]+=h;
				szz[j][i][k]+=h;
			}
		}
		break;
	}
}
