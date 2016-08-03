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
 * Implementation acoording to Komatitsch, D. and Martin, R.(2007): "An unsplit convolutional perfectly matched
 * layer improved at grazing incidence for the seismic wave equation", geophysics, Vol.72, No.5
 * similar to fdveps (2D) 
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s_CPML_elastic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
		float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz, float ***  pi, float ***  u,
		float ***  uipjp, float ***  ujpkp, float ***  uipkp,
		float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		float * K_z, float * a_z, float * b_z, float * K_z_half, float * a_z_half, float * b_z_half,
		float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy, float *** psi_vxz, float *** psi_vyz, float *** psi_vzz){


	extern float DT, DX, DY, DZ;
	extern int MYID, LOG, FDCOEFF, FDORDER;
	extern FILE *FP;
	extern int FREE_SURF;
	extern int NPROCX, NPROCY, NPROCZ, POS[4];
	extern int FW, NY, NZ;
	extern int OUTNTIMESTEPINFO;

	int i, j, k, h1;
	double time=0.0, time1=0.0, time2=0.0;
	float vxx=0.0,vxy=0.0,vxz=0.0,vyx=0.0,vyy=0.0,vyz=0.0,vzx=0.0,vzy=0.0,vzz=0.0;
	float vxyyx,vyzzy,vxzzx,vxxyyzz,vyyzz,vxxzz,vxxyy;
	float g,f,fipjp,fjpkp,fipkp;
	float b1=1.0, b2=0.0;

	//float dthalbe;
	//dthalbe = DT/2.0;

	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) time1=MPI_Wtime();


	if(FDCOEFF==1){ /* Taylor coefficients*/

		switch (FDORDER){
		case 2 :
			b1=1.0; b2=0.0;
			break;
		case 4 :
			b1=9.0/8.0; b2=-1.0/24.0;
			break;
		default:
			b1=9.0/8.0; b2=-1.0/24.0;
			//warning(" Warning in function update_s_ssg_CPML_elastic : CPML Update of FDORDER > 4 not implemented yet, coefficients of FDORDER==4 are used instead!");
			//warning is issued in checkfd_ssg, section "ABSORBING BOUNDARY"
			break;
		}
	}

	if(FDCOEFF==2){ /* Holberg coefficients E=0.1 %*/

		switch (FDORDER){
		case 2 :
			b1=1.0; b2=0.0;
			break;
		case 4 :
			b1=1.1382; b2=-0.046414;
			break;
		default:
			b1=1.1382; b2=-0.046414;
			//warning(" Warning in function update_s_ssg_CPML_elastic : CPML Update of FDORDER > 4 not implemented yet, coefficients of FDORDER==4 are used instead!");
			//warning is issued in checkfd_ssg, section "ABSORBING BOUNDARY"
			break;
		}
	}



	if (POS[1]==0){
		for (j=1;j<=NY;j++){
			for (i=1;i<=FW;i++){
				for (k=1;k<=NZ;k++){

					vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
					vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
					vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
					vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
					vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
					vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
					vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
					vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
					vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

					psi_vxx[j][i][k] = b_x[i] * psi_vxx[j][i][k] + a_x[i] * vxx;
					vxx = vxx / K_x[i] + psi_vxx[j][i][k];
					psi_vyx[j][i][k] = b_x_half[i] * psi_vyx[j][i][k] + a_x_half[i] * vyx;
					vyx = vyx / K_x_half[i] + psi_vyx[j][i][k];
					psi_vzx[j][i][k] = b_x_half[i] * psi_vzx[j][i][k] + a_x_half[i] * vzx;
					vzx = vzx / K_x_half[i] + psi_vzx[j][i][k];


					if((POS[2]==0 && FREE_SURF==0) && (j<=FW)){

						psi_vxy[j][i][k] = b_y_half[j] * psi_vxy[j][i][k] + a_y_half[j] * vxy;
						vxy = vxy / K_y_half[j] + psi_vxy[j][i][k];
						psi_vyy[j][i][k] = b_y[j] * psi_vyy[j][i][k] + a_y[j] * vyy;
						vyy = vyy / K_y[j] + psi_vyy[j][i][k];
						psi_vzy[j][i][k] = b_y_half[j] * psi_vzy[j][i][k] + a_y_half[j] * vzy;
						vzy = vzy / K_y_half[j] + psi_vzy[j][i][k]; }

					if((POS[2]==NPROCY-1) && (j>=ny2+1)){
						h1 = (j-ny2+FW);

						psi_vxy[h1][i][k] = b_y_half[h1] * psi_vxy[h1][i][k] + a_y_half[h1] * vxy;
						vxy = vxy / K_y_half[h1] + psi_vxy[h1][i][k];
						psi_vyy[h1][i][k] = b_y[h1] * psi_vyy[h1][i][k] + a_y[h1] * vyy;
						vyy = vyy / K_y[h1] + psi_vyy[h1][i][k];
						psi_vzy[h1][i][k] = b_y_half[h1] * psi_vzy[h1][i][k] + a_y_half[h1] * vzy;
						vzy = vzy / K_y_half[h1] + psi_vzy[h1][i][k]; }

					if((POS[3]==0) && (k<=FW)){

						psi_vxz[j][i][k] = b_z_half[k] * psi_vxz[j][i][k] + a_z_half[k] * vxz;
						vxz = vxz / K_z_half[k] + psi_vxz[j][i][k];
						psi_vyz[j][i][k] = b_z_half[k] * psi_vyz[j][i][k] + a_z_half[k] * vyz;
						vyz = vyz / K_z_half[k] + psi_vyz[j][i][k];
						psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
						vzz = vzz / K_z[k] + psi_vzz[j][i][k];}

					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){

						h1 = (k-nz2+FW);

						psi_vxz[j][i][h1] = b_z_half[h1] * psi_vxz[j][i][h1] + a_z_half[h1] * vxz;
						vxz = vxz / K_z_half[h1] + psi_vxz[j][i][h1];
						psi_vyz[j][i][h1] = b_z_half[h1] * psi_vyz[j][i][h1] + a_z_half[h1] * vyz;
						vyz = vyz / K_z_half[h1] + psi_vyz[j][i][h1];
						psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
						vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];}




					fipjp=uipjp[j][i][k]*DT;
					fjpkp=ujpkp[j][i][k]*DT;
					fipkp=uipkp[j][i][k]*DT;
					g=pi[j][i][k];
					f=2.0*u[j][i][k];

					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;

					sxy[j][i][k]+=(fipjp*vxyyx);
					syz[j][i][k]+=(fjpkp*vyzzy);
					sxz[j][i][k]+=(fipkp*vxzzx);
					sxx[j][i][k]+=DT*((g*vxxyyzz)-(f*vyyzz));
					syy[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxzz));
					szz[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxyy));

				}
			}
		}
	}

	if(POS[1]==NPROCX-1){
		for (j=1;j<=NY;j++){
			for (i=nx2+1;i<=nx2+FW;i++){
				for (k=1;k<=NZ;k++){

					vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
					vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
					vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
					vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
					vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
					vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
					vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
					vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
					vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

					h1 = i-nx2+FW;

					psi_vxx[j][h1][k] = b_x[h1] * psi_vxx[j][h1][k] + a_x[h1] * vxx;
					vxx = vxx /K_x[h1] + psi_vxx[j][h1][k];
					psi_vyx[j][h1][k] = b_x_half[h1] * psi_vyx[j][h1][k] + a_x_half[h1] * vyx;
					vyx = vyx  /K_x_half[h1] + psi_vyx[j][h1][k];
					psi_vzx[j][h1][k] = b_x_half[h1] * psi_vzx[j][h1][k] + a_x_half[h1] * vzx;
					vzx = vzx / K_x_half[h1]  +psi_vzx[j][h1][k];

					if((POS[2]==0 && FREE_SURF==0) && (j<=FW)){

						psi_vxy[j][i][k] = b_y_half[j] * psi_vxy[j][i][k] + a_y_half[j] * vxy;
						vxy = vxy / K_y_half[j] + psi_vxy[j][i][k];
						psi_vyy[j][i][k] = b_y[j] * psi_vyy[j][i][k] + a_y[j] * vyy;
						vyy = vyy / K_y[j] + psi_vyy[j][i][k];
						psi_vzy[j][i][k] = b_y_half[j] * psi_vzy[j][i][k] + a_y_half[j] * vzy;
						vzy = vzy / K_y_half[j] + psi_vzy[j][i][k]; }


					if((POS[2]==NPROCY-1) && (j>=ny2+1)){
						h1 = (j-ny2+FW);

						psi_vxy[h1][i][k] = b_y_half[h1] * psi_vxy[h1][i][k] + a_y_half[h1] * vxy;
						vxy = vxy / K_y_half[h1] + psi_vxy[h1][i][k];
						psi_vyy[h1][i][k] = b_y[h1] * psi_vyy[h1][i][k] + a_y[h1] * vyy;
						vyy = vyy / K_y[h1] + psi_vyy[h1][i][k];
						psi_vzy[h1][i][k] = b_y_half[h1] * psi_vzy[h1][i][k] + a_y_half[h1] * vzy;
						vzy = vzy / K_y_half[h1] + psi_vzy[h1][i][k]; }


					if((POS[3]==0) && (k<=FW)){

						psi_vxz[j][i][k] = b_z_half[k] * psi_vxz[j][i][k] + a_z_half[k] * vxz;
						vxz = vxz / K_z_half[k] + psi_vxz[j][i][k];
						psi_vyz[j][i][k] = b_z_half[k] * psi_vyz[j][i][k] + a_z_half[k] * vyz;
						vyz = vyz / K_z_half[k] + psi_vyz[j][i][k];
						psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
						vzz = vzz / K_z[k] + psi_vzz[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){

						h1 = (k-nz2+FW);

						psi_vxz[j][i][h1] = b_z_half[h1] * psi_vxz[j][i][h1] + a_z_half[h1] * vxz;
						vxz = vxz / K_z_half[h1] + psi_vxz[j][i][h1];
						psi_vyz[j][i][h1] = b_z_half[h1] * psi_vyz[j][i][h1] + a_z_half[h1] * vyz;
						vyz = vyz / K_z_half[h1] + psi_vyz[j][i][h1];
						psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
						vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];}




					fipjp=uipjp[j][i][k]*DT;
					fjpkp=ujpkp[j][i][k]*DT;
					fipkp=uipkp[j][i][k]*DT;
					g=pi[j][i][k];
					f=2.0*u[j][i][k];

					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;

					sxy[j][i][k]+=(fipjp*vxyyx);
					syz[j][i][k]+=(fjpkp*vyzzy);
					sxz[j][i][k]+=(fipkp*vxzzx);
					sxx[j][i][k]+=DT*((g*vxxyyzz)-(f*vyyzz));
					syy[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxzz));
					szz[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxyy));


				}
			}
		}
	}

	if((POS[2]==0 && FREE_SURF==0)){
		for (j=1;j<=FW;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=1;k<=NZ;k++){

					vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
					vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
					vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
					vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
					vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
					vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
					vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
					vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
					vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

					psi_vxy[j][i][k] = b_y_half[j] * psi_vxy[j][i][k] + a_y_half[j] * vxy;
					vxy = vxy / K_y_half[j] + psi_vxy[j][i][k];
					psi_vyy[j][i][k] = b_y[j] * psi_vyy[j][i][k] + a_y[j] * vyy;
					vyy = vyy / K_y[j] + psi_vyy[j][i][k];
					psi_vzy[j][i][k] = b_y_half[j] * psi_vzy[j][i][k] + a_y_half[j] * vzy;
					vzy = vzy / K_y_half[j] + psi_vzy[j][i][k];

					if((POS[3]==0) && (k<=FW)){

						psi_vxz[j][i][k] = b_z_half[k] * psi_vxz[j][i][k] + a_z_half[k] * vxz;
						vxz = vxz / K_z_half[k] + psi_vxz[j][i][k];
						psi_vyz[j][i][k] = b_z_half[k] * psi_vyz[j][i][k] + a_z_half[k] * vyz;
						vyz = vyz / K_z_half[k] + psi_vyz[j][i][k];
						psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
						vzz = vzz / K_z[k] + psi_vzz[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){

						h1 = (k-nz2+FW);

						psi_vxz[j][i][h1] = b_z_half[h1] * psi_vxz[j][i][h1] + a_z_half[h1] * vxz;
						vxz = vxz / K_z_half[h1] + psi_vxz[j][i][h1];
						psi_vyz[j][i][h1] = b_z_half[h1] * psi_vyz[j][i][h1] + a_z_half[h1] * vyz;
						vyz = vyz / K_z_half[h1] + psi_vyz[j][i][h1];
						psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
						vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];}


					fipjp=uipjp[j][i][k]*DT;
					fjpkp=ujpkp[j][i][k]*DT;
					fipkp=uipkp[j][i][k]*DT;
					g=pi[j][i][k];
					f=2.0*u[j][i][k];

					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;

					sxy[j][i][k]+=(fipjp*vxyyx);
					syz[j][i][k]+=(fjpkp*vyzzy);
					sxz[j][i][k]+=(fipkp*vxzzx);
					sxx[j][i][k]+=DT*((g*vxxyyzz)-(f*vyyzz));
					syy[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxzz));
					szz[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxyy));


				}
			}
		}
	}

	if(POS[2]==NPROCY-1){
		for (j=ny2+1;j<=ny2+FW;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=1;k<=NZ;k++){

					vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
					vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
					vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
					vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
					vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
					vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
					vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
					vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
					vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

					h1 = (j-ny2+FW);

					psi_vxy[h1][i][k] = b_y_half[h1] * psi_vxy[h1][i][k] + a_y_half[h1] * vxy;
					vxy = vxy / K_y_half[h1] + psi_vxy[h1][i][k];
					psi_vyy[h1][i][k] = b_y[h1] * psi_vyy[h1][i][k] + a_y[h1] * vyy;
					vyy = vyy / K_y[h1] + psi_vyy[h1][i][k];
					psi_vzy[h1][i][k] = b_y_half[h1] * psi_vzy[h1][i][k] + a_y_half[h1] * vzy;
					vzy = vzy / K_y_half[h1] + psi_vzy[h1][i][k];

					if((POS[3]==0) && (k<=FW)){

						psi_vxz[j][i][k] = b_z_half[k] * psi_vxz[j][i][k] + a_z_half[k] * vxz;
						vxz = vxz / K_z_half[k] + psi_vxz[j][i][k];
						psi_vyz[j][i][k] = b_z_half[k] * psi_vyz[j][i][k] + a_z_half[k] * vyz;
						vyz = vyz / K_z_half[k] + psi_vyz[j][i][k];
						psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
						vzz = vzz / K_z[k] + psi_vzz[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){

						h1 = (k-nz2+FW);

						psi_vxz[j][i][h1] = b_z_half[h1] * psi_vxz[j][i][h1] + a_z_half[h1] * vxz;
						vxz = vxz / K_z_half[h1] + psi_vxz[j][i][h1];
						psi_vyz[j][i][h1] = b_z_half[h1] * psi_vyz[j][i][h1] + a_z_half[h1] * vyz;
						vyz = vyz / K_z_half[h1] + psi_vyz[j][i][h1];
						psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
						vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];}


					fipjp=uipjp[j][i][k]*DT;
					fjpkp=ujpkp[j][i][k]*DT;
					fipkp=uipkp[j][i][k]*DT;
					g=pi[j][i][k];
					f=2.0*u[j][i][k];

					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;

					sxy[j][i][k]+=(fipjp*vxyyx);
					syz[j][i][k]+=(fjpkp*vyzzy);
					sxz[j][i][k]+=(fipkp*vxzzx);
					sxx[j][i][k]+=DT*((g*vxxyyzz)-(f*vyyzz));
					syy[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxzz));
					szz[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxyy));


				}
			}
		}
	}


	if(POS[3]==0){
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=1;k<=FW;k++){

					vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
					vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
					vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
					vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
					vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
					vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
					vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
					vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
					vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

					psi_vxz[j][i][k] = b_z_half[k] * psi_vxz[j][i][k] + a_z_half[k] * vxz;
					vxz = vxz / K_z_half[k] + psi_vxz[j][i][k];
					psi_vyz[j][i][k] = b_z_half[k] * psi_vyz[j][i][k] + a_z_half[k] * vyz;
					vyz = vyz / K_z_half[k] + psi_vyz[j][i][k];
					psi_vzz[j][i][k] = b_z[k] * psi_vzz[j][i][k] + a_z[k] * vzz;
					vzz = vzz / K_z[k] + psi_vzz[j][i][k];



					fipjp=uipjp[j][i][k]*DT;
					fjpkp=ujpkp[j][i][k]*DT;
					fipkp=uipkp[j][i][k]*DT;
					g=pi[j][i][k];
					f=2.0*u[j][i][k];

					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;

					sxy[j][i][k]+=(fipjp*vxyyx);
					syz[j][i][k]+=(fjpkp*vyzzy);
					sxz[j][i][k]+=(fipkp*vxzzx);
					sxx[j][i][k]+=DT*((g*vxxyyzz)-(f*vyyzz));
					syy[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxzz));
					szz[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxyy));

				}
			}
		}
	}


	if(POS[3]==NPROCZ-1){		
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz2+1;k<=nz2+FW;k++){

					vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
					vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
					vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
					vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
					vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
					vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
					vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
					vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
					vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;

					h1 = (k-nz2+FW);

					psi_vxz[j][i][h1] = b_z_half[h1] * psi_vxz[j][i][h1] + a_z_half[h1] * vxz;
					vxz = vxz / K_z_half[h1] + psi_vxz[j][i][h1];
					psi_vyz[j][i][h1] = b_z_half[h1] * psi_vyz[j][i][h1] + a_z_half[h1] * vyz;
					vyz = vyz / K_z_half[h1] + psi_vyz[j][i][h1];
					psi_vzz[j][i][h1] = b_z[h1] * psi_vzz[j][i][h1] + a_z[h1] * vzz;
					vzz = vzz / K_z[h1] + psi_vzz[j][i][h1];





					fipjp=uipjp[j][i][k]*DT;
					fjpkp=ujpkp[j][i][k]*DT;
					fipkp=uipkp[j][i][k]*DT;
					g=pi[j][i][k];
					f=2.0*u[j][i][k];

					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;

					sxy[j][i][k]+=(fipjp*vxyyx);
					syz[j][i][k]+=(fjpkp*vyzzy);
					sxz[j][i][k]+=(fipkp*vxzzx);
					sxx[j][i][k]+=DT*((g*vxxyyzz)-(f*vyyzz));
					syy[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxzz));
					szz[j][i][k]+=DT*((g*vxxyyzz)-(f*vxxyy));

				}
			}
		}
	}


	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0){
			time2=MPI_Wtime();
			time=time2-time1;
			fprintf(FP," Real time for CPML stress tensor update: \t %4.2f s.\n",time);
		}
	return time;

}
