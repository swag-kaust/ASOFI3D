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
 *   by a staggered grid finite difference scheme of 4th order accuracy in space
 *   and second order accuracy in time
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,
		float *** syz, float *** sxz, float  ***  rho,  float  *** rjp, float  *** rkp, float  *** rip,
		float **  srcpos_loc, float ** signals, int nsrc, float *** absorb_coeff, int * stype, float *** svx, float *** svy, float *** svz,
        float *** svx_2, float *** svy_2, float *** svz_2, float *** svx_3, float *** svy_3, float *** svz_3,
        float *** svx_4, float *** svy_4, float *** svz_4){


	extern float DT, DX, DY, DZ, SOURCE_ALPHA, SOURCE_BETA;
	double time=0.0, time1=0.0, time2=0.0;
	extern int MYID, FDORDER, FDORDER_TIME, LOG, ABS_TYPE, FDCOEFF; /* variable "FW" removed, not in use */
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;

	int i, j, k, l;
	float  amp, alpha_rad, beta_rad;
	float b1, b2, b3, b4, b5, b6, dx, dy, dz;
	float sxx_x, sxy_y, sxz_z, syy_y, sxy_x, syz_z;
	float szz_z, sxz_x, syz_y;
    float c1, c2, c3, c4; /* Coefficients for Adam Bashforth */

    float **svx_j,**svy_j,**svz_j;
    float **svx_j_2,**svy_j_2,**svz_j_2;
    float **svx_j_3,**svy_j_3,**svz_j_3;
    float **svx_j_4,**svy_j_4,**svz_j_4;
    
    float *svx_j_i,*svy_j_i,*svz_j_i;
    float *svx_j_i_2,*svy_j_i_2,*svz_j_i_2;
    float *svx_j_i_3,*svy_j_i_3,*svz_j_i_3;
    float *svx_j_i_4,*svy_j_i_4,*svz_j_i_4;
    
    
	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) time1=MPI_Wtime();

    
    
    switch (FDORDER_TIME){
            
        case 2:
            
            switch (FDORDER){
                    
                case 2 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=1.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.00100; } /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*b1*(sxx[j][i+1][k]-sxx[j][i][k]);
                                sxy_y = dy*b1*(sxy[j][i][k]-sxy[j-1][i][k]);
                                sxz_z = dz*b1*(sxz[j][i][k]-sxz[j][i][k-1]); /* backward operator */
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((sxx_x + sxy_y +sxz_z)/rip[j][i][k]);
                                
                                syy_y = dy*b1*(syy[j+1][i][k]-syy[j][i][k]);
                                sxy_x = dx*b1*(sxy[j][i][k]-sxy[j][i-1][k]);
                                syz_z = dz*b1*(syz[j][i][k]-syz[j][i][k-1]);
                                
                                
                                vy[j][i][k]+= ((syy_y + sxy_x + syz_z)/rjp[j][i][k]);
                                
                                szz_z = dz*b1*(szz[j][i][k+1]-szz[j][i][k]);
                                sxz_x = dx*b1*(sxz[j][i][k]-sxz[j][i-1][k]);
                                syz_y = dy*b1*(syz[j][i][k]-syz[j-1][i][k]);
                                
                                
                                vz[j][i][k]+= ((szz_z + sxz_x + syz_y)/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 4 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
                                
                                
                                vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
                                
                                
                                vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 6 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3]));
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3]));
                                
                                
                                vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k]));
                                
                                
                                vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 8 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3])+
                                            b4*(sxz[j][i][k+3]-sxz[j][i][k-4]));
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4]));
                                
                                
                                vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k]));
                                
                                
                                vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 10 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    
                    b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor Coefficients*/
                    if(FDCOEFF==2){
                        b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
                                            b5*(sxx[j][i+5][k]-sxx[j][i-4][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k])+
                                            b5*(sxy[j+4][i][k]-sxy[j-5][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3])+
                                            b4*(sxz[j][i][k+3]-sxz[j][i][k-4])+
                                            b5*(sxz[j][i][k+4]-sxz[j][i][k-5]));
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k])+
                                            b5*(syy[j+5][i][k]-syy[j-4][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k])+
                                            b5*(sxy[j][i+4][k]-sxy[j][i-5][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4])+
                                            b5*(syz[j][i][k+4]-syz[j][i][k-5]));
                                
                                
                                vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3])+
                                            b5*(szz[j][i][k+5]-szz[j][i][k-4]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k])+
                                            b5*(sxz[j][i+4][k]-sxz[j][i-5][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k])+
                                            b5*(syz[j+4][i][k]-syz[j-5][i][k]));
                                
                                
                                vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];
                                
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
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
                                            b5*(sxx[j][i+5][k]-sxx[j][i-4][k])+
                                            b6*(sxx[j][i+6][k]-sxx[j][i-5][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k])+
                                            b5*(sxy[j+4][i][k]-sxy[j-5][i][k])+
                                            b6*(sxy[j+5][i][k]-sxy[j-6][i][k]));
                                
                                sxz_z = dz*(b1*(sxy[j][i][k]-sxy[j][i][k-1])+
                                            b2*(sxy[j][i][k+1]-sxy[j][i][k-2])+
                                            b3*(sxy[j][i][k+2]-sxy[j][i][k-3])+
                                            b4*(sxy[j][i][k+3]-sxy[j][i][k-4])+
                                            b5*(sxy[j][i][k+4]-sxy[j][i][k-5])+
                                            b6*(sxy[j][i][k+5]-sxy[j][i][k-6]));
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k])+
                                            b5*(syy[j+5][i][k]-syy[j-4][i][k])+
                                            b6*(syy[j+6][i][k]-syy[j-5][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k])+
                                            b5*(sxy[j][i+4][k]-sxy[j][i-5][k])+
                                            b6*(sxy[j][i+5][k]-sxy[j][i-6][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4])+
                                            b5*(syz[j][i][k+4]-syz[j][i][k-5])+
                                            b6*(syz[j][i][k+5]-syz[j][i][k-6]));
                                
                                
                                
                                vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3])+
                                            b5*(szz[j][i][k+5]-szz[j][i][k-4])+
                                            b6*(szz[j][i][k+6]-szz[j][i][k-5]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k])+
                                            b5*(sxz[j][i+4][k]-sxz[j][i-5][k])+
                                            b6*(sxz[j][i+5][k]-sxz[j][i-6][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k])+
                                            b5*(syz[j+4][i][k]-syz[j-5][i][k])+
                                            b6*(syz[j+5][i][k]-syz[j-6][i][k]));
                                
                                
                                vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];
                                
                            }
                        }
                    }
                    
                    break;
                    
            }
            break; /* break for FDORDER_TIME=2 */
            
        case 3:
            
            c1=25.0/24.0; c2=-1.0/12.0; c3=1.0/24.0;
            
            switch (FDORDER){
                    
                case 2 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=1.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.00100; } /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*b1*(sxx[j][i+1][k]-sxx[j][i][k]);
                                sxy_y = dy*b1*(sxy[j][i][k]-sxy[j-1][i][k]);
                                sxz_z = dz*b1*(sxz[j][i][k]-sxz[j][i][k-1]); /* backward operator */
                                
                                syy_y = dy*b1*(syy[j+1][i][k]-syy[j][i][k]);
                                sxy_x = dx*b1*(sxy[j][i][k]-sxy[j][i-1][k]);
                                syz_z = dz*b1*(syz[j][i][k]-syz[j][i][k-1]);
                                
                                szz_z = dz*b1*(szz[j][i][k+1]-szz[j][i][k]);
                                sxz_x = dx*b1*(sxz[j][i][k]-sxz[j][i-1][k]);
                                syz_y = dy*b1*(syz[j][i][k]-syz[j-1][i][k]);
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 4 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k)))/rkp[j][i][k]);
                                
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 6 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k]));
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 8 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3])+
                                            b4*(sxz[j][i][k+3]-sxz[j][i][k-4]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k]));
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 10 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    
                    b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor Coefficients*/
                    if(FDCOEFF==2){
                        b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
                                            b5*(sxx[j][i+5][k]-sxx[j][i-4][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k])+
                                            b5*(sxy[j+4][i][k]-sxy[j-5][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3])+
                                            b4*(sxz[j][i][k+3]-sxz[j][i][k-4])+
                                            b5*(sxz[j][i][k+4]-sxz[j][i][k-5]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k])+
                                            b5*(syy[j+5][i][k]-syy[j-4][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k])+
                                            b5*(sxy[j][i+4][k]-sxy[j][i-5][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4])+
                                            b5*(syz[j][i][k+4]-syz[j][i][k-5]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3])+
                                            b5*(szz[j][i][k+5]-szz[j][i][k-4]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k])+
                                            b5*(sxz[j][i+4][k]-sxz[j][i-5][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k])+
                                            b5*(syz[j+4][i][k]-syz[j-5][i][k]));
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k)))/rkp[j][i][k]);
                                
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
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
                                            b5*(sxx[j][i+5][k]-sxx[j][i-4][k])+
                                            b6*(sxx[j][i+6][k]-sxx[j][i-5][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k])+
                                            b5*(sxy[j+4][i][k]-sxy[j-5][i][k])+
                                            b6*(sxy[j+5][i][k]-sxy[j-6][i][k]));
                                
                                sxz_z = dz*(b1*(sxy[j][i][k]-sxy[j][i][k-1])+
                                            b2*(sxy[j][i][k+1]-sxy[j][i][k-2])+
                                            b3*(sxy[j][i][k+2]-sxy[j][i][k-3])+
                                            b4*(sxy[j][i][k+3]-sxy[j][i][k-4])+
                                            b5*(sxy[j][i][k+4]-sxy[j][i][k-5])+
                                            b6*(sxy[j][i][k+5]-sxy[j][i][k-6]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k])+
                                            b5*(syy[j+5][i][k]-syy[j-4][i][k])+
                                            b6*(syy[j+6][i][k]-syy[j-5][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k])+
                                            b5*(sxy[j][i+4][k]-sxy[j][i-5][k])+
                                            b6*(sxy[j][i+5][k]-sxy[j][i-6][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4])+
                                            b5*(syz[j][i][k+4]-syz[j][i][k-5])+
                                            b6*(syz[j][i][k+5]-syz[j][i][k-6]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3])+
                                            b5*(szz[j][i][k+5]-szz[j][i][k-4])+
                                            b6*(szz[j][i][k+6]-szz[j][i][k-5]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k])+
                                            b5*(sxz[j][i+4][k]-sxz[j][i-5][k])+
                                            b6*(sxz[j][i+5][k]-sxz[j][i-6][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k])+
                                            b5*(syz[j+4][i][k]-syz[j-5][i][k])+
                                            b6*(syz[j+5][i][k]-syz[j-6][i][k]));
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
            }
            break; /* break for FDORDER_TIME=3 */
            
            
        case 4:
            
            c1=13.0/12.0; c2=-5.0/24.0; c3=1.0/6.0; c4=-1.0/24.0;
            
            switch (FDORDER){
                    
                case 2 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=1.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.00100; } /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        svx_j_4=*(svx_4+j);svy_j_4=*(svy_4+j);svz_j_4=*(svz_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            svx_j_i_4=*(svx_j_4+i);svy_j_i_4=*(svy_j_4+i);svz_j_i_4=*(svz_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*b1*(sxx[j][i+1][k]-sxx[j][i][k]);
                                sxy_y = dy*b1*(sxy[j][i][k]-sxy[j-1][i][k]);
                                sxz_z = dz*b1*(sxz[j][i][k]-sxz[j][i][k-1]); /* backward operator */
                                
                                syy_y = dy*b1*(syy[j+1][i][k]-syy[j][i][k]);
                                sxy_x = dx*b1*(sxy[j][i][k]-sxy[j][i-1][k]);
                                syz_z = dz*b1*(syz[j][i][k]-syz[j][i][k-1]);
                                
                                szz_z = dz*b1*(szz[j][i][k+1]-szz[j][i][k]);
                                sxz_x = dx*b1*(sxz[j][i][k]-sxz[j][i-1][k]);
                                syz_y = dy*b1*(syz[j][i][k]-syz[j-1][i][k]);
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k))+c4*(*(svx_j_i_4+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k))+c4*(*(svy_j_i_4+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k))+c4*(*(svz_j_i_4+k)))/rkp[j][i][k]);
                            
                            }
                        }
                    }
                    
                    break;
                    
                case 4 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        svx_j_4=*(svx_4+j);svy_j_4=*(svy_4+j);svz_j_4=*(svz_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            svx_j_i_4=*(svx_j_4+i);svy_j_i_4=*(svy_j_4+i);svz_j_i_4=*(svz_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k))+c4*(*(svx_j_i_4+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k))+c4*(*(svy_j_i_4+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k))+c4*(*(svz_j_i_4+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 6 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        svx_j_4=*(svx_4+j);svy_j_4=*(svy_4+j);svz_j_4=*(svz_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            svx_j_i_4=*(svx_j_4+i);svy_j_i_4=*(svy_j_4+i);svz_j_i_4=*(svz_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k]));
                                

                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k))+c4*(*(svx_j_i_4+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k))+c4*(*(svy_j_i_4+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k))+c4*(*(svz_j_i_4+k)))/rkp[j][i][k]);
                            }
                        }
                    }
                    
                    break;
                    
                case 8 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients*/
                    if(FDCOEFF==2){
                        b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        svx_j_4=*(svx_4+j);svy_j_4=*(svy_4+j);svz_j_4=*(svz_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            svx_j_i_4=*(svx_j_4+i);svy_j_i_4=*(svy_j_4+i);svz_j_i_4=*(svz_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3])+
                                            b4*(sxz[j][i][k+3]-sxz[j][i][k-4]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k]));
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k))+c4*(*(svx_j_i_4+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k))+c4*(*(svy_j_i_4+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k))+c4*(*(svz_j_i_4+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
                case 10 :
                    
                    dx=DT/DX;
                    dy=DT/DY;
                    dz=DT/DZ;
                    
                    
                    b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor Coefficients*/
                    if(FDCOEFF==2){
                        b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        svx_j_4=*(svx_4+j);svy_j_4=*(svy_4+j);svz_j_4=*(svz_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            svx_j_i_4=*(svx_j_4+i);svy_j_i_4=*(svy_j_4+i);svz_j_i_4=*(svz_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
                                            b5*(sxx[j][i+5][k]-sxx[j][i-4][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k])+
                                            b5*(sxy[j+4][i][k]-sxy[j-5][i][k]));
                                
                                sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+
                                            b2*(sxz[j][i][k+1]-sxz[j][i][k-2])+
                                            b3*(sxz[j][i][k+2]-sxz[j][i][k-3])+
                                            b4*(sxz[j][i][k+3]-sxz[j][i][k-4])+
                                            b5*(sxz[j][i][k+4]-sxz[j][i][k-5]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k])+
                                            b5*(syy[j+5][i][k]-syy[j-4][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k])+
                                            b5*(sxy[j][i+4][k]-sxy[j][i-5][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4])+
                                            b5*(syz[j][i][k+4]-syz[j][i][k-5]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3])+
                                            b5*(szz[j][i][k+5]-szz[j][i][k-4]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k])+
                                            b5*(sxz[j][i+4][k]-sxz[j][i-5][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k])+
                                            b5*(syz[j+4][i][k]-syz[j-5][i][k]));
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k))+c4*(*(svx_j_i_4+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k))+c4*(*(svy_j_i_4+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k))+c4*(*(svz_j_i_4+k)))/rkp[j][i][k]);
                                
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
                        svx_j=*(svx+j);svy_j=*(svy+j);svz_j=*(svz+j);
                        svx_j_2=*(svx_2+j);svy_j_2=*(svy_2+j);svz_j_2=*(svz_2+j);
                        svx_j_3=*(svx_3+j);svy_j_3=*(svy_3+j);svz_j_3=*(svz_3+j);
                        svx_j_4=*(svx_4+j);svy_j_4=*(svy_4+j);svz_j_4=*(svz_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            svx_j_i=*(svx_j+i);svy_j_i=*(svy_j+i);svz_j_i=*(svz_j+i);
                            svx_j_i_2=*(svx_j_2+i);svy_j_i_2=*(svy_j_2+i);svz_j_i_2=*(svz_j_2+i);
                            svx_j_i_3=*(svx_j_3+i);svy_j_i_3=*(svy_j_3+i);svz_j_i_3=*(svz_j_3+i);
                            svx_j_i_4=*(svx_j_4+i);svy_j_i_4=*(svy_j_4+i);svz_j_i_4=*(svz_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+
                                            b2*(sxx[j][i+2][k]-sxx[j][i-1][k])+
                                            b3*(sxx[j][i+3][k]-sxx[j][i-2][k])+
                                            b4*(sxx[j][i+4][k]-sxx[j][i-3][k])+
                                            b5*(sxx[j][i+5][k]-sxx[j][i-4][k])+
                                            b6*(sxx[j][i+6][k]-sxx[j][i-5][k]));
                                
                                sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+
                                            b2*(sxy[j+1][i][k]-sxy[j-2][i][k])+
                                            b3*(sxy[j+2][i][k]-sxy[j-3][i][k])+
                                            b4*(sxy[j+3][i][k]-sxy[j-4][i][k])+
                                            b5*(sxy[j+4][i][k]-sxy[j-5][i][k])+
                                            b6*(sxy[j+5][i][k]-sxy[j-6][i][k]));
                                
                                sxz_z = dz*(b1*(sxy[j][i][k]-sxy[j][i][k-1])+
                                            b2*(sxy[j][i][k+1]-sxy[j][i][k-2])+
                                            b3*(sxy[j][i][k+2]-sxy[j][i][k-3])+
                                            b4*(sxy[j][i][k+3]-sxy[j][i][k-4])+
                                            b5*(sxy[j][i][k+4]-sxy[j][i][k-5])+
                                            b6*(sxy[j][i][k+5]-sxy[j][i][k-6]));
                                
                                syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+
                                            b2*(syy[j+2][i][k]-syy[j-1][i][k])+
                                            b3*(syy[j+3][i][k]-syy[j-2][i][k])+
                                            b4*(syy[j+4][i][k]-syy[j-3][i][k])+
                                            b5*(syy[j+5][i][k]-syy[j-4][i][k])+
                                            b6*(syy[j+6][i][k]-syy[j-5][i][k]));
                                
                                sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+
                                            b2*(sxy[j][i+1][k]-sxy[j][i-2][k])+
                                            b3*(sxy[j][i+2][k]-sxy[j][i-3][k])+
                                            b4*(sxy[j][i+3][k]-sxy[j][i-4][k])+
                                            b5*(sxy[j][i+4][k]-sxy[j][i-5][k])+
                                            b6*(sxy[j][i+5][k]-sxy[j][i-6][k]));
                                
                                syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+
                                            b2*(syz[j][i][k+1]-syz[j][i][k-2])+
                                            b3*(syz[j][i][k+2]-syz[j][i][k-3])+
                                            b4*(syz[j][i][k+3]-syz[j][i][k-4])+
                                            b5*(syz[j][i][k+4]-syz[j][i][k-5])+
                                            b6*(syz[j][i][k+5]-syz[j][i][k-6]));
                                
                                szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+
                                            b2*(szz[j][i][k+2]-szz[j][i][k-1])+
                                            b3*(szz[j][i][k+3]-szz[j][i][k-2])+
                                            b4*(szz[j][i][k+4]-szz[j][i][k-3])+
                                            b5*(szz[j][i][k+5]-szz[j][i][k-4])+
                                            b6*(szz[j][i][k+6]-szz[j][i][k-5]));
                                
                                sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+
                                            b2*(sxz[j][i+1][k]-sxz[j][i-2][k])+
                                            b3*(sxz[j][i+2][k]-sxz[j][i-3][k])+
                                            b4*(sxz[j][i+3][k]-sxz[j][i-4][k])+
                                            b5*(sxz[j][i+4][k]-sxz[j][i-5][k])+
                                            b6*(sxz[j][i+5][k]-sxz[j][i-6][k]));
                                
                                
                                syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+
                                            b2*(syz[j+1][i][k]-syz[j-2][i][k])+
                                            b3*(syz[j+2][i][k]-syz[j-3][i][k])+
                                            b4*(syz[j+3][i][k]-syz[j-4][i][k])+
                                            b5*(syz[j+4][i][k]-syz[j-5][i][k])+
                                            b6*(syz[j+5][i][k]-syz[j-6][i][k]));
                                
                                
                                
                                *(svx_j_i+k)=sxx_x + sxy_y +sxz_z;
                                *(svy_j_i+k)=syy_y + sxy_x + syz_z;
                                *(svz_j_i+k)=szz_z + sxz_x + syz_y;
                                
                                
                                /* updating components of particle velocities */
                                vx[j][i][k]+= ((c1*(*(svx_j_i+k))+c2*(*(svx_j_i_2+k))+c3*(*(svx_j_i_3+k))+c4*(*(svx_j_i_4+k)))/rip[j][i][k]);
                                vy[j][i][k]+= ((c1*(*(svy_j_i+k))+c2*(*(svy_j_i_2+k))+c3*(*(svy_j_i_3+k))+c4*(*(svy_j_i_4+k)))/rjp[j][i][k]);
                                vz[j][i][k]+= ((c1*(*(svz_j_i+k))+c2*(*(svz_j_i_2+k))+c3*(*(svz_j_i_3+k))+c4*(*(svz_j_i_4+k)))/rkp[j][i][k]);
                                
                            }
                        }
                    }
                    
                    break;
                    
            }
            break; /* break for FDORDER_TIME=4 */
            
    }
    
    /* Adding body force components to corresponding particle velocities */
    for (l=1;l<=nsrc;l++) {
        i=(int)srcpos_loc[1][l];
        j=(int)srcpos_loc[2][l];
        k=(int)srcpos_loc[3][l];
        //amp=signals[l][nt]; // unscaled force amplitude
        amp=(DT*signals[l][nt])/(DX*DY*DZ);// scaled force amplitude with F= 1N
        
        //fprintf(FP," amp at timestep nt %i = %5.5e with DX=%5.2f DY=%5.2f DZ=%5.2f  DT=%5.8f\n",nt,amp,DX,DY,DZ,DT);
        
        switch (stype[l]){
            case 2 :
                vx[j][i][k] += amp/rip[j][i][k];  /* single force in x */
                break;
            case 3 :
                vy[j][i][k] += amp/rjp[j][i][k];  /* single force in y / vertical direction */
                break;
            case 4 :
                vz[j][i][k] += amp/rkp[j][i][k];  /* single force in z */
                break;
            case 5 :
                alpha_rad=SOURCE_ALPHA*PI/180; /* custom force */
                beta_rad=SOURCE_BETA*PI/180;
                vx[j][i][k]+=cos(alpha_rad)*sin(beta_rad)*amp/rip[j][i][k];
                vy[j][i][k]+=cos(beta_rad)*amp/rjp[j][i][k]; /*vertical component*/
                vz[j][i][k]+=sin(alpha_rad)*sin(beta_rad)*amp/rkp[j][i][k];
                break;
            default:
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
                    sxy[j][i][k]*=absorb_coeff[j][i][k];
                    syz[j][i][k]*=absorb_coeff[j][i][k];
                    sxz[j][i][k]*=absorb_coeff[j][i][k];
                    sxx[j][i][k]*=absorb_coeff[j][i][k];
                    syy[j][i][k]*=absorb_coeff[j][i][k];
                    szz[j][i][k]*=absorb_coeff[j][i][k];
                    
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
