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
 *   stress free surface condition
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_acoustic(int ndepth,  float *** pi, float *** sxx, float *** vx, float *** vy, float *** vz){

	int i, k ,j, m, fdoh;
	float  vxx, vyy, vzz;
	float  h, g;
	

	extern int NX, NZ, FDCOEFF, FDORDER;
	extern float DT, DX, DY, DZ;
        register float b1, b2, b3, b4, b5, b6, dx, dy, dz;
	
        fdoh=FDORDER/2;

	j=ndepth;     /* The free surface is located exactly in y=(ndepth-1/2)*dh meter!! */

	
	switch (FDORDER){
	case 2 :
		dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;
		
	for (k=1;k<=NZ;k++){
		for (i=1;i<=NX;i++){


			/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
			sxx[j][i][k]=0.0;
			
			for (m=1; m<=fdoh; m++) {
			  sxx[j-m][i][k]=-sxx[j+m][i][k];
			}
			
                        vxx=dx*(vx[j][i][k]-vx[j][i-1][k]);
			vyy=dy*(vy[j][i][k]-vy[j-1][i][k]);
			vzz=dz*(vz[j][i][k]-vz[j][i][k-1]);


			/* partially updating sxx and szz in the same way*/
			g=pi[j][i][k];
			h=-(g*(vxx+vyy+vzz));
			sxx[j][i][k]+=h;


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
	
	for (k=1;k<=NZ;k++){
		for (i=1;i<=NX;i++){


			/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
			sxx[j][i][k]=0.0;
			
			for (m=1; m<=fdoh; m++) {
			  sxx[j-m][i][k]=-sxx[j+m][i][k];
			}
			
                                vxx=dx*(b1*(vx[j][i][k]-vx[j][i-1][k])+
				    b2*(vx[j][i+1][k]-vx[j][i-2][k]));
				vyy=dy*(b1*(vy[j][i][k]-vy[j-1][i][k])+
				    b2*(vy[j+1][i][k]-vy[j-2][i][k]));
				vzz=dz*(b1*(vz[j][i][k]-vz[j][i][k-1])+
				    b2*(vz[j][i][k+1]-vz[j][i][k-2]));


			/* partially updating sxx and szz in the same way*/
			g=pi[j][i][k];
			h=-(g*(vxx+vyy+vzz));
			sxx[j][i][k]+=h;

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
	
	for (k=1;k<=NZ;k++){
		for (i=1;i<=NX;i++){


			/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
			sxx[j][i][k]=0.0;
			
			for (m=1; m<=fdoh; m++) {
			  sxx[j-m][i][k]=-sxx[j+m][i][k];
			}
			

				vxx=dx*(b1*(vx[j][i][k]-vx[j][i-1][k])+
				    b2*(vx[j][i+1][k]-vx[j][i-2][k])+
				    b3*(vx[j][i+2][k]-vx[j][i-3][k]));
				vyy=dy*(b1*(vy[j][i][k]-vy[j-1][i][k])+
				    b2*(vy[j+1][i][k]-vy[j-2][i][k])+
				    b3*(vy[j+2][i][k]-vy[j-3][i][k]));
				vzz=dz*(b1*(vz[j][i][k]-vz[j][i][k-1])+
				    b2*(vz[j][i][k+1]-vz[j][i][k-2])+
				    b3*(vz[j][i][k+2]-vz[j][i][k-3]));


			/* partially updating sxx and szz in the same way*/
			g=pi[j][i][k];
			h=-(g*(vxx+vyy+vzz));
			sxx[j][i][k]+=h;


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
	
	for (k=1;k<=NZ;k++){
		for (i=1;i<=NX;i++){


			/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
			sxx[j][i][k]=0.0;
			
			for (m=1; m<=fdoh; m++) {
			  sxx[j-m][i][k]=-sxx[j+m][i][k];
			}
			
                                vxx=dx*(b1*(vx[j][i][k]-vx[j][i-1][k])+
				    b2*(vx[j][i+1][k]-vx[j][i-2][k])+
				    b3*(vx[j][i+2][k]-vx[j][i-3][k])+
				    b4*(vx[j][i+3][k]-vx[j][i-4][k]));
				vyy=dy*(b1*(vy[j][i][k]-vy[j-1][i][k])+
				    b2*(vy[j+1][i][k]-vy[j-2][i][k])+
				    b3*(vy[j+2][i][k]-vy[j-3][i][k])+
				    b4*(vy[j+3][i][k]-vy[j-4][i][k]));
				vzz=dz*(b1*(vz[j][i][k]-vz[j][i][k-1])+
				    b2*(vz[j][i][k+1]-vz[j][i][k-2])+
				    b3*(vz[j][i][k+2]-vz[j][i][k-3])+
				    b4*(vz[j][i][k+3]-vz[j][i][k-4]));


			/* partially updating sxx and szz in the same way*/
			g=pi[j][i][k];
			h=-(g*(vxx+vyy+vzz));
			sxx[j][i][k]+=h;


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
	
	for (k=1;k<=NZ;k++){
		for (i=1;i<=NX;i++){


			/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
			sxx[j][i][k]=0.0;
			
			for (m=1; m<=fdoh; m++) {
			  sxx[j-m][i][k]=-sxx[j+m][i][k];
			}
			
			        vxx=dx*(b1*(vx[j][i][k]-vx[j][i-1][k])+
				    b2*(vx[j][i+1][k]-vx[j][i-2][k])+
				    b3*(vx[j][i+2][k]-vx[j][i-3][k])+
				    b4*(vx[j][i+3][k]-vx[j][i-4][k])+
				    b5*(vx[j][i+4][k]-vx[j][i-5][k]));
				vyy=dy*(b1*(vy[j][i][k]-vy[j-1][i][k])+
				    b2*(vy[j+1][i][k]-vy[j-2][i][k])+
				    b3*(vy[j+2][i][k]-vy[j-3][i][k])+
				    b4*(vy[j+3][i][k]-vy[j-4][i][k])+
				    b5*(vy[j+4][i][k]-vy[j-5][i][k]));
				vzz=dz*(b1*(vz[j][i][k]-vz[j][i][k-1])+
				    b2*(vz[j][i][k+1]-vz[j][i][k-2])+
				    b3*(vz[j][i][k+2]-vz[j][i][k-3])+
				    b4*(vz[j][i][k+3]-vz[j][i][k-4])+
				    b5*(vz[j][i][k+4]-vz[j][i][k-5]));



			/* partially updating sxx and szz in the same way*/
			g=pi[j][i][k];
			h=-(g*(vxx+vyy+vzz));
			sxx[j][i][k]+=h;


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
	
	for (k=1;k<=NZ;k++){
		for (i=1;i<=NX;i++){


			/*Mirroring the components of the stress tensor to make
					  a stress free surface (method of imaging, Levander, 1988)*/
			sxx[j][i][k]=0.0;
			
			for (m=1; m<=fdoh; m++) {
			  sxx[j-m][i][k]=-sxx[j+m][i][k];
			}
			

				vxx=dx*(b1*(vx[j][i][k]-vx[j][i-1][k])+
				    b2*(vx[j][i+1][k]-vx[j][i-2][k])+
				    b3*(vx[j][i+2][k]-vx[j][i-3][k])+
				    b4*(vx[j][i+3][k]-vx[j][i-4][k])+
				    b5*(vx[j][i+4][k]-vx[j][i-5][k])+
				    b6*(vx[j][i+5][k]-vx[j][i-6][k]));
				vyy=dy*(b1*(vy[j][i][k]-vy[j-1][i][k])+
				    b2*(vy[j+1][i][k]-vy[j-2][i][k])+
				    b3*(vy[j+2][i][k]-vy[j-3][i][k])+
				    b4*(vy[j+3][i][k]-vy[j-4][i][k])+
				    b5*(vy[j+4][i][k]-vy[j-5][i][k])+
				    b6*(vy[j+5][i][k]-vy[j-6][i][k]));
				vzz=dz*(b1*(vz[j][i][k]-vz[j][i][k-1])+
				    b2*(vz[j][i][k+1]-vz[j][i][k-2])+
				    b3*(vz[j][i][k+2]-vz[j][i][k-3])+
				    b4*(vz[j][i][k+3]-vz[j][i][k-4])+
				    b5*(vz[j][i][k+4]-vz[j][i][k-5])+
				    b6*(vz[j][i][k+5]-vz[j][i][k-6]));


			/* partially updating sxx and szz in the same way*/
			g=pi[j][i][k];
			h=-(g*(vxx+vyy+vzz));
			sxx[j][i][k]+=h;


		}
	}
	
       break;	
       }
}
