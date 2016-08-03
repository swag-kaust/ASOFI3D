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
/* ----------------------------------------------------------------------
 * Initialization of the wave field with zero values (zero wavefield)
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"
void zero_elastic_CPML( int NX, int NY, int NZ, float *** vx, float *** vy, float *** vz,float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz, float *** rxx, float *** ryy, float *** rzz, float *** rxy, float *** ryz, float *** rxz, float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y, float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z, float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy, float *** psi_vxz, float *** psi_vyz, float *** psi_vzz,float *** rxx_2, float *** ryy_2,
                       float *** rzz_2, float *** rxy_2, float *** ryz_2, float *** rxz_2,float *** rxx_3, float *** ryy_3,
                       float *** rzz_3, float *** rxy_3, float *** ryz_3, float *** rxz_3,float *** rxx_4, float *** ryy_4,
                       float *** rzz_4, float *** rxy_4, float *** ryz_4, float *** rxz_4){


    extern int FDORDER, ABS_TYPE, FW, POS[4],L,FDORDER_TIME;
    int nx1, ny1, nz1, nx2, ny2, nz2, a,b,l, i, j, k;

    l=1;
    if(ABS_TYPE==1 && FDORDER==2){l=2;}
	
	if(POS[2]==0){
	  a=0;
	  b=1;}
	else{
	  a=1;
	  b=l;}
	

    ny1=a-b*FDORDER/2;
    ny2=NY+l*FDORDER/2;
    nx1=1-l*FDORDER/2;
    nx2=NX+l*FDORDER/2;
    nz1=1-l*FDORDER/2;
    nz2=NZ+l*FDORDER/2;

	  

    for (j=ny1;j<=ny2;j++){
      for (i=nx1;i<=nx2;i++){
	for (k=nz1;k<=nz2;k++){
	  vx[j][i][k]=0.0;
	  vy[j][i][k]=0.0;
	  vz[j][i][k]=0.0;
	  sxy[j][i][k]=0.0;
	  syz[j][i][k]=0.0;
	}
      }
    }	
	
    ny1=1-l*FDORDER/2;
   
    for (j=ny1;j<=ny2;j++){
      for (i=nx1;i<=nx2;i++){
	for (k=nz1;k<=nz2;k++){	  
	  sxx[j][i][k]=0.0;
	  sxz[j][i][k]=0.0;
	  syy[j][i][k]=0.0;
	  szz[j][i][k]=0.0;
	}
      }
    }		
	
	
	
	
    if(L){	
	for (j=1;j<=NY;j++){
	  for (i=1;i<=NX;i++){
	    for (k=1;k<=NZ;k++){
	      rxx[j][i][k]=0.0;
	      rxy[j][i][k]=0.0;
	      rxz[j][i][k]=0.0;
	      ryy[j][i][k]=0.0;
	      ryz[j][i][k]=0.0;
	      rzz[j][i][k]=0.0;
            if (FDORDER_TIME>2) {
                rxx_2[j][i][k]=0.0;
                ryy_2[j][i][k]=0.0;
                rzz_2[j][i][k]=0.0;
                rxy_2[j][i][k]=0.0;
                ryz_2[j][i][k]=0.0;
                rxz_2[j][i][k]=0.0;
                rxx_3[j][i][k]=0.0;
                ryy_3[j][i][k]=0.0;
                rzz_3[j][i][k]=0.0;
                rxy_3[j][i][k]=0.0;
                ryz_3[j][i][k]=0.0;
                rxz_3[j][i][k]=0.0;
                if (FDORDER_TIME==4) {
                    rxx_4[j][i][k]=0.0;
                    ryy_4[j][i][k]=0.0;
                    rzz_4[j][i][k]=0.0;
                    rxy_4[j][i][k]=0.0;
                    ryz_4[j][i][k]=0.0;
                    rxz_4[j][i][k]=0.0;
                }
            }
	    }
	  }
	}
    }
    
    for (j=1;j<=NY;j++){
      for (i=1;i<=NX;i++){
	for (k=1;k<=2*FW;k++){
	 psi_sxz_z[j][i][k]=0.0;
	 psi_syz_z[j][i][k]=0.0;
	 psi_szz_z[j][i][k]=0.0;
	 psi_vxz[j][i][k]=0.0;
	 psi_vyz[j][i][k]=0.0;
	 psi_vzz[j][i][k]=0.0;
	}
      }
    }

     for (j=1;j<=NY;j++){
      for (i=1;i<=2*FW;i++){
	for (k=1;k<=NZ;k++){
	 psi_sxx_x[j][i][k]=0.0;
	 psi_sxy_x[j][i][k]=0.0;
	 psi_sxz_x[j][i][k]=0.0;
	 psi_vxx[j][i][k]=0.0;
	 psi_vyx[j][i][k]=0.0;
	 psi_vzx[j][i][k]=0.0;
	}
      }
    }   

for (j=1;j<=2*FW;j++){
      for (i=1;i<=NX;i++){
	for (k=1;k<=NZ;k++){
	 psi_sxy_y[j][i][k]=0.0;
	 psi_syy_y[j][i][k]=0.0;
	 psi_syz_y[j][i][k]=0.0;
	 psi_vxy[j][i][k]=0.0;
	 psi_vyy[j][i][k]=0.0;
	 psi_vzy[j][i][k]=0.0;
	}
      }
    }   
}
