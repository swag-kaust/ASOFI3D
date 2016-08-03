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
 *  Initialization of the wave field with zero values (zero wavefield)
 *  
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_elastic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *** vx, float *** vy, float *** vz,
float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz,
                  float *** vxyyx,float *** vyzzy,float *** vxzzx,float *** vxxyyzz,float *** vyyzz,float *** vxxzz,float *** vxxyy,
                  float *** vxyyx_2,float *** vyzzy_2,float *** vxzzx_2,float *** vxxyyzz_2,float *** vyyzz_2,float *** vxxzz_2,float *** vxxyy_2,
                  float *** vxyyx_3,float *** vyzzy_3,float *** vxzzx_3,float *** vxxyyzz_3,float *** vyyzz_3,float *** vxxzz_3,float *** vxxyy_3,
                  float *** vxyyx_4,float *** vyzzy_4,float *** vxzzx_4,float *** vxxyyzz_4,float *** vyyzz_4,float *** vxxzz_4,float *** vxxyy_4,
                  float *** svx, float *** svy, float *** svz,
                  float *** svx_2, float *** svy_2, float *** svz_2, float *** svx_3, float *** svy_3, float *** svz_3,
                  float *** svx_4, float *** svy_4, float *** svz_4){


    extern int FDORDER_TIME;
	register int i, j, k;
    
    for (j=ny1;j<=ny2;j++){
        for (i=nx1;i<=nx2;i++){
            for (k=nz1;k<=nz2;k++){
                vx[j][i][k]=0.0;
                vy[j][i][k]=0.0;
                vz[j][i][k]=0.0;
                sxx[j][i][k]=0.0;
                syy[j][i][k]=0.0;
                szz[j][i][k]=0.0;
                sxy[j][i][k]=0.0;
                syz[j][i][k]=0.0;
                sxz[j][i][k]=0.0;
                
                if(FDORDER_TIME > 2){
                    vxyyx[j][i][k]=0.0;
                    vyzzy[j][i][k]=0.0;
                    vxzzx[j][i][k]=0.0;
                    vxxyyzz[j][i][k]=0.0;
                    vyyzz[j][i][k]=0.0;
                    vxxzz[j][i][k]=0.0;
                    vxxyy[j][i][k]=0.0;
                    
                    vxyyx_2[j][i][k]=0.0;
                    vyzzy_2[j][i][k]=0.0;
                    vxzzx_2[j][i][k]=0.0;
                    vxxyyzz_2[j][i][k]=0.0;
                    vyyzz_2[j][i][k]=0.0;
                    vxxzz_2[j][i][k]=0.0;
                    vxxyy_2[j][i][k]=0.0;
                    
                    vxyyx_3[j][i][k]=0.0;
                    vyzzy_3[j][i][k]=0.0;
                    vxzzx_3[j][i][k]=0.0;
                    vxxyyzz_3[j][i][k]=0.0;
                    vyyzz_3[j][i][k]=0.0;
                    vxxzz_3[j][i][k]=0.0;
                    vxxyy_3[j][i][k]=0.0;
                    
                    
                    svx[j][i][k]=0.0;
                    svy[j][i][k]=0.0;
                    svz[j][i][k]=0.0;
                    
                    svx_2[j][i][k]=0.0;
                    svy_2[j][i][k]=0.0;
                    svz_2[j][i][k]=0.0;
                    
                    svx_3[j][i][k]=0.0;
                    svy_3[j][i][k]=0.0;
                    svz_3[j][i][k]=0.0;
                    
                    if (FDORDER_TIME==4){
                        vxyyx_4[j][i][k]=0.0;
                        vyzzy_4[j][i][k]=0.0;
                        vxzzx_4[j][i][k]=0.0;
                        vxxyyzz_4[j][i][k]=0.0;
                        vyyzz_4[j][i][k]=0.0;
                        vxxzz_4[j][i][k]=0.0;
                        vxxyy_4[j][i][k]=0.0;
                        
                        svx_4[j][i][k]=0.0;
                        svy_4[j][i][k]=0.0;
                        svz_4[j][i][k]=0.0;
                    }
                }
            }
        }
    }
    
}
