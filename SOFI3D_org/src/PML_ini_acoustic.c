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
 *   Definition of PML and non-PML domains 
 *   new version of the PML_ini program
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void PML_ini_acoustic(int * xa, int * xb, int * ya, int * yb, int * za, int * zb){


	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BLOCK, FW;
	
	int ifw, i, x0, y0, z0, x1, y1, z1;
	

	ifw = FW;
	i=1;
        x0=1; y0=1; z0=1;
        x1=NX; y1=NY; z1=NZ;


/* PML domain in x-direction*/

	if (POS[1]==0){
	
	xa[i] = x0; xb[i] = ifw; x0=ifw+1;
	ya[i] = y0; yb[i] = y1;
	za[i] = z0; zb[i] = z1; 
	
        i=i+1;

	}

        if (POS[1]==NPROCX-1){

        xa[i] = NX-ifw+1; xb[i] = x1; x1=NX-ifw;
	ya[i] = y0; yb[i] = y1;
	za[i] = z0; zb[i] = z1; 
	
        i=i+1;
        
        }


/* PML domain in y-direction*/

        if (POS[2]==0){
	
	xa[i] = x0; xb[i] = x1; 
	ya[i] = y0; yb[i] = ifw;  y0=ifw+1;
	za[i] = z0; zb[i] = z1; 
	
        i=i+1;

	}

        if (POS[2]==NPROCY-1){
	
	xa[i] = x0; xb[i] = x1; 
	ya[i] = NY-ifw+1; yb[i] = y1;  y1=NY-ifw;
	za[i] = z0; zb[i] = z1; 
	
        i=i+1;

	}


/* PML domain in z-direction*/

        if (POS[3]==0){
	
	xa[i] = x0; xb[i] = x1; 
	ya[i] = y0; yb[i] = y1; 
	za[i] = z0; zb[i] = ifw;  z0=ifw+1;
	
        i=i+1;

	}

        if (POS[3]==NPROCZ-1){
	
	xa[i] = x0; xb[i] = x1; 
	ya[i] = y0; yb[i] = y1;  
	za[i] = NZ-ifw+1; zb[i] = z1;  z1=NZ-ifw;
	
        i=i+1;

	}

       BLOCK=i-1;


/* calculation area */

	xa[0] = x0; xb[0] = x1;
	ya[0] = y0; yb[0] = y1;
	za[0] = z0; zb[0] = z1; 

}
	
