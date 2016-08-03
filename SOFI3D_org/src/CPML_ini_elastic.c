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
 *  Definition of CPML boundary domains
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void CPML_ini_elastic(int * xb, int * yb, int * zb){


	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, FW, FREE_SURF;

        xb[0]=1; 
	xb[1]=NX;
	yb[0]=1; 
	yb[1]=NY;
	zb[0]=1;
        zb[1]=NZ;


/* PML domain in x-direction*/

	if (POS[1]==0){
	xb[0]=FW+1;
	}

        if (POS[1]==NPROCX-1){
        xb[1]=NX-FW;
	}


/* PML domain in y-direction*/

        if (POS[2]==0 && FREE_SURF==0){
	yb[0]=FW+1;
	}

        if (POS[2]==NPROCY-1){
	yb[1]=NY-FW;
	}


/* PML domain in z-direction*/

        if (POS[3]==0){
	zb[0]=FW+1;
	}

        if (POS[3]==NPROCZ-1){
	zb[1]=NZ-FW;
	}

}

	
