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
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void PML_ini(int * xa, int * xb, int * ya, int * yb, int * za, int * zb){


	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BLOCK, FW;
	
	int ifw;
	
/* define PML and non-PML domains */
	
	/* corners */
	/* -----------------------------------------------------------------------------*/
        /* first corner */
	ifw = FW;
	
	if ((POS[1]==0)&&(POS[2]==0)&&(POS[3]==0)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = ifw+1; xb[2] = NX;
	ya[2] = 1; yb[2] = NY;
	za[2] = 1; zb[2] = ifw; 
	
        /* third block */
	xa[3] = ifw+1; xb[3] = NX;
	ya[3] = 1; yb[3] = ifw;
	za[3] = ifw+1; zb[3] = NZ; 
	
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}

	/* second corner */
	if ((POS[1]==0)&&(POS[2]==0)&&(POS[3]==NPROCZ-1)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = ifw+1; xb[2] = NX;
	ya[2] = 1; yb[2] = NY;
	za[2] = NZ - ifw; zb[2] = NZ; 
	
        /* third block */
	xa[3] = ifw+1; xb[3] = NX;
	ya[3] = 1; yb[3] = ifw;
	za[3] = 1; zb[3] = NZ - ifw - 1; 
	
	/* calculation area */
	xa[4] = ifw+1; xb[4] = NX;
	ya[4] = ifw+1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* third corner */
	if ((POS[3]==NPROCZ-1)&&(POS[1]==NPROCX-1)&&(POS[2]==0)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = 1; xb[2] = NX - ifw - 1;
	ya[2] = 1; yb[2] = NY;
	za[2] = NZ - ifw; zb[2] = NZ; 
	
        /* third block */
	xa[3] = 1; xb[3] = NX - ifw - 1;
	ya[3] = 1; yb[3] = ifw;
	za[3] = 1; zb[3] = NZ - ifw - 1; 
	
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ - ifw -1; 
	
	}
	
	/* fourth corner */
	if ((POS[3]==0)&&(POS[1]==NPROCX-1)&&(POS[2]==0)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = 1; xb[2] = NX - ifw - 1;
	ya[2] = 1; yb[2] = NY;
	za[2] = 1; zb[2] = ifw; 
	
        /* third block */
	xa[3] = 1; xb[3] = NX - ifw - 1;
	ya[3] = 1; yb[3] = ifw;
	za[3] = ifw + 1; zb[3] = NZ; 
	
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}

        /* fifth corner */
	if ((POS[1]==0)&&(POS[3]==0)&&(POS[2]==NPROCY-1)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = ifw+1; xb[2] = NX;
	ya[2] = 1; yb[2] = NY;
	za[2] = 1; zb[2] = ifw; 
	
        /* third block */
	xa[3] = ifw+1; xb[3] = NX;
	ya[3] = NY - ifw; yb[3] = NY;
	za[3] = ifw+1; zb[3] = NZ; 
	
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = ifw - 1; zb[4] = NZ; 
	
	}
	
	/* sixth corner */
	if ((POS[1]==0)&&(POS[2]==NPROCY-1)&&(POS[3]==NPROCZ-1)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = ifw+1; xb[2] = NX;
	ya[2] = 1; yb[2] = NY;
	za[2] = NZ - ifw; zb[2] = NZ; 
	
        /* third block */
	xa[3] = ifw+1; xb[3] = NX;
	ya[3] = NY - ifw; yb[3] = NY;
	za[3] = 1; zb[3] = NZ - ifw - 1; 
	
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* seventh corner */
	if ((POS[3]==NPROCZ-1)&&(POS[1]==NPROCX-1)&&(POS[2]==NPROCY-1)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = 1; xb[2] = NX - ifw - 1;
	ya[2] = 1; yb[2] = NY;
	za[2] = NZ - ifw; zb[2] = NZ; 
	
        /* third block */
	xa[3] = 1; xb[3] = NX - ifw - 1;
	ya[3] = NY - ifw; yb[3] = NY;
	za[3] = 1; zb[3] = NZ - ifw - 1; 
	
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* eighth corner */
	if ((POS[3]==0)&&(POS[1]==NPROCX-1)&&(POS[2]==NPROCY-1)){
	
	BLOCK = 3;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
	
	/* second block */
	xa[2] = 1; xb[2] = NX - ifw - 1;
	ya[2] = 1; yb[2] = NY;
	za[2] = 1; zb[2] = ifw; 
	
        /* third block */
	xa[3] = 1; xb[3] = NX - ifw - 1;
	ya[3] = NY - ifw; yb[3] = NY;
	za[3] = ifw + 1; zb[3] = NZ; 
	
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}
	
	/* faces */
	/* ----------------------------------------------------------------------------- */
	/* 1st face */
	
	if ((POS[1]==0)&&(POS[2]>0)&&(POS[2]<NPROCY-1)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 1;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
		
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
	/* 2nd face */
	if ((POS[1]==NPROCX-1)&&(POS[2]>0)&&(POS[2]<NPROCY-1)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 1;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
	/* 3rd face */
	if ((POS[3]==NPROCZ-1)&&(POS[2]>0)&&(POS[2]<NPROCY-1)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 1;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = NZ - ifw; zb[1] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* 4th face */
	if ((POS[3]==0)&&(POS[2]>0)&&(POS[2]<NPROCY-1)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 1;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = ifw; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}
	
	/* 5th face */
	if ((POS[2]==0)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 1;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = ifw;
	za[1] = 1; zb[1] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
        /* 6th face */
	if ((POS[2]==NPROCY-1)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 1;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = NY - ifw; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
	/* Edges */
	/* ------------------------------------------------------------------ */
	/* 1st edge */
	if ((POS[2]==NPROCY-1)&&(POS[3]==0)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = ifw; 

        /* second block */
	xa[2] = 1; xb[2] = NX;
	ya[2] = NY - ifw; yb[2] = NY;
	za[2] = ifw + 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}
	
	/* 2nd edge */
	if ((POS[2]==0)&&(POS[3]==0)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = ifw; 

        /* second block */
	xa[2] = 1; xb[2] = NX;
	ya[2] = 1; yb[2] = ifw;
	za[2] = ifw + 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}
	
	/* 3rd edge */
	if ((POS[1]==0)&&(POS[3]==0)&&(POS[2]>0)&&(POS[2]<NPROCY-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = ifw; 

        /* second block */
	xa[2] = 1; xb[2] = ifw;
	ya[2] = 1; yb[2] = NY;
	za[2] = ifw + 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}
	
	/* 4th edge */
	if ((POS[1]==NPROCX-1)&&(POS[3]==0)&&(POS[2]>0)&&(POS[2]<NPROCY-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = ifw; 

        /* second block */
	xa[2] = NX - ifw; xb[2] = NX;
	ya[2] = 1; yb[2] = NY;
	za[2] = ifw + 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = 1; yb[4] = NY;
	za[4] = ifw + 1; zb[4] = NZ; 
	
	}
	
	/* 5th edge */
	if ((POS[2]==NPROCY-1)&&(POS[3]==NPROCZ-1)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = NZ - ifw; zb[1] = NZ; 

        /* second block */
	xa[2] = 1; xb[2] = NX;
	ya[2] = NY - ifw; yb[2] = NY;
	za[2] = 1; zb[2] = NZ - ifw - 1; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* 6th edge */
	if ((POS[2]==0)&&(POS[3]==NPROCZ-1)&&(POS[1]>0)&&(POS[1]<NPROCX-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = NZ - ifw; zb[1] = NZ; 

        /* second block */
	xa[2] = 1; xb[2] = NX;
	ya[2] = 1; yb[2] = ifw;
	za[2] = 1; zb[2] = NZ - ifw - 1; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = ifw  + 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* 7th edge */
	if ((POS[1]==NPROCX-1)&&(POS[3]==NPROCZ-1)&&(POS[2]>0)&&(POS[2]<NPROCY-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = NZ - ifw; zb[1] = NZ; 

        /* second block */
	xa[2] = NX - ifw; xb[2] = NX;
	ya[2] = 1; yb[2] = NY;
	za[2] = 1; zb[2] = NZ - ifw - 1; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* 8th edge */
	if ((POS[1]==0)&&(POS[3]==NPROCZ-1)&&(POS[2]>0)&&(POS[2]<NPROCY-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = NZ - ifw; zb[1] = NZ; 

        /* second block */
	xa[2] = 1; xb[2] = ifw;
	ya[2] = 1; yb[2] = NY;
	za[2] = 1; zb[2] = NZ - ifw - 1; 
		
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ - ifw - 1; 
	
	}
	
	/* 9th edge */
	if ((POS[1]==0)&&(POS[2]==0)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 

        /* second block */
	xa[2] = ifw + 1; xb[2] = NX;
	ya[2] = 1; yb[2] = ifw;
	za[2] = 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
	/* 10th edge */
	if ((POS[1]==NPROCX-1)&&(POS[2]==0)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 

        /* second block */
	xa[2] = 1; xb[2] = NX - ifw - 1;
	ya[2] = 1; yb[2] = ifw;
	za[2] = 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = ifw + 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
	/* 11th edge */
	if ((POS[1]==NPROCX-1)&&(POS[2]==NPROCY-1)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = NX - ifw; xb[1] = NX;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 

        /* second block */
	xa[2] = 1; xb[2] = NX - ifw - 1;
	ya[2] = NY - ifw; yb[2] = NY;
	za[2] = 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = 1; xb[4] = NX - ifw - 1;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = 1; zb[4] = NZ; 
	
	}
	
	/* 12th edge */
	if ((POS[1]==0)&&(POS[2]==NPROCY-1)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 2;
	
	/* first block */
	xa[1] = 1; xb[1] = ifw;
	ya[1] = 1; yb[1] = NY;
	za[1] = 1; zb[1] = NZ; 

        /* second block */
	xa[2] = ifw + 1; xb[2] = NX ;
	ya[2] = NY - ifw; yb[2] = NY;
	za[2] = 1; zb[2] = NZ; 
		
	/* calculation area */
	xa[4] = ifw + 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY - ifw - 1;
	za[4] = 1; zb[4] = NZ; 
	}
	
	/* no boundary */
	/* -----------------------------------------------------------------*/
	if ((POS[1]<NPROCX-1)&&(POS[1]>0)&&(POS[2]<NPROCY-1)&&(POS[2]>0)&&(POS[3]>0)&&(POS[3]<NPROCZ-1)){
	
	BLOCK = 0;
			
	/* calculation area */
	xa[4] = 1; xb[4] = NX;
	ya[4] = 1; yb[4] = NY;
	za[4] = 1; zb[4] = NZ; 
	}

}
