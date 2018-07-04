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

	
