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
 * This is function initsour.
   Purpose: Computing Position of PE which includes the source position,
   and computing source point in the local grid of that PE
   
----------------------------------------------------------------------*/

#include "fd.h"

int initsour(int nxs,int nys, int  nzs, int *nxsl,int *nysl, int  *nzsl )
{

	extern  int	IENDX, IENDY, IENDZ;
	extern  int   MYID, POS[4];
	extern float PLANE_WAVE_DEPTH;
	extern FILE *FP;

	int npsp;

	/* init of the source coordinates and using nsps as root processor*/
	npsp = -1;

	if ((POS[1]==((nxs-1)/IENDX)) && (POS[2]==((nys-1)/IENDY)) && (POS[3]==((nzs-1)/IENDZ)))
		npsp = MYID;
		
	if (npsp == MYID){
		*nxsl=(nxs)-POS[1]*IENDX;
		*nysl=(nys)-POS[2]*IENDY;
		*nzsl=(nzs)-POS[3]*IENDZ;
		if (PLANE_WAVE_DEPTH <= 0.0) {
			fprintf(FP,"\n **Message from initsource (printed by PE %d):\n",MYID);
			fprintf(FP," PE which includes source (npsp) is me.\n");
			fprintf(FP," Gridpoint of source position within subarray of PE %d :\n",npsp);
			fprintf(FP," nxs= %d \t nys= %d \t nzs= %d \n",*nxsl, *nysl, *nzsl);
		}	

	}
		

	return npsp;

}
