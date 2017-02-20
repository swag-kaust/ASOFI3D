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
 *   write local model to file              
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void writemod(char modfile[STRING_SIZE], float *** rho, int format){


	/* extern variables */
	/*extern int MYID;*/
	extern int NX, NY, NZ, POS[4], IDX, IDY, IDZ;


	int i, j, k;
	FILE *fpmod;
	char file[STRING_SIZE];

	/*printf("\n\n PE %d is writing model to \n",MYID);*/
	sprintf(file,"%s.%i%i%i",modfile,POS[1],POS[2],POS[3]);
	/*printf("\t%s\n\n", file);*/
	fpmod=fopen(file,"w");
	for (k=1;k<=NZ;k+=IDZ)
	for (i=1;i<=NX;i+=IDX)
	for (j=1;j<=NY;j+=IDY)
		writedsk(fpmod,rho[j][i][k],format);
				
	fclose(fpmod);


}


