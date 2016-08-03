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
/*-------------------------------------------------------------
  * Cat seismograms (collect seismogram portions from each PE for collective output)
  *
  *------------------------------------------------------------- */

#include "fd.h"

void	catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns) {

	//extern FILE *FP;

	int		i, j, k;
	float		**fulldata2;

	//fprintf(FP,"\n **Message from function catseis:\n");
	//fprintf(FP,"\n Allocating memory \n");
	/* temporary global data array for MPI-exchange */
	fulldata2 = fmatrix(1,ntr_glob,1,ns);

	k = 0;	/* trace counter for local data array */

	//fprintf(FP," Start loop over ntr_glob = %d traces with each ns = %d samples \n",ntr_glob,ns);
	/* loop over global traces: copy traces of local array	*/
	/* to appropriate locations in the global array		*/
	for(i=1;i<=ntr_glob;i++)
	{
		if (recswitch[i]) {
			k++;
			for(j=1;j<=ns;j++)	fulldata2[i][j] = data[k][j];
		}
	}

	//fprintf(FP," Start MPI_Allreduce \n");
	MPI_Allreduce(&fulldata2[1][1], &fulldata[1][1], ntr_glob*ns, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	//fprintf(FP," De-allocating memory \n");
	free_matrix(fulldata2, 1,ntr_glob,1,ns);
}
