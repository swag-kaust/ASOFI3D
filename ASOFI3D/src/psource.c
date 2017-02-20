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
 *   generate P-wave source at source nodes
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, float *** sxx, float *** syy, float *** szz,
float **  srcpos_loc, float ** signals, int nsrc, int * stype){


	extern float DX, DY, DZ;
	//extern float DT;
	extern int NT;

	int i, j, k, l;
	float amp=0.0;



	/* adding source wavelet to stress components 
	   (explosive source) at source points */


	for (l=1;l<=nsrc;l++) {
	    if(stype[l]==1){
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		k=(int)srcpos_loc[3][l];

		//amp=signals[l][nt]; //unscaled explosive source
		//amp=(signals[l][nt])/(DX*DY*DZ); //scaled explosive source with respect to spatial discretization, seismic Moment = 1 Nm

		//scaled explosive source with respect to spatial and temporal discretization, seismic Moment = 1 Nm
		// -> additional temporal derivation
		if(nt==1){amp=signals[l][nt+1]/(2.0*DX*DY*DZ);}
		if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/(2.0*DX*DY*DZ);}
		if(nt==NT){amp=-signals[l][nt-1]/(2.0*DX*DY*DZ);}

		sxx[j][i][k]+=amp;
		syy[j][i][k]+=amp;
		szz[j][i][k]+=amp;
	    }
	}

}
