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
 *   loop over snapshotfiles which have to be merged.                                   
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"      /* definition of global variables  */


int main(int argc, char **argv){

int nsnap;
char *fileinp="";
//FILE *FP;

fileinp = argv[1];
printf(" ***********************************************************\n");
printf(" This is program SNAPMERGE. \n");
printf(" Merge of snapshot files from the parallel  \n 3-D Viscoelastic Finite Difference Modelling      \n");
printf("                                                            \n");
printf(" written by  T. Bohlen                          \n");
printf(" Geophysical Institute, Department of Physics,         \n");
printf(" Institute of Technology, Karlsruhe, Germany         \n");
printf(" http://www.gpi.kit.edu \n");
printf(" ***********************************************************\n");
printf("\n");
//printf(" Syntax if excecuted from ./par directory: ../bin/snapmerge in_and_out/sofi3D.inp \n");
printf(" Syntax example if excecuted from ./par directory: ../bin/snapmerge in_and_out/sofi3D.json \n");
printf(" Input file for the snapmerge process from command line : %s \n",fileinp);

if ((FP=fopen(fileinp,"r"))==NULL) err(" Opening input file failed.");
else printf(" Opening input file was successful.\n\n");

/* read parameters from parameter-file */

//read json formated input file
read_par_json(stdout, fileinp);
fclose(FP);


NXG=NX;
NYG=NY;	
NZG=NZ;	
NX = NXG/NPROCX;
NY = NYG/NPROCY;
NZ = NZG/NPROCZ;

nsnap=1+iround((TSNAP2-TSNAP1)/TSNAPINC);

/*printf("NX = %i, NY = %d, NZ = %d",NX,NY,NZ);*/
FP=stdout;

	switch(SNAP){
	case 1 : /*particle velocity*/
		merge(nsnap,1);
		merge(nsnap,2);
		merge(nsnap,3);
		break;
	case 2 : /*pressure */
		merge(nsnap,6);
		break;
	case 4 : /*particle velocity and pressure*/
		merge(nsnap,1);
		merge(nsnap,2);
		merge(nsnap,3);	
		merge(nsnap,6);	
	case 3 :/*energy*/
		merge(nsnap,4);
		merge(nsnap,5);
		break;
	default :
		warning(" snapmerge: cannot identify content of snapshot !");
		break;

	}	
return 0;	

}
