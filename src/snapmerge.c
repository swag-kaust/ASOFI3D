/*------------------------------------------------------------------------
 *   loop over snapshotfiles which have to be merged.                                   
 *
 *  ----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

#include "fd.h"
#include "globvar.h"      /* definition of global variables  */

void _usage() {
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
printf(" Syntax example if excecuted from ./par directory: ../bin/snapmerge in_and_out/sofi3D.json \n");
}

int main(int argc, char **argv) {
int nsnap;
char *fileinp="";
//FILE *FP;


_usage();
if (argc != 2) {
    exit(1);
}
fileinp = argv[1];
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

nsnap=1+floor((TSNAP2-TSNAP1)/TSNAPINC);
fprintf(FP,"nsnap = %d\n", nsnap);
printf("Number of snapshots to be saved: nsnap = %d\n", nsnap);

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
