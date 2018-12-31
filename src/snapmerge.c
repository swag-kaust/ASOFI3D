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
printf(" Syntax example if excecuted from ./par directory: ../bin/snapmerge in_and_out/asofi3D.json \n");
}

static int snapmerge(int argc, char **argv) {
int nsnap;
char *fileinp="";
extern FILE *FP;
extern float TSNAP1;

FILE *fp_param;

_usage();
if (argc != 2) {
    // !!!!!!!!!!!!!!!!!!!! TODO: Add explanation
    exit(1);
}
fileinp = argv[1];

printf("\n");
printf("Input file for the snapmerge from command line: %s\n", fileinp);

fp_param=fopen(fileinp, "r");
if (fp_param != NULL) {
    printf("Opening input file was successful.\n");
    fclose(fp_param);
	fp_param = NULL;
} else {
    err("Opening input file failed.");
}

FP = stdout;

/* read parameters from parameter-file */

//read json formated input file
read_par_json(FP, fileinp);

NXG=NX;
NYG=NY;	
NZG=NZ;	
NX = NXG/NPROCX;
NY = NYG/NPROCY;
NZ = NZG/NPROCZ;

if (TSNAP2 > TIME) {
    fprintf(FP, "WARNING: TSNAP2 = %f is larger than TIME = %f. "
            "Set TSNAP2 = TIME\n", TSNAP2, TIME);
    TSNAP2 = TIME;
}

nsnap=1+floor((TSNAP2-TSNAP1)/TSNAPINC);
fprintf(FP, "Number of snapshots to be saved: nsnap = %d\n", nsnap);

/*printf("NX = %i, NY = %d, NZ = %d",NX,NY,NZ);*/

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
		printf("[%s] Parameter SNAP set to zero in the input file, "
               "therefore, no snapshots were produced during simulation\n",
               __func__);
		break;
	}
return 0;	

}

int main(int argc, char **argv) {
    snapmerge(argc, argv);
}
