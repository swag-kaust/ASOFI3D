#include "fd.h"


/*
 * Write local part of the model quantity `q` to the file `modfile`.
 * 'Local part' means 'belonging to the current MPI process'.
 *
 * Parameters
 * ----------
 * modfile :
 *     Prefix of the name of the model file, e.g., 'model/test_'.
 * q :
 *     Quantity being written to file (density, stiffness elements, etc.).
 * format :
 *     File format.
 *
 * See also
 * --------
 * writedsk    To see possible file formats.
 *
 */


void writemod(char modfile[STRING_SIZE], float ***q, int format){
	// External (global) variables.
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
		writedsk(fpmod, q[j][i][k],format);
				
	fclose(fpmod);


}


