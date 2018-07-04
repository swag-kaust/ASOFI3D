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


