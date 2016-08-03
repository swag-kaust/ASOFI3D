/* $Id: benchmark_acoustic.c,v 1.1 2007/11/03 01:07:37 koehn Exp $ */
/*
 *   jastram model
 */

#include "fd.h"

void model_acoustic(float  ***  rho, float ***  pi){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
 
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], MYID;
	extern char  MFILE[STRING_SIZE];
	extern float DX, DY;
	
	/* local variables */
	float piv;
	float Vp, Rho;
	int i, j, k, ii, jj, kk;
	char filename[STRING_SIZE];
		
	/*-----------------------------------------------------------------------*/

 	for (j=1;j<=NYG;j++){
		for (i=1;i<=NXG;i++){
			for (k=1;k<=NZG;k++){
				
				
				Vp=3500.0; Rho=2000.0;
	
				piv=Vp*Vp*Rho;

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
		}
	}	

	


	

	sprintf(filename,"%s.fdmpi.pi",MFILE);

	writemod(filename,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
		
}



