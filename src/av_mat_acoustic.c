/*------------------------------------------------------------------------
 * Averaging of material parameters for acoustic wave modelling
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"


void av_mat_acoustic(float *** rho, float  *** rjp, float  *** rkp, float  *** rip ){


	extern int NX, NY, NZ, MYID;
	extern FILE *FP;
	double time1=0.0, time2=0.0;
	int i, j, k;
	
	if (MYID==0){
		fprintf(FP,"\n\n **Message from av_mat_acoustic (written by PE %d):",MYID);
		fprintf(FP,"\n Averaging of material parameters ... \n");
		time1=MPI_Wtime();
	}	

	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){

	       
	       
	       /* arithmetic averaging of density */

	       rjp[j][i][k]=0.5*(rho[j][i][k]+rho[j+1][i][k+1]);
	       rkp[j][i][k]=0.5*(rho[j][i][k]+rho[j][i][k+1]);
	       rip[j][i][k]=0.5*(rho[j][i][k]+rho[j][i+1][k]);


			}
		}
	}



	if (MYID==0){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}



}




















