/*------------------------------------------------------------------------
 * Averaging of material parameters
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"


void av_mat(float *** rho, 
        float *** C44, float *** C55, float *** C66,
		float *** taus,
		float  *** C66ipjp, float *** C44jpkp, float *** C55ipkp, float *** tausipjp,
		float  *** tausjpkp, float  *** tausipkp, float  *** rjp, float  *** rkp, float  *** rip ) {


	extern int NX, NY, NZ, MYID, L;
	extern FILE *FP;
	double time1=0.0, time2=0.0;
	int i, j, k;

	if (MYID==0){
		fprintf(FP,"\n\n **Message from av_mat (written by PE %d):",MYID);
		fprintf(FP,"\n Averaging of material parameters ... \n");
		time1=MPI_Wtime();
	}

	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){


				/* harmonic averaging of shear modulus */
				C66ipjp[j][i][k]=4.0/((1.0/C66[j][i][k])+(1.0/C66[j][i+1][k])+(1.0/C66[j+1][i+1][k])+(1.0/C66[j+1][i][k]));
				C44jpkp[j][i][k]=4.0/((1.0/C44[j][i][k])+(1.0/C44[j][i][k+1])+(1.0/C44[j+1][i][k+1])+(1.0/C44[j+1][i][k]));
				C55ipkp[j][i][k]=4.0/((1.0/C55[j][i][k])+(1.0/C55[j][i][k+1])+(1.0/C55[j][i+1][k+1])+(1.0/C55[j][i+1][k]));

				/* arithmetic averaging of TAU for S-waves and density */

				if (L){
					tausipjp[j][i][k]=0.25*(taus[j][i][k]+taus[j][i+1][k]+taus[j+1][i+1][k]+taus[j+1][i][k]);
					tausjpkp[j][i][k]=0.25*(taus[j][i][k]+taus[j+1][i][k]+taus[j+1][i][k+1]+taus[j][i][k+1]);
					tausipkp[j][i][k]=0.25*(taus[j][i][k]+taus[j][i+1][k]+taus[j][i+1][k+1]+taus[j][i][k+1]);
				}

				rjp[j][i][k]=0.5*(rho[j][i][k]+rho[j+1][i][k]);
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




















