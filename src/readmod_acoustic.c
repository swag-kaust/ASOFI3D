/* ----------------------------------------------------------------------
 *   homogeneous space with a tunnel
 *
 ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"


void readmod_acoustic(float  ***  rho, float ***  pi, int ishot){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], MYID;
	extern char  MFILE[STRING_SIZE];
	extern FILE *FP;

	/* local variables */
	float piv;
	float Vp, Rho;
	float * rhobuffer, * vpbuffer;
	int i, j, k, ii, jj, kk, N, ind;
	FILE  *fp_vp=NULL, *fp_rho=NULL;
	char filename[STRING_SIZE];
	char fname_vp[STRING_SIZE], fname_rho[STRING_SIZE];
	double 	time1=0.0, time2=0.0, time3=0.0, time4=0.0;

	/*-----------------------------------------------------------------------*/



	sprintf(fname_vp,"%s_shot%d.vp",MFILE,ishot);
	sprintf(fname_rho,"%s_shot%d.rho",MFILE,ishot);

	if (MYID==0) {
		printf(" PE %d: opening model file %s ...",MYID,fname_vp);
		fp_vp=fopen(fname_vp,"rb");
		if (fp_vp==NULL) err(" could not open file vp-file ");

		printf(" PE %d: opening model file %s ...",MYID,fname_rho);
		fp_rho=fopen(fname_rho,"rb");
		if (fp_rho==NULL) err(" could not open file rho-file ");
	}



	/* loop over global grid */

	if (MYID==0) time1=MPI_Wtime();

	N=NYG*NXG;
	rhobuffer=(float *)malloc(N*4);
	if (rhobuffer==NULL) err(" allocation failure rhobuffer!");
	vpbuffer=(float *)malloc(N*4);
	if (vpbuffer==NULL) err(" allocation failure vpbuffer!");

    size_t N_actual;

	for (k=1;k<=NZG;k++){
		if (MYID==0) {
			time2=MPI_Wtime();
			fprintf(FP," reading k=%d\n",k);
			fprintf(FP," real time : %4.2f s.\n",time2-time1);
			time1=time2;

			N_actual = fread(vpbuffer, 4, N, fp_vp);
            if (N_actual != (size_t) N) {
                warning("Possible error while reading model for 'Vp'");
            }

			N_actual = fread(rhobuffer, 4, N, fp_rho);
            if (N_actual != (size_t) N) {
                warning("Possible error while reading model for 'Vp'");
            }

		}

		if (MYID==0){
			fprintf(FP," broadcasting buffer arrays %d bytes...",2*N*4);
			time3=MPI_Wtime();
		}
		MPI_Bcast(vpbuffer,N,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(rhobuffer,N,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		time4=MPI_Wtime();
		if (MYID==0){
			time4=MPI_Wtime();
			fprintf(FP," real time : %4.2f s.\n",time4-time3);
		}
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				ind=(i-1)*NYG +j;
				Vp=vpbuffer[ind-1];
				Rho=rhobuffer[ind-1];

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

	if (MYID==0) {
		fclose(fp_vp);
		fclose(fp_rho);
	}

	free(rhobuffer);
	free(vpbuffer);



	sprintf(filename,"%s.SOFI3D.rho",MFILE);

	writemod(filename,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);


}



