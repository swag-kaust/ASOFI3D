/*
 *   homogeneous space with a tunnel
 *   last update 05.02.99, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
 
	extern float DT, DH, *FL, TS, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	extern FILE *FP;
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float y, depth;
	float *pts, ts, tp, sumu, sumpi;
	float * rhobuffer, * vpbuffer, * vsbuffer;
	int i, j, k, l, ii, jj, kk, N, ind;
        FILE *fp, *fp_vp, *fp_vs, *fp_rho;
	char modfile[STRING_SIZE];
	char fname_vp[STRING_SIZE], fname_vs[STRING_SIZE], fname_rho[STRING_SIZE];
	double 	time1, time2, time3, time4;
		
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];

	sprintf(fname_vp,"%s.vp",MFILE);
	sprintf(fname_vs,"%s.vs",MFILE);
 	sprintf(fname_rho,"%s.rho",MFILE);
	
	
      	if (MYID==0) {
		printf(" PE %d: opening model file %s ...",MYID,fname_vp);
        	fp_vp=fopen(fname_vp,"rb");
        	if (fp_vp==NULL) err(" could not open file vp-file ");

      		printf(" PE %d: opening model file %s ...",MYID,fname_vs);
        	fp_vs=fopen(fname_vs,"rb");
        	if (fp_vs==NULL) err(" could not open file vs-file ");

      		printf(" PE %d: opening model file %s ...",MYID,fname_rho);
        		fp_rho=fopen(fname_rho,"rb");
        	if (fp_rho==NULL) err(" could not open file rho-file ");
	}


	tp=TAU; ts=TAU;
	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

	/* loop over global grid */
		
	if (MYID==0) time1=MPI_Wtime();

	N=NXG*NZG;
	rhobuffer=(float *)malloc(N*4);	
	if (rhobuffer==NULL) err(" allocation failure rhobuffer!");
	vpbuffer=(float *)malloc(N*4);
	if (vpbuffer==NULL) err(" allocation failure vpbuffer!");
	vsbuffer=(float *)malloc(N*4);
	if (vsbuffer==NULL) err(" allocation failure vsbuffer!");
	
 	for (j=1;j<=NYG;j++){
		if (MYID==0) {
			time2=MPI_Wtime();
			fprintf(FP," reading j=%d\n",j);
			fprintf(FP," real time : %4.2f s.\n",time2-time1);
			time1=time2;
				
			fread(vpbuffer,4,N,fp_vp);
			fread(vsbuffer,4,N,fp_vs);
    			fread(rhobuffer,4,N,fp_rho);
			
		}

		if (MYID==0){ 
			fprintf(FP," broadcasting buffer arrays %d bytes...",3*N*4);
			time3=MPI_Wtime();
		}
		MPI_Bcast(vpbuffer,N,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(vsbuffer,N,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(rhobuffer,N,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		time4=MPI_Wtime();
		if (MYID==0){
			time4=MPI_Wtime();
			fprintf(FP," real time : %4.2f s.\n",time4-time3);
		}
	
		for (i=1;i<=NXG;i++){
			for (k=1;k<=NZG;k++){
				ind=(i-1)*NXG +k;
				Vp=vpbuffer[ind-1];
				Vs=vsbuffer[ind-1];
				Rho=rhobuffer[ind-1]*1000.0;
				
				/*Vp=1900.0; Vs=150.0; Rho=1000.0;*/

			

				muv=Vs*Vs*Rho/(1.0+sumu);
				piv=Vp*Vp*Rho/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					taus[jj][ii][kk]=ts;
					taup[jj][ii][kk]=tp;
					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
		}
	}	

        if (MYID==0) {
		fclose(fp_vp);
        	fclose(fp_vs);
        	fclose(fp_rho);
	}
	
	free(rhobuffer);
	free(vpbuffer);
	free(vsbuffer);

/*	sprintf(modfile,"model/tomo3_rho.bin");
	writemod(modfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);

	sprintf(modfile,"model/tomo3_u.bin");
	writemod(modfile,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);

	sprintf(modfile,"model/tomo3_pi.bin");
	writemod(modfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
*/
	free_vector(pts,1,L);
	
}



