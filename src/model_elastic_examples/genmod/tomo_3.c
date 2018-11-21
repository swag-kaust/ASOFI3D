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
	int i, j, k, l, ii, jj, kk;
        FILE *fp, *fp_vp, *fp_vs, *fp_rho;
	char modfile[STRING_SIZE];
	char fname_vp[STRING_SIZE], fname_vs[STRING_SIZE], fname_rho[STRING_SIZE];
	double 	time1, time2;
		
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
	
	
	if (MYID==0) printf(" PE %d: opening model file %s ...",MYID,fname_vp);
        fp_vp=fopen(fname_vp,"rb");
        if (fp_vp==NULL) err(" could not open file vp-file ");

       	if (MYID==0) printf(" PE %d: opening model file %s ...",MYID,fname_vs);
        fp_vs=fopen(fname_vs,"rb");
        if (fp_vs==NULL) err(" could not open file vs-file ");

      	if (MYID==0) printf(" PE %d: opening model file %s ...",MYID,fname_rho);
        fp_rho=fopen(fname_rho,"rb");
        if (fp_rho==NULL) err(" could not open file rho-file ");


	tp=TAU; ts=TAU;
	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

	/* loop over global grid */
		
	if (MYID==0) time1=MPI_Wtime();
		
 	for (j=1;j<=NYG;j++){
		if (MYID==0) {
			time2=MPI_Wtime();
			fprintf(FP," reading j=%d\n",j);
			fprintf(FP," real time : %4.2f s.\n",time2-time1);
			time1=time2;
		}
		for (i=1;i<=NXG;i++){
			for (k=1;k<=NZG;k++){
				fread(&Vp,4,1,fp_vp);
				fread(&Vs,4,1,fp_vs);
    				fread(&Rho,4,1,fp_rho);
				Rho*=1000.0;
				
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

        fclose(fp_vp);
        fclose(fp_vs);
        fclose(fp_rho);

	sprintf(modfile,"model/tomo3.bin");

/*	writemod(modfile,u,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
*/
	free_vector(pts,1,L);
	
}



