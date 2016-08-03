/*
 *   homogeneous space with a tunnel
 *   last update 05.02.99, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
 
	extern float DT, DH, *FL, TS;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs;
	float y;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
        FILE *fp;
	char modfile[STRING_SIZE];
	extern FILE *FP;
	
	
	/* parameters for layer 2 */
	const float vp=2600.0, vpvs=1.73, Rho=2300.0, h1=70.0, Qp=10000.0, Qs=10000.0;

	
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];

			
        fp=fopen(MFILE,"rb");
        if (fp==NULL) err(" could not open file vp-file ");


	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		if (MYID==0) fprintf(FP," raeding k=%d\n",k);
		for (i=1;i<=NXG;i++){
		
 			for (j=1;j<=NYG;j++){
					
				y=(float)j*DH;
				tp=2.0/Qp; ts=2.0/Qs;
					
				fread(&Vp,4,1,fp);
				
				if (y<=h1) Vp=vp;
					
				Vs=vp/vpvs;


				sumu=0.0; 
				sumpi=0.0;
				for (l=1;l<=L;l++){
					sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
					sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
				}

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

        fclose(fp);
 	sprintf(modfile,"model/mallik.bin");

	writemod(modfile,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);

	free_vector(pts,1,L);
	
}



