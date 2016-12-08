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
	float Vp, Vs, Rho;
	float y, depth;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
        FILE *fp;
	char modfile[STRING_SIZE];
	

	/* mantle */
	const float vp2=8000.0, vs2=4000.0, rho2=3190.0, Qp2=1000.0, Qs2=1000.0;

	/* crust */
	const float vp1=6200.0, vs1=3000.0, rho1=2900.0, Qp1=500.0, Qs1=500.0;

	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];


       	printf(" PE %d: opening file %s ...",MYID,MFILE);

        fp=fopen(MFILE,"r");
        if (fp==NULL) err(" could not open file MFILE ");



	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
		
			fscanf(fp,"%e", &depth);
			depth*=1000.0;
			
			for (j=1;j<=NYG;j++){
				
				y=(float)j*DH;

				if (y<=depth){
					Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;}
				else{
					Vp=vp2; Vs=vs2; Rho=rho2; tp=2.0/Qp2; ts=2.0/Qs2;}
				


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

/*	sprintf(modfile,"model/iceland.bin");

	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
*/
	free_vector(pts,1,L);
	
}



