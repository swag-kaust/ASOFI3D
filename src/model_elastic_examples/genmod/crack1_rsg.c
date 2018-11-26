/*
 *   homogeneous space with a tunnel
 *   last update 05.02.99, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DH, *FL, TAU, TS;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float x, y, z;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	

	/* background */
	const float vp0=6000.0, vs0=3000.0, rho0=2000.0;
	/* position of chrack */
	const float vpt=0.0, vst=0.0, rhot=1e-3;
	const float xc=5.0, yc=5.0, zc=5.0, xd=15.0, yd=10.0, zd=10.0;

	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

	ts=TAU;  
	tp=TAU;

	ws=2.0*PI*FL[1];

	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				
				Vp=vp0; Vs=vs0; Rho=rho0;
				
				
				x=(float)i*DH;
				y=(float)j*DH;
				z=(float)k*DH;
				

				if ((x>xc)&&(x<(xc+xd))&&(y>yc)&&(y<(yc+yd))&&(z>zc)&&(z<(zc+zd))){
					Vp=vpt; Vs=vst; Rho=rhot;}
				

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

	sprintf(modfile,"model/model.bin");

	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);

	free_vector(pts,1,L);
}



