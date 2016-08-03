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
//	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float x,y,z, dvs, r2;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	

	/* crust */
	const float h1=20.0e3, vp1=7100.0, vs1=3890.0, rho1=3000.0, Qp1=1000.0, Qs1=1000.0;

	/* mantle */
	const float vp2=8000.0, vs2=4380.0, rho2=3150.0, Qp2=1000.0, Qs2=1000.0;
	
	/* plume */
	const float xc=(float)NXG*DH/2.0, zc=(float)NZG*DH/2.0, ybase=400.0e3, ytop=40.0e3;
	const float ll=175.0e3, dvsmax=-4.2*vs2/100.0, vpvs=1.8257;

	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];



	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
		
			for (j=1;j<=NYG;j++){
				
				x=(float)i*DH;
				y=(float)j*DH;
				z=(float)k*DH;

				if (y<=h1){
					Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;}
					
				else if (y<=ytop){
					Vp=vp2; Vs=vs2; Rho=rho2; tp=2.0/Qp2; ts=2.0/Qs2;}
				
				else if (y<=ybase){
					r2=((x-xc)*(x-xc))+((z-zc)*(z-zc));
					dvs=dvsmax*exp(-r2/(ll*ll));

					Vs=vs2+dvs;
					Vp=Vs*vpvs;
					Rho=rho2; tp=2.0/Qp2; ts=2.0/Qs2;}
				else {
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


	sprintf(modfile,"model/plume3.bin");

	writemod(modfile,u,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);

	free_vector(pts,1,L);
	
}



