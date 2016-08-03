/*
 *   homogeneous space with a tunnel
 *   last update 05.02.99, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
 
	extern float DT, DH, *FL, TAU, TS, REC_ARRAY_DEPTH, REFREC[4];
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float x, y, z, d;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	

	 
	/* background */
	const float vp1=5700.0, vs1=3400.0, rho1=2200.0, Qp1=500.0, Qs1=500.0;

	/* tunnel */
	const float vpt=0.0, vst=1.0e-4, rhot=1.25, Qpt=1000.0, Qst=1000.0;

	float axis_y=(float)NGY*DH/2.0, axis_x=(float)NGX*DH, 
	            axis_z=(float)NGZ*DH/2.0, radius=10;

	/* reflector */
	/*const float vp2=4000.0, vs2=2400.0, rho2=1800.0, Qp2=100.0, Qs2=100.0;
	const float xs=20.0;*/
	
	/*excevation damage zone*/
	const float vpedz=1800.0*sqrt(3), vsedz=1800.0, rhoedz=1500.0, Qpt=1000.0, Qst=1000.0;

	float axis_y=(float)NYG*DH/2.0, axis_x=(float)NXG*DH, 
	            axis_z=(float)NZG*DH/2.0, radiusedz=10;
	
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
				
				/* background */
				Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;
				
				
				x=(float)i*DH;
				y=(float)j*DH;
				z=(float)k*DH;
				
			
				
				/* create tunnel */
				/* distance from axis of the tunnel */

				d=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z)));
				if ((d<radius)){ 
					Vp=vpt; Vs=vst; Rho=rhot;tp=2.0/Qpt; ts=2.0/Qst;}
					
				if ((d>radius)&&(d<radius+radiusedz)){ 
					Vp=vpedz+(vpt-vpedz)/10*(d-10); Vs=vsedz+(vst-vsedz)/10*(d-10);
					Rho=rhoedz+(rhot-rhoedz)/10*(d-10);tp=2.0/Qpt; ts=2.0/Qst;}	
					
				/* reflector */
				/*if (x>=x2+xs){
					Vp=vp2; Vs=vs2; Rho=rho2;tp=2.0/Qp2; ts=2.0/Qs2;}*/
	
							
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



	writemod(MFILE,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(MFILE,3);

	free_vector(pts,1,L);
	
}



