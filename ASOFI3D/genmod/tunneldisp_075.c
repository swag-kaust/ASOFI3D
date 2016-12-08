/*
 *   homogeneous half space 
 *   last update 02.11.02, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DH, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	/*extern char  MFILE[STRING_SIZE];*/
	char  outfile[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float *pts, ts, tp, sumu, sumpi, d, x,y,z;
	int i, j, k, l, ii, jj, kk;
	/*char modfile[STRING_SIZE];*/
	
	/* background */
	const float vp1=5700.0, vs1=3400.0, rho1=2200.0, Qp1=500.0, Qs1=500.0;

	/* tunnel */
	const float vpt=0.0, vst=1.0e-4, rhot=1.25, Qpt=1000.0, Qst=1000.0;

	const float axis_y=(float)NYG*DH/2.0, axis_z=(float)NZG*DH/2.0, radius=3.75;
		    
	/*excevation damage zone*/
	const float vpedz=1800.0*sqrt(3.0), vsedz=1800.0, rhoedz=1500.0, Qpedz=500.0,
	Qsedz=500.0;

	const float radiusedz=10.0;
	
	
	
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
					
				x=(float)i*DH;
				y=(float)j*DH;
				z=(float)k*DH;
				
				/* background */
				Vp=vp1; 
				Vs=vs1; 
				Rho=rho1;
				tp=2.0/Qp1; 
				ts=2.0/Qs1;
						
												
				/*excevation damage zone*/	
				/*tube*/				
				
				d=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z)));
				
				/*if ((d>radius)&&(d<(radius+radiusedz))){ 
					Vp=vpedz+((vp1-vpedz)*(d-radius)/radiusedz); 
					
					
					Vs=vsedz+((vs1-vsedz)*(d-radius)/radiusedz);		
					
					Rho=rhoedz+((rho1-rhoedz)*(d-radius)/radiusedz);
					
					tp=2.0/Qpedz; 
					ts=2.0/Qsedz;}	
				*/
				
				/* tunnel */
				if ((d<radius)){ 
					Vp=vpt; 
					Vs=vst; 
					Rho=rhot;
					tp=2.0/Qpt; 
					ts=2.0/Qst;} 
				
					
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


	/*
	sprintf(outfile,"model/test.taus");
	writemod(outfile,taus,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);
	*/
	/*
	sprintf(outfile,"model/test.taup");
	writemod(outfile,taup,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	sprintf(outfile,"model/test.u");
	writemod(outfile,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	*/
	sprintf(outfile,"model/tunneldisp_075.rho");
	writemod(outfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);
	
	/*
	sprintf(outfile,"model/test.pi");
	writemod(outfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	*/



	free_vector(pts,1,L);
	
	
}



