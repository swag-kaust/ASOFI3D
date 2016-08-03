/*
 *   homogeneous half space 
 *   last update 02.11.02, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DX, DY, DZ, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	/*extern char  MFILE[STRING_SIZE];*/
	char  outfile[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float *pts, ts, tp, sumu, sumpi, d,x,y,z;
	/*float y_fault, alpha_fault, x_fault;*/
	int i, j, k, l, ii, jj, kk;
	/*char modfile[STRING_SIZE];*/
	
	/* background */
	const float vphard=6000.0, vshard=3000.0, rhohard=2200.0, Qphard=500.0, Qshard=500.0;
	

	/* tunnel */
	/*const float vpt=0.0, vst=1.0e-4, rhot=1.25, Qpt=1000.0, Qst=1000.0, tunnel_length=60;*/

	const float axis_y=(float)NYG*DH/2.0, axis_z=(float)NZG*DH/2.0, radius=5.0;
		    
	/*reflector */
	/*const float vp_fault=4000.0, vs_fault=2400.0, rho_fault=1800.0, Qp_fault=100.0, Qs_fault=100.0;
	const float x_fault90=100.0, alpha_f=60;*/
	
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
				
				/* calculating actual meters from gridpoints */					
				x=(float)i*DX;
				y=(float)j*DY;
				z=(float)k*DZ;
				
				/* background hardrock*/
				Vp=vphard; Vs=vshard; Rho=rhohard;
				tp=2.0/Qphard; ts=2.0/Qshard;
				
				
				/* inserting tunnel 
				if ((d<radius) && x<tunnel_length){ 
					Vp=vpt; Vs=vst; Rho=rhot;
					tp=2.0/Qpt; ts=2.0/Qst;}*/
					
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



        /*output of model according to IDX,IDY,IDZ from inputfile*/
	sprintf(outfile,"model/trogmodel.piv");
	writemod(outfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	
	sprintf(outfile,"model/trogmodel.rho");
	writemod(outfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	
	free_vector(pts,1,L);
	
	
}



