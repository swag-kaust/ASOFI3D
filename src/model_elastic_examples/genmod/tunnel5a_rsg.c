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
	float x, y, z, d, dip_rad, yr;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	

	/* background  (beneath reflector)*/
	const float vp2=4000.0, vs2=2400.0, rho2=1800.0, Qp2=100.0, Qs2=100.0;

	/* tunnel */
	const float vpt=0.0, vst=0.0, rhot=1e-3, Qpt=1000.0, Qst=1000.0;
	const float x1=0.0, x2=(float)NXG*DH, axis_y=25.0, axis_z=25.0, radius=5.0;

	/* reflector */
	/*intersection point with tunnel */
	const float xs=90.0, ys=30.0, dip=35.0; /* dip in degrees from vertical */
	/* Parameters around tunnel */
	const float vp1=5700.0, vs1=3400.0, rho1=2200.0, Qp1=500.0, Qs1=500.0;

	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];

	dip_rad=dip*PI/180.0;

	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				
				/* background */
				Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;
				
				
				x=(float)i*DH;
				y=(float)j*DH;
				z=(float)k*DH;
				
					
				/* dipping reflector */
				yr=tan(dip_rad)*(xs-x);
				if (y>yr+axis_y+radius){
					Vp=vp2; Vs=vs2; Rho=rho2;tp=2.0/Qp2; ts=2.0/Qs2;}
				


				/* create tunnel */
				/* distance from axis of the tunnel */
/*				d=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z)));
				if ((d<radius)&&(x>=x1)&&(x<=x2)){ 
					Vp=vpt; Vs=vst; Rho=rhot;tp=2.0/Qpt; ts=2.0/Qst;}
*/					
							
							
				
				
				
/*				if ((x>x1)&&(x<x2)&&(y>axis_y-radius)&&(y<axis_y+radius)&&(z>axis_z-radius)&&(z<axis_z+radius))
				{ Vp=vpt; Vs=vst; Rho=rhot; tp=2.0/Qpt; ts=2.0/Qst;}
*/
				if (y<axis_y+radius){ 
					Vp=vpt; Vs=vst; Rho=rhot;tp=2.0/Qpt; ts=2.0/Qst;}

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


/*	sprintf(modfile,"model/model.bin");

	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
*/
	free_vector(pts,1,L);
	
}



