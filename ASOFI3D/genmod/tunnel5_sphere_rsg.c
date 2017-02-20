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
	float x, y, z, d, dip_rad, yr;
	float xd1, yd1, zd1, xd2, yd2, zd2, xd3, yd3, zd3, r1, r2, r3; 
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	

	/* spheres */
	const float vp2=4000.0, vs2=2400.0, rho2=1800.0, Qp2=100.0, Qs2=100.0, R=2.5;
	
	const float sx1=10.0, sy1=40.0, sz1=REFREC[3];
	const float sx2=25.0, sy2=70.0, sz2=REFREC[3];
	const float sx3=50.0, sy3=55.0, sz3=REFREC[3];
	 

	/* tunnel */
	const float vpt=0.0, vst=0.0, rhot=1.25, Qpt=1000.0, Qst=1000.0;

	float x1=0.0, x2=(float)NXG*DH, axis_y=REFREC[2], 
	            axis_z=REFREC[3], radius=REC_ARRAY_DEPTH;

	/* reflector */
	/*intersection point with tunnel */
	const float xs=70.0,  dip=35.0; /* dip in degrees from vertical */
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


	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				
				/* background */
				Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;
				
				
				x=(float)i*DH;
				y=(float)j*DH;
				z=(float)k*DH;
				
				xd1=x-sx1; yd1=y-sy1; zd1=z-sz1;
				xd2=x-sx2; yd2=y-sy2; zd2=z-sz2;
				xd3=x-sx3; yd3=y-sy3; zd3=z-sz3;
				
				r1=sqrt(xd1*xd1+yd1*yd1+zd1*zd1);
				r2=sqrt(xd2*xd2+yd2*yd2+zd2*zd2);
				r3=sqrt(xd3*xd3+yd3*yd3+zd3*zd3);
				
				if ((r1<=R)||(r2<=R)||(r3<=R)){
					Vp=vp2; Vs=vs2; Rho=rho2;
				} 
				
				
				/* create tunnel */
				/* distance from axis of the tunnel */

				d=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z)));
				if ((d<radius)&&(x>=x1)&&(x<=x2)){ 
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


	sprintf(modfile,"model/tunnel5_sphere.bin");

/*	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
*/
	free_vector(pts,1,L);
	
}



