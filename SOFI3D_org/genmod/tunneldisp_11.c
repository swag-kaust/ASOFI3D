/*
 *   homogeneous half space 
 *   last update 02.11.02, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU, TS, DX, DY, DZ;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho, y, x, z;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	char outfile[STRING_SIZE];
	float d, linfactor;
	
	/*tunnel parameters*/	
	const float axis_y=(float)NYG*DY/2.0, axis_z=(float)NZG*DZ/2.0, radius=2.6;
        const float vpt=0.0, vst=1.0e-4, rhot=1.25, Qpt=1000.0, Qst=1000.0;

        /*excevation damage zone parameters in % of original values vp=vp-vp*edz */
        const float vpedz=0.07, vsedz=0.07, rhoedz=0.07, Qpedz=0.0, Qsedz=0.0;

	/* background parameters */
	const float vp1=6000.0, vs1=3000.0, rho1=2000.0;


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

				x=(float)i*DX;
                                y=(float)j*DY;
                                z=(float)k*DZ;

				Vp=vp1; Vs=vs1; Rho=rho1;
 				
				d=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z)));

                                if ((d>radius)&&(d<(3*radius))){
                                        linfactor=(-d/radius+3)*0.5;
                                        Vp=vp1-vp1*(vpedz)*linfactor;
                                        Vs=vs1-vs1*(vsedz)*linfactor;
                                        Rho=rho1-rho1*(rhoedz)*linfactor;

                                        /*Vp=vpedz+((vp1-vpedz)*(d-radius)/radiusedz);
                                        Vs=vsedz+((vs1-vsedz)*(d-radius)/radiusedz);
                                        Rho=rhoedz+((rho1-rhoedz)*(d-radius)/radiusedz);*/

                                        /*if ((x=10)&&(y=40))   {
                                                printf("vp =%5.0f with d=%5.2f and -0.5*(d/radius-3)=%5.2f \n",Vp,d,linfactor);
                                        }
                                        tp=2.0/Qpedz;
                                        ts=2.0/Qsedz;*/
                                }

                                /* tunnel*/
                                if ((d<radius)){
                                        Vp=vpt; Vs=vst; Rho=rhot;
                                        /*tp=2.0/Qpt; ts=2.0/Qst;*/}
 
  				
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


	sprintf(outfile,"model/tunneldisp_10.rho");
        writemod(outfile,rho,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(outfile,3);

        sprintf(outfile,"model/tunneldisp_10.vp");
        writemod(outfile,pi,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(outfile,3);

        sprintf(outfile,"model/tunneldisp_10.vs");
        writemod(outfile,u,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(outfile,3);

	free_vector(pts,1,L);
}



