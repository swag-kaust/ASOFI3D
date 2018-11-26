/*
 *  Milet
 *   last update 16.07.2003, T. Bohlen
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
	float *pts, ts, tp, sumu, sumpi, val;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
        int  izgrad, ival;
        FILE *fp;
 	


 
/*----------------------SPECIFY HERE MATERIAL PROPERTIES --------------   */

 
        const float VS0=90, VP0=300.0; 		/* velocities of weathering layer */
	const float VS1=150.0, VP1=600.0; 	/* velocities of first layer */
        const float RHO0=1800.0;			/* density of all layers*/
        const float RHO1=2000.0;			/* density of all layers*/

        const float VSR=590.0;        /* refractor velocities */
        const float VPR=950.0;
        const float RHOR=2200.0;			

        const float zgrad=1.6;		/* depth of weathering layer */

	const float Q=10000.0;
	
/*-------------------------------------------------------------------------*/

	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];

      fp=fopen(MFILE,"r");
      if (fp==NULL) err(" Can't open refractor file ! ");
	
	izgrad=iround(zgrad/DH);

	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
		
			fread(&val, sizeof(float), 1, fp);
		   	val=fabs(val); 
          		ival=iround(val/DH);
			
			for (j=1;j<=NYG;j++){
			
				
				/* background */
				Vp=VPR; Vs=VSR; Rho=RHOR; tp=2.0/Q; ts=2.0/Q;
				
				if (j<=ival){
				Vp=VP1; Vs=VS1; Rho=RHO1; tp=2.0/Q; ts=2.0/Q;}

				if (j<=izgrad){
				Vp=VP0; Vs=VS0; Rho=RHO0; tp=2.0/Q; ts=2.0/Q;}

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


/*	sprintf(modfile,"model/milet3D.bin");

	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
*/
	free_vector(pts,1,L);
	fclose(fp);
	
}



