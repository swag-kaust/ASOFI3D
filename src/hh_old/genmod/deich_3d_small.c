/*
*   Dike-Model-3D
*   D.Koehn
*   Kiel, 8.11.2003
*
*   update: 16.11.2003 - added "groundwater-line-polynom"
*           20.12.2003 - fixed polynomial coefficient bug
*           18.03.2004 - NR don't run on Altix - calculate Polynomial Coefficients with MATLAB !
*           28.03.2004 - included dike-shape-function f3
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
	float Vp, Vs, Rho, xi, yi, zi;
	float *pts, ts, tp, sumu, sumpi;
	float f1, f2;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	
        
	/* some mathematical constants*/
        const float pim = 4.0*atan(1.0);
	const float pic = pim/180.0;
	
	/* geometrical parameters*/
	const float alpha=36.32*pic, delta=24.52*pic;
	const float l1=1.48, l2=1.0, l3=2.485, h=1.19;
	const float dp=0.05, da=0.04, T=0.3;
	

        /* parameters for air */
        /*const float Vp0=0.0, Vs0=0.0001, Rho0=1.25;*/
        const float Vp0=0.0, Vs0=0.0001, Rho0=1.25;

	/* parameters for water*/
	//const float Vp1=1500.0, Vs1=1.0e-4, Rho1=1000.0;

	/* parameters for Lehmkern  */
	const float Vp2=500.0, Vs2=250.0, Rho2=1800.0;

        /* parameters for Plexiglas (PMMA)*/
	/*const float Vp2w=3150.0, Vs2w=1570.0, Rho2w=1180.0;*/
	/*const float Vp2w=0.0, Vs2w=0.0001, Rho2w=1.25;*/
	const float Vp2w=0.0, Vs2w=0.0001, Rho2w=1.25;
	
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
				
				xi=(float)i*DH;
				yi=(float)j*DH;
				zi=(float)k*DH;
				
				/* fill model with clay*/
                                Vp=Vp2; Vs=Vs2; Rho=Rho2;

                                /* area lambda 1 */
                                /* ------------ */
                                /* calculate "dike-shape-function - f1"        */
				/* to describe the boundary of the dike on landside. */

				f1=(-tan(alpha)*xi)+h+(da);

                                /* air */
                                if((yi<=f1)){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* area lambda 2 */
				/* ------------ */

                                /* calculate "dike-shape-functions-f2" */
				/* to describe the boundary of the dike on seaside. */

				f2 = (tan(delta) * xi) - (tan(delta)*((2*dp)+l1+l2+l3+da)) + h;
                                
			          
                                /* air */
                                if(yi<f2){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* boundaries*/
				
				/* fill front with air */
				if(zi<(da+dp)){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* front boundary */
				if((zi<(da+dp))&&(zi>da)){
                                Vp=Vp2w; Vs=Vs2w; Rho=Rho2w;
				}
				
				/* background boundary */
				/* fill back with plexiglas*/
				if(zi>(da+dp+T)){
                                Vp=Vp2w; Vs=Vs2w; Rho=Rho2w;
				}
				
				
				/* lower boundary */
				if((yi>(da+h-dp))){
                                Vp=Vp2w; Vs=Vs2w; Rho=Rho2w;
				}
				
				/* left boundary */
                                if((xi<=da+dp)){
                                Vp=Vp2w; Vs=Vs2w; Rho=Rho2w;
				}
				
				/* right boundary */
                                if((xi>=da+dp+l1+l2+l3)){
                                Vp=Vp2w; Vs=Vs2w; Rho=Rho2w;
				}

				/* upper boundary*/
				if((yi<=da)){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* lower boundary */
				if((yi>=(h+da))){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* left boundary */
                                if((xi<=da)){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* right boundary */
                                if((xi>da+dp+l1+l2+l3+dp)){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* fill back with air */
				if(zi>(da+dp+T+dp)){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
				}
				
				/* fill back with air */
				if(k<4){
                                Vp=Vp0; Vs=Vs0; Rho=Rho0;
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



