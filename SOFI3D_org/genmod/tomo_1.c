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
	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float y, depth;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
        FILE *fp, *fp_vp, *fp_vs, *fp_rho;
	char modfile[STRING_SIZE];
	char fname_vp[STRING_SIZE], fname_vs[STRING_SIZE], fname_rho[STRING_SIZE];
	

	/* water */
	const float h1=20.0, vp1=1450.0, vs1=0.0, rho1=1000.0, Qp1=10000.0, Qs1=10000.0;

	/* second layer is read from file */
	const float h2=10.0, vp2=1550.0, vs2=150.0, rho2=1400.0, Qp2=500.0, Qs2=500.0;
	
	/* third layer */
	const float h3=15.0, vp3=1600.0, vs3=350.0, rho3=1600.0, Qp3=500.0, Qs3=500.0;

	/* half space */
	const float vp4=1640.0, vs4=500.0, rho4=1900.0, Qp4=500.0, Qs4=500.0;

	
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];

	sprintf(fname_vp,"%s.vp",MFILE);
	sprintf(fname_vs,"%s.vs",MFILE);
 	sprintf(fname_rho,"%s.rho",MFILE);
	
	
/*      	printf(" PE %d: opening model file %s ...",MYID,fname_vp);
        fp_vp=fopen(fname_vp,"rb");
        if (fp_vp==NULL) err(" could not open file vp-file ");

      	printf(" PE %d: opening model file %s ...",MYID,fname_vs);
        fp_vs=fopen(fname_vs,"rb");
        if (fp_vs==NULL) err(" could not open file vs-file ");

      	printf(" PE %d: opening model file %s ...",MYID,fname_rho);
        fp_rho=fopen(fname_rho,"rb");
        if (fp_rho==NULL) err(" could not open file rho-file ");
*/

	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
		
 			for (j=1;j<=NYG;j++){
				
				y=(float)j*DH;

				if (y<=h1){
					Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;}
				else if (y<=h1+h2){
/*					fread(&Vp,4,1,fp_vp);
					fread(&Vs,4,1,fp_vs);
    					fread(&Rho,4,1,fp_rho);
*/
					Vp=vp2; Vs=vs2; Rho=rho2;
					tp=2.0/Qp2; ts=2.0/Qs2;
				}
				else if (y<=h1+h2+h3){
					Vp=vp3; Vs=vs3; Rho=rho3; tp=2.0/Qp3; ts=2.0/Qs3;}
				else {
					Vp=vp4; Vs=vs4; Rho=rho4; tp=2.0/Qp4; ts=2.0/Qs4;}
				
				


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

 /*       fclose(fp_vp);
        fclose(fp_vs);
        fclose(fp_rho);
*/
	sprintf(modfile,"model/tomo_1.bin");

	writemod(modfile,u,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);

	free_vector(pts,1,L);
	
}



