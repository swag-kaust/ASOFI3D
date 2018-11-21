/*
 *   3-D reservoir model 1
 *   last update 06.02.00, T. Bohlen
 */

#include "fd.h"

void model(float  ***  prho, float ***  ppi, float ***  pu){

	/* extern variables */
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], MYID;
	extern float DH;

	/* local variables */
	float rho, vp, vs, y;
	int i, j, k, ii, jj, kk;
	FILE *fp;
	char modfile[STRING_SIZE];



	/*----------------------SPECIFY HERE MATERIAL PROPERTIES ------------------*/
	/* paramemeters for computation of density and Vs */
	const float vp0=3000.0, y0=90.0;

	/*-----------------------------------------------------------------------*/


	fp=fopen("random3d2.vp","r");
	if (fp==NULL) err(" Could not open model file ! ");

	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
				y=(float)j*DH;
				fread(&vp, sizeof(float), 1, fp);
				/*printf("k= %d \t i= %d \t j= %d \t vp= %e \n",k,i,j,vp);*/
				if (y<y0) vp=vp0;
				/*!!!!!!!!*/
				/*vp=vp0;*/
				/*!!!!!!!!*/
				/*vs=-314.59 + 0.61*vp;*/  /*unconsolidated sandstone, Blangy (1992) */
				/*rho=1498.0 + 0.22*vp;*/
				
				vs=2000.0;
				rho=2100.0;
								
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){

					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					pu[jj][ii][kk]=vs*vs*rho;
					prho[jj][ii][kk]=1.0/rho;
					ppi[jj][ii][kk]=vp*vp*rho;
				}
			}
		}
	}
	
	fclose(fp);	
	sprintf(modfile,"model.bin");

	/* each PE writes his model to disk */
/*	writemod(modfile,prho,3); */


	/*  model files are then merged into one file by PE 0 */
	/*if (MYID==0) mergemod(modfile,3); */
}



