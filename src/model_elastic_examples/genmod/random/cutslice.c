/*
 *  Cut vertical 2-D slice out of 3-D volume 
 *   last update 05.01.01, T. Bohlen
 */

/* files to include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "util.c"

int main(int argc, char **argv){


	const int NX=240, NY=600, NZ=240;
	const int zslice=120;

	/* local variables */
	
	float *** volume, std, vp;
	double sumvp, avg, dN;
	int i, j, k, N;
	char volfile[74], slicefile[74];
	FILE *fpv, *fps;


	printf(" trying memory allocation for data storage ...");
	volume = f3tensor(1,NY,1,NX,1,NZ);
	printf(" succesfull.\n\n");

	sprintf(volfile,"random3d2.vp");
	sprintf(slicefile,"random2d2.vp");

	fpv=fopen(volfile,"rb");

	printf(" read 3-D volume %s ...",volfile);
	fflush(stdout);
	for (k=1;k<=NZ;k++)
		for (i=1;i<=NX;i++)
			for (j=1;j<=NY;j++){
				fread(&vp,sizeof(float), 1, fpv);
				volume[j][i][k]=vp;			
				}

	printf (" done.\n\n\n");
	fclose(fpv);

	printf(" write vertical slice at %d to %s ...",zslice, slicefile);
	fps=fopen(slicefile,"wb");
	sumvp=0.0;	
	for (k=zslice;k<=zslice;k++)
		for (i=1;i<=NX;i++)
			for (j=1;j<=NY;j++){
				vp=volume[j][i][k];
				sumvp+=(double)vp;
				fwrite(&vp,sizeof(float), 1, fps);
				}
	
	printf (" done.\n\n\n");
	fclose(fps);


	/* compute statistics of 2-D fluctuations */


	N=NX*NY;
	dN=(double)N;
	avg=sumvp/dN;
	sumvp=0.0;	
	for (k=zslice;k<=zslice;k++)
		for (i=1;i<=NX;i++)
			for (j=1;j<=NY;j++){
				vp=volume[j][i][k];
				sumvp+=(avg-vp)*(avg-vp);
			}
					
			
			
	printf(" statistics of 2-D velocity fluctuations: \n");
	std=sqrt(sumvp/dN);
	printf(" average= %e \n",avg);
	printf(" standard deviation = %e \n",std);

	printf(" Use: \n ximage n1=%d n2=%d < %s \n to visualize model.\n\n",
	         NY,NX,slicefile); 


	free_f3tensor(volume,1,NY,1,NX,1,NZ);

		return 0;

}



