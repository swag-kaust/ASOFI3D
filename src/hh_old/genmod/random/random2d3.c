/*
 *  Generate 3-D random medium
 *   last update 22.01.00, T. Bohlen
 */

/* files to include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "fftn.c"
#include "bessik.c"
#include "beschb.c"
#include "chebev.c"
#include "gammln.c"
#include "nrutil.c"

#ifndef M_PI
# define M_PI	3.14159265358979323846264338327950288
#endif

double drand48(void);
void srand48(long seedval);

int main(int argc, char **argv){


	const int NX=240, NY=600;
	const float DH=2.0;
	const float vm=3000.0, sigma_proz=5.0;
	const float vcut=5.0*sigma_proz*vm/100.0;
	float nu=0.15;  /* Hurst coefficient*/
	float cl=45.0; /* correlation length */

	/* local variables */
	float *A1r, *A1i;
	float *x, *y, ri, rk, rip, rkp;
	float a1r, a1i, rand_max, r, a, phi, sigma, mem, max=0.0, min=0.0;
	float vp, vpupl, vplol, vpmax=vm, vpmin=vm, g;
	double sumvp=0.0, avg, dN, std;
	int i, j, N;
	int dims [2];
	char modfile[74];
	FILE *fp;

	rand_max=(float)RAND_MAX;
	sigma=sigma_proz*vm/100.0;

	dims[0]=NY;
	dims[1]=NX;

	vpupl=vm+vcut;
	vplol=vm-vcut;
	printf(" generating 2-D random media model...\n");
	printf(" size: NX=%d \t NY=%d  gridpoints\n",NX,NY);
	printf(" average velocity: %5.2f m/s\n",vm);
	printf(" standard deviation: %5.2f m/s\n",sigma);
	printf(" cutting velocities lower than %5.2f and higher than %5.2f m/s\n",vplol,vpupl);
	printf(" Correlation length: a=%5.3f meter \n",cl);
	printf(" Hurst coefficient: nu=%5.3f \n",nu);
	printf(" grid spacing: %5.2f meter \n",DH);

	printf("\n\n\n");
	printf(" starting ...\n\n");
	N=(NX)*(NY);
	mem=4.0*(2.0*(N+1)+NX+1+NY+1)*pow(2,-20);
	printf(" trying memory allocation for data storage (%5.2f  mega bytes) ...",mem);
	A1r = (float *) calloc (N+1,sizeof(float));
	A1i = (float *) calloc (N+1,sizeof(float));
	x = (float *) calloc ((NX+1),sizeof(float));
	y = (float *) calloc ((NY+1),sizeof(float));

	if (A1r == NULL || A1i == NULL || x == NULL || x == NULL) {
		fprintf (stderr, "Unable to allocate memory for data storage.\n");
		return 1;
	}

	printf(" succesfull.\n\n");

	/* calculate autocorrelation function (ACF) for exponential medium */

	/*vectors for spatial directions */
	for (i=1;i<=NX;i++) x[i]=(-(NX/2)+i)*DH;
	for (j=1;j<=NY;j++) y[j]=(-(NY/2)+j)*DH;

	/* ACF */
	printf(" calculating von Karman autocorrelation function for \n");
	printf(" correlation length=%5.3f, nu=%5.3f. \n\n",cl,nu);
	fflush(stdout);

	g=pow(2.0,1-nu)/exp(gammln(nu));
	printf("gamma(%f)=%f\n",nu,exp(gammln(nu)));

				
	for (i=1;i<=NX;i++)
		for (j=1;j<=NY;j++){
			r=sqrt(x[i]*x[i]+y[j]*y[j])/cl;
			if (r==0.0) r=1e-10;
			/* exponential ACF */
			/*a=exp(-r);*/

			/* gaussian ACF */
			/*a=exp(-r*r);*/

			/* von-Karman ACF */
			bessik(r,nu,&ri,&rk,&rip,&rkp);
			a=g*pow(r,nu)*rk;

			/*printf("i=%d \t j=%d \t r=%f \t ACF=%f\n",i,j,r,a);*/

			A1r[(i-1)*NY+(j)]=a;
			A1i[(i-1)*NY+(j)]=0.0;
			if (a>max) max=a;			
			if (a<min) min=a;			
				
				
		}

	printf (" done.\n\n\n");
	printf(" maximum of ACF= %5.2f \t minimum of ACF= %5.2f \n\n",max,min);

	/* 3D forward fft */
	printf (" starting forward 2-D fft...");
	fflush(stdout);
	fftnf(2, dims, &A1r[1], &A1i[1], 1, 0.0);
	printf (" done.\n\n\n");

	printf (" muliply amplitude spectrum with random phase fluctuations...");
	fflush(stdout);

	srand48(N);
		for (i=1;i<=N;i++){
		a1r= A1r[i];
		a1i= A1i[i];
		a=sqrt(sqrt(a1r*a1r+a1i*a1i));
		r=(float)drand48();
		/* printf("r(%d)=%f\n",i,r);*/
		phi=-M_PI*(2.0*r-1.0);
		A1r[i]=a*cos(phi);
		A1i[i]=a*sin(phi);
	}

	printf ("done.\n\n\n");

	/* 2D inverse fft */
	printf (" starting inverse 2-D fft...");
	fflush(stdout);
	fftnf(2, dims, &A1r[1], &A1i[1], -1,(double) N);
	printf (" done.\n\n\n");


	printf(" scaling the fluctuations to the desired standard deviation \n");
	printf(" (sigma = %5.2f percent (%5.2f m/s)) ... ",sigma_proz,sigma);
	fflush(stdout);

	for (i=1;i<=N;i++){
		/*A1r[i]=vm+(A1r[i]*sqrt(2.0)*sigma*sqrt((double)(N)));*/
		A1r[i]=vm+(A1r[i]*sigma*sqrt(2)*sqrt((double)(NX*NY)));
	}

	printf (" done.\n\n\n");



	/* compute statistics of fluctuations */


	sumvp=0.0;	
	for (i=1;i<=N;i++){				
		vp=A1r[i];
		if (vp>vpupl) vp=vpupl;
		if (vp<vplol) vp=vplol;
		if (vp>vpmax) vpmax=vp;
		if (vp<vpmin) vpmin=vp;
		sumvp+=(double)vp;
	}

	dN=(double)N;
	avg=sumvp/dN;

	sumvp=0.0;	
	for (i=1;i<=N;i++){	
		vp=A1r[i];
		sumvp+=(avg-vp)*(avg-vp);
	}
					
			
			
	printf(" statistics of velocity fluctuations after inverse fft: \n");
	std=sqrt(sumvp/dN);
	printf(" average= %e \n",avg);
	printf(" standard deviation = %e \n",std);
	printf(" maximum velocity = %f \n", vpmax);
	printf(" minmum velocity = %f \n",vpmin);


	sprintf(modfile,"random2d3.vp");
	printf(" write velocity fluctuations to ");
	printf(modfile); 
	printf("\n");
	fflush(stdout);
	fp=fopen(modfile,"wb");

	for (i=1;i<=NX;i++)
		for (j=1;j<=NY;j++){
			vp=A1r[(i-1)*NY+(j)];
			fwrite(&vp,sizeof(float), 1, fp);
		}


	fclose(fp);
	printf (" done.\n\n\n");
	
	printf(" Use: \n ximage n1=%d n2=%d < %s \n to visualize model.\n\n",
	         NY,NX,modfile); 

	fft_free();
	free(A1i);
	free(A1r);

	free(x);
	free(y);

	return 0;

}



