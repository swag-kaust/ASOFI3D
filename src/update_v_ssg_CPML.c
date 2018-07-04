#include "fd.h"
#include "data_structures.h"


/**
 * Correct particle velocity using CPML boundary condition.
 * 
 * CPML (convolutional perfectly matched layer) is an absorbing boundary
 * condition.
 * Correction is done for the 4th order spatial FD sheme.
 *
 * Parameters:
 * nx1, nx2, ny1, ny2, nz1, nz2 :
 *     Grid dimensions.
 * nt :
 *     Time step.
 * v :
 *     Velocity field.
 * s :
 *     Stress tensor.
 * rho, rjp, rkp, rip :
 *     Density and its shifts on the staggered grid.
 * srcpos_loc, signals, nsrc, stype :
 *     Source parameters. NOT USED.
 * absorb_coeff :
 *     Absorption coefficient ?. NOT USED.
 * K_x, a_x, b_x, K_x_half, a_x_half, b_x_half :
 *     Parameters of the Perfectly Matched Layer (PML) along x-axis.
 * K_y, a_y, b_y, K_y_half, a_y_half, b_y_half :
 *     Parameters of the Perfectly Matched Layer (PML) along y-axis.
 * K_z, a_z, b_z, K_z_half, a_z_half, b_z_half :
 *     Parameters of the Perfectly Matched Layer (PML) along z-axis.
 * psi_sxx_x, psi_sxy_x, psi_sxz_x, psi_sxy_y, psi_syy_y,
 * psi_syz_y, psi_sxz_z, psi_syz_z, psi_szz_z :
 *     Memory variables.
 *
 * References
 * ----------
 * Komatitsch, D. and Martin, R., 2007
 * An unsplit convolutional perfectly matched layer
 * improved at grazing incidence for the seismic wave equation
 * Geophysics, Vol. 72, No. 5
 * https://doi.org/10.1190/1.2757586
 */
double update_v_CPML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, Velocity *v,
		Tensor3d *s,
		float  ***  rho,  float  *** rjp, float  *** rkp, float  *** rip,
		float **  srcpos_loc, float ** signals, int nsrc, float *** absorb_coeff, int * stype,
		float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		float * K_z, float * a_z, float * b_z, float * K_z_half, float * a_z_half, float * b_z_half,
		float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y, float *** psi_syy_y,
		float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z) {

        float ***vx = v->x;
        float ***vy = v->y;
        float ***vz = v->z;

        float ***sxx = s->xx;
        float ***syy = s->yy;
        float ***szz = s->zz;
        float ***sxy = s->xy;
        float ***syz = s->yz;
        float ***sxz = s->xz;

	extern float DT, DX, DY, DZ;
	double time=0.0, time1=0.0, time2=0.0;
	extern int MYID, LOG, FDCOEFF, FDORDER;
	extern FILE *FP;
	extern int FREE_SURF;
	extern int NPROCX, NPROCY, NPROCZ, POS[4];
	extern int FW, NY, NZ;
	extern int OUTNTIMESTEPINFO;

	int i, j, k, h1;
	float b1=1.0, b2=0.0, dx, dy, dz;
	float sxx_x=0.0, sxy_y=0.0, sxz_z=0.0, syy_y=0.0, sxy_x=0.0, syz_z=0.0;
	float szz_z=0.0, sxz_x=0.0, syz_y=0.0;

	dx=DT/DX;
	dy=DT/DY;
	dz=DT/DZ;



	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) time1=MPI_Wtime();


	if(FDCOEFF==1){ /* Taylor coefficients*/

		switch (FDORDER){
		case 2 :
			b1=1.0; b2=0.0;
			break;
		case 4 :
			b1=9.0/8.0; b2=-1.0/24.0;
			break;
		default:
			b1=9.0/8.0; b2=-1.0/24.0;
			//warning(" Warning in function update_v_ssg_CPML : CPML Update of FDORDER > 4 not implemented yet, coefficients of FDORDER==4 are used instead!");
			//warning is issued in checkfd_ssg, section "ABSORBING BOUNDARY"
			break;
		}
	}


	if(FDCOEFF==2){ /* Holberg coefficients E=0.1 %*/

		switch (FDORDER){
		case 2 :
			b1=1.0; b2=0.0;
			break;
		case 4 :
			b1=1.1382; b2=-0.046414;
			break;
		default:
			b1=1.1382; b2=-0.046414;
			//warning(" Warning in function update_v_ssg_CPML : CPML Update of FDORDER > 4 not implemented yet, coefficients of FDORDER==4 are used instead!");
			//warning is issued in checkfd_ssg, section "ABSORBING BOUNDARY"
			break;
		}
	}



	/* boundaries in x-direction */

	if (POS[1]==0){
#ifdef _OPENACC
#pragma acc parallel 
#pragma acc loop independent
#endif
		for (j=1;j<=NY;j++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
			for (i=1;i<=FW;i++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
				for (k=1;k<=NZ;k++){

					sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
					sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
					sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
					sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
					syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
					syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
					sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
					syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
					szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));


					psi_sxx_x[j][i][k] = b_x_half[i] * psi_sxx_x[j][i][k] + a_x_half[i] * sxx_x;
					sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i][k];
					psi_sxy_x[j][i][k] = b_x[i] * psi_sxy_x[j][i][k] + a_x[i] * sxy_x;
					sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i][k];
					psi_sxz_x[j][i][k] = b_x[i] * psi_sxz_x[j][i][k] + a_x[i] * sxz_x;
					sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i][k];

					if((POS[2]==0 && FREE_SURF==0) && (j<=FW)){

						psi_sxy_y[j][i][k] = b_y[j] * psi_sxy_y[j][i][k] + a_y[j] * sxy_y;
						sxy_y = sxy_y / K_y[j] + psi_sxy_y[j][i][k];
						psi_syy_y[j][i][k] = b_y_half[j] * psi_syy_y[j][i][k] + a_y_half[j] * syy_y;
						syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i][k];
						psi_syz_y[j][i][k] = b_y[j] * psi_syz_y[j][i][k] + a_y[j] * syz_y;
						syz_y = syz_y / K_y[j] + psi_syz_y[j][i][k]; 	 }


					if((POS[2]==NPROCY-1) && (j>=ny2+1)){
						h1 = (j-ny2+FW);
						psi_sxy_y[h1][i][k] = b_y[h1] * psi_sxy_y[h1][i][k] + a_y[h1] * sxy_y;
						sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i][k];
						psi_syy_y[h1][i][k] = b_y_half[h1] * psi_syy_y[h1][i][k] + a_y_half[h1] * syy_y;
						syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i][k];
						psi_syz_y[h1][i][k] = b_y[h1] * psi_syz_y[h1][i][k] + a_y[h1] * syz_y;
						syz_y = syz_y / K_y[h1] + psi_syz_y[h1][i][k];
					}


					if((POS[3]==0) && (k<=FW)){
						psi_sxz_z[j][i][k] = b_z[k] * psi_sxz_z[j][i][k] + a_z[k] * sxz_z;
						sxz_z = sxz_z / K_z[k] + psi_sxz_z[j][i][k];
						psi_syz_z[j][i][k] = b_z[k] * psi_syz_z[j][i][k] + a_z[k] * syz_z;
						syz_z = syz_z / K_z[k] + psi_syz_z[j][i][k];
						psi_szz_z[j][i][k] = b_z_half[k] * psi_szz_z[j][i][k] + a_z_half[k] * szz_z;
						szz_z = szz_z / K_z_half[k] + psi_szz_z[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){
						h1 = (k-nz2+FW);
						psi_sxz_z[j][i][h1] = b_z[h1] * psi_sxz_z[j][i][h1] + a_z[h1] * sxz_z;
						sxz_z = sxz_z / K_z[h1] + psi_sxz_z[j][i][h1];
						psi_syz_z[j][i][h1] = b_z[h1] * psi_syz_z[j][i][h1] + a_z[h1] * syz_z;
						syz_z = syz_z / K_z[h1] + psi_syz_z[j][i][h1];
						psi_szz_z[j][i][h1] = b_z_half[h1] * psi_szz_z[j][i][h1] + a_z_half[h1] * szz_z;
						szz_z = szz_z / K_z_half[h1] + psi_szz_z[j][i][h1];}


					vx[j][i][k]+= (sxx_x + sxy_y + sxz_z)/rip[j][i][k];
					vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
					vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];

				}
			}
		}
	}

	if(POS[1]==NPROCX-1){
#ifdef _OPENACC
#pragma acc parallel 
#pragma acc loop independent
#endif
		for (j=1;j<=NY;j++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
			for (i=nx2+1;i<=nx2+FW;i++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
				for (k=1;k<=NZ;k++){

					sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
					sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
					sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
					sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
					syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
					syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
					sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
					syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
					szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));

					h1 = (i-nx2+FW);

					psi_sxx_x[j][h1][k] = b_x_half[h1] * psi_sxx_x[j][h1][k] + a_x_half[h1] * sxx_x;
					sxx_x = sxx_x / K_x_half[h1] + psi_sxx_x[j][h1][k];
					psi_sxy_x[j][h1][k] = b_x[h1] * psi_sxy_x[j][h1][k] + a_x[h1] * sxy_x;
					sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1][k];
					psi_sxz_x[j][h1][k] = b_x[h1] * psi_sxz_x[j][h1][k] + a_x[h1] * sxz_x;
					sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1][k];



					if((POS[2]==0 && FREE_SURF==0) && (j<=FW)){
						psi_sxy_y[j][i][k] = b_y[j] * psi_sxy_y[j][i][k] + a_y[j] * sxy_y;
						sxy_y = sxy_y / K_y[j] + psi_sxy_y[j][i][k];
						psi_syy_y[j][i][k] = b_y_half[j] * psi_syy_y[j][i][k] + a_y_half[j] * syy_y;
						syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i][k];
						psi_syz_y[j][i][k] = b_y[j] * psi_syz_y[j][i][k] + a_y[j] * syz_y;
						syz_y = syz_y / K_y[j] + psi_syz_y[j][i][k]; 	 }


					if((POS[2]==NPROCY-1) && (j>=ny2+1)){
						h1 = (j-ny2+FW);
						psi_sxy_y[h1][i][k] = b_y[h1] * psi_sxy_y[h1][i][k] + a_y[h1] * sxy_y;
						sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i][k];
						psi_syy_y[h1][i][k] = b_y_half[h1] * psi_syy_y[h1][i][k] + a_y_half[h1] * syy_y;
						syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i][k];
						psi_syz_y[h1][i][k] = b_y[h1] * psi_syz_y[h1][i][k] + a_y[h1] * syz_y;
						syz_y = syz_y / K_y[h1] + psi_syz_y[h1][i][k];
					}


					if((POS[3]==0) && (k<=FW)){
						psi_sxz_z[j][i][k] = b_z[k] * psi_sxz_z[j][i][k] + a_z[k] * sxz_z;
						sxz_z = sxz_z / K_z[k] + psi_sxz_z[j][i][k];
						psi_syz_z[j][i][k] = b_z[k] * psi_syz_z[j][i][k] + a_z[k] * syz_z;
						syz_z = syz_z / K_z[k] + psi_syz_z[j][i][k];
						psi_szz_z[j][i][k] = b_z_half[k] * psi_szz_z[j][i][k] + a_z_half[k] * szz_z;
						szz_z = szz_z / K_z_half[k] + psi_szz_z[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){
						h1 = (k-nz2+FW);
						psi_sxz_z[j][i][h1] = b_z[h1] * psi_sxz_z[j][i][h1] + a_z[h1] * sxz_z;
						sxz_z = sxz_z / K_z[h1] + psi_sxz_z[j][i][h1];
						psi_syz_z[j][i][h1] = b_z[h1] * psi_syz_z[j][i][h1] + a_z[h1] * syz_z;
						syz_z = syz_z / K_z[h1] + psi_syz_z[j][i][h1];
						psi_szz_z[j][i][h1] = b_z_half[h1] * psi_szz_z[j][i][h1] + a_z_half[h1] * szz_z;
						szz_z = szz_z / K_z_half[h1] + psi_szz_z[j][i][h1];}

					vx[j][i][k]+= (sxx_x + sxy_y + sxz_z)/rip[j][i][k];
					vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
					vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];

				}
			}
		}
	}

	if((POS[2]==0 && FREE_SURF==0)){
#ifdef _OPENACC
#pragma acc parallel 
#pragma acc loop independent
#endif
		for (j=1;j<=FW;j++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
			for (i=nx1;i<=nx2;i++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
				for (k=1;k<=NZ;k++){

					sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
					sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
					sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
					sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
					syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
					syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
					sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
					syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
					szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));

					psi_sxy_y[j][i][k] = b_y[j] * psi_sxy_y[j][i][k] + a_y[j] * sxy_y;
					sxy_y = sxy_y / K_y[j] + psi_sxy_y[j][i][k];

					psi_syy_y[j][i][k] = b_y_half[j] * psi_syy_y[j][i][k] + a_y_half[j] * syy_y;
					syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i][k];

					psi_syz_y[j][i][k] = b_y[j] * psi_syz_y[j][i][k] + a_y[j] * syz_y;
					syz_y = syz_y / K_y[j] + psi_syz_y[j][i][k];

					if((POS[3]==0) && (k<=FW)){
						psi_sxz_z[j][i][k] = b_z[k] * psi_sxz_z[j][i][k] + a_z[k] * sxz_z;
						sxz_z = sxz_z / K_z[k] + psi_sxz_z[j][i][k];
						psi_syz_z[j][i][k] = b_z[k] * psi_syz_z[j][i][k] + a_z[k] * syz_z;
						syz_z = syz_z / K_z[k] + psi_syz_z[j][i][k];
						psi_szz_z[j][i][k] = b_z_half[k] * psi_szz_z[j][i][k] + a_z_half[k] * szz_z;
						szz_z = szz_z / K_z_half[k] + psi_szz_z[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){
						h1 = (k-nz2+FW);
						psi_sxz_z[j][i][h1] = b_z[h1] * psi_sxz_z[j][i][h1] + a_z[h1] * sxz_z;
						sxz_z = sxz_z / K_z[h1] + psi_sxz_z[j][i][h1];
						psi_syz_z[j][i][h1] = b_z[h1] * psi_syz_z[j][i][h1] + a_z[h1] * syz_z;
						syz_z = syz_z / K_z[h1] + psi_syz_z[j][i][h1];
						psi_szz_z[j][i][h1] = b_z_half[h1] * psi_szz_z[j][i][h1] + a_z_half[h1] * szz_z;
						szz_z = szz_z / K_z_half[h1] + psi_szz_z[j][i][h1];}


					vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
					vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
					vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];

				}
			}
		}
	}

	if(POS[2]==NPROCY-1){
#ifdef _OPENACC
#pragma acc parallel 
#pragma acc loop independent
#endif
		for (j=ny2+1;j<=ny2+FW;j++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
			for (i=nx1;i<=nx2;i++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
				for (k=1;k<=NZ;k++){

					sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
					sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
					sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
					sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
					syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
					syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
					sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
					syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
					szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));

					h1 = (j-ny2+FW);

					psi_sxy_y[h1][i][k] = b_y[h1] * psi_sxy_y[h1][i][k] + a_y[h1] * sxy_y;
					sxy_y = sxy_y / K_y[h1] + psi_sxy_y[h1][i][k];
					psi_syy_y[h1][i][k] = b_y_half[h1] * psi_syy_y[h1][i][k] + a_y_half[h1] * syy_y;
					syy_y = syy_y / K_y_half[h1] + psi_syy_y[h1][i][k];
					psi_syz_y[h1][i][k] = b_y[h1] * psi_syz_y[h1][i][k] + a_y[h1] * syz_y;
					syz_y = syz_y / K_y[h1] + psi_syz_y[h1][i][k];

					if((POS[3]==0) && (k<=FW)){
						psi_sxz_z[j][i][k] = b_z[k] * psi_sxz_z[j][i][k] + a_z[k] * sxz_z;
						sxz_z = sxz_z / K_z[k] + psi_sxz_z[j][i][k];
						psi_syz_z[j][i][k] = b_z[k] * psi_syz_z[j][i][k] + a_z[k] * syz_z;
						syz_z = syz_z / K_z[k] + psi_syz_z[j][i][k];
						psi_szz_z[j][i][k] = b_z_half[k] * psi_szz_z[j][i][k] + a_z_half[k] * szz_z;
						szz_z = szz_z / K_z_half[k] + psi_szz_z[j][i][k];}


					if((POS[3]==NPROCZ-1) && (k>=nz2+1)){
						h1 = (k-nz2+FW);
						psi_sxz_z[j][i][h1] = b_z[h1] * psi_sxz_z[j][i][h1] + a_z[h1] * sxz_z;
						sxz_z = sxz_z / K_z[h1] + psi_sxz_z[j][i][h1];
						psi_syz_z[j][i][h1] = b_z[h1] * psi_syz_z[j][i][h1] + a_z[h1] * syz_z;
						syz_z = syz_z / K_z[h1] + psi_syz_z[j][i][h1];
						psi_szz_z[j][i][h1] = b_z_half[h1] * psi_szz_z[j][i][h1] + a_z_half[h1] * szz_z;
						szz_z = szz_z / K_z_half[h1] + psi_szz_z[j][i][h1];}

					vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
					vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
					vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];

				} 
			}
		}
	}

	/* boundaries in z-direction */

	if(POS[3]==0){
#ifdef _OPENACC
#pragma acc parallel 
#pragma acc loop independent
#endif
		for (j=ny1;j<=ny2;j++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
			for (i=nx1;i<=nx2;i++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
				for (k=1;k<=FW;k++){

					sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
					sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
					sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
					sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
					syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
					syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
					sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
					syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
					szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));


					psi_sxz_z[j][i][k] = b_z[k] * psi_sxz_z[j][i][k] + a_z[k] * sxz_z;
					sxz_z = sxz_z / K_z[k] + psi_sxz_z[j][i][k];
					psi_syz_z[j][i][k] = b_z[k] * psi_syz_z[j][i][k] + a_z[k] * syz_z;
					syz_z = syz_z / K_z[k] + psi_syz_z[j][i][k];
					psi_szz_z[j][i][k] = b_z_half[k] * psi_szz_z[j][i][k] + a_z_half[k] * szz_z;
					szz_z = szz_z / K_z_half[k] + psi_szz_z[j][i][k];

					vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
					vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
					vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];

				}
			}
		}
	}

	if(POS[3]==NPROCZ-1){
#ifdef _OPENACC
#pragma acc parallel 
#pragma acc loop independent
#endif
		for (j=ny1;j<=ny2;j++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
			for (i=nx1;i<=nx2;i++){
#ifdef _OPENACC
#pragma acc loop independent
#endif
				for (k=nz2+1;k<=nz2+FW;k++){
					sxx_x = dx*(b1*(sxx[j][i+1][k]-sxx[j][i][k])+b2*(sxx[j][i+2][k]-sxx[j][i-1][k]));
					sxy_x = dx*(b1*(sxy[j][i][k]-sxy[j][i-1][k])+b2*(sxy[j][i+1][k]-sxy[j][i-2][k]));
					sxz_x = dx*(b1*(sxz[j][i][k]-sxz[j][i-1][k])+b2*(sxz[j][i+1][k]-sxz[j][i-2][k]));
					sxy_y = dy*(b1*(sxy[j][i][k]-sxy[j-1][i][k])+b2*(sxy[j+1][i][k]-sxy[j-2][i][k]));
					syy_y = dy*(b1*(syy[j+1][i][k]-syy[j][i][k])+b2*(syy[j+2][i][k]-syy[j-1][i][k]));
					syz_y = dy*(b1*(syz[j][i][k]-syz[j-1][i][k])+b2*(syz[j+1][i][k]-syz[j-2][i][k]));
					sxz_z = dz*(b1*(sxz[j][i][k]-sxz[j][i][k-1])+b2*(sxz[j][i][k+1]-sxz[j][i][k-2]));
					syz_z = dz*(b1*(syz[j][i][k]-syz[j][i][k-1])+b2*(syz[j][i][k+1]-syz[j][i][k-2]));
					szz_z = dz*(b1*(szz[j][i][k+1]-szz[j][i][k])+b2*(szz[j][i][k+2]-szz[j][i][k-1]));

					h1 = (k-nz2+FW);

					psi_sxz_z[j][i][h1] = b_z[h1] * psi_sxz_z[j][i][h1] + a_z[h1] * sxz_z;
					sxz_z = sxz_z / K_z[h1] + psi_sxz_z[j][i][h1];
					psi_syz_z[j][i][h1] = b_z[h1] * psi_syz_z[j][i][h1] + a_z[h1] * syz_z;
					syz_z = syz_z / K_z[h1] + psi_syz_z[j][i][h1];
					psi_szz_z[j][i][h1] = b_z_half[h1] * psi_szz_z[j][i][h1] + a_z_half[h1] * szz_z;
					szz_z = szz_z / K_z_half[h1] + psi_szz_z[j][i][h1];

					vx[j][i][k]+= (sxx_x + sxy_y +sxz_z)/rip[j][i][k];
					vy[j][i][k]+= (syy_y + sxy_x + syz_z)/rjp[j][i][k];
					vz[j][i][k]+= (szz_z + sxz_x + syz_y)/rkp[j][i][k];
				}
			}
		}
	}


	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
			time2=MPI_Wtime();
			time=time2-time1;
			fprintf(FP," Real time for CPML particle velocity update: \t %4.2f s.\n",time);
		}
	return time;

}
