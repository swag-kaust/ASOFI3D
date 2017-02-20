/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI3D.
 * 
 * SOFI3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI3D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   Generates elastic model properties (vp,vs,density) on the fly
 *   if L>0 damping model (qp and qs) are generated, too
 *
 *   depending on model dimension in vertical direction and local variable "h"
 *   this function can generate a
 *   	-> homogeneneous full space
 *   	-> layer over half space
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void model_visco(float  ***  rho, float ***  pi, float ***  u,
		float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float DT, *FL, TAU, FREF, DY;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern int WRITE_MODELFILES;
	extern char  MFILE[STRING_SIZE];
	//extern FILE *FP;
	extern float TS;

	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rhov, Qp, Qs;
	float *pts=NULL, sumu=0.0, sumpi=0.0;
	float *** pwavemod=NULL, *** swavemod=NULL;
	float *** qpmod=NULL, *** qsmod=NULL;
	float y;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];

	/*-----------------material property definition -------------------------*/

	/* parameters for layer 1 */
	const float vp1=3500.0, vs1=2000.0, rho1=2000.0, h=100000.0, qp1=20.0, qs1=10.0;

	/* parameters for layer 2 */
	const float vp2=5700.0, vs2=3400.0, rho2=2500.0, qp2=200.0, qs2=100.0;
	
	

	/*internal switch for writing all models to file (WRITE_MODELFILES=1)
	 * or just density (WRITE_MODELFILES=0)
	 * BE AWARE that the output of additional models besides density
	 * cause extra but temporal memory allocation of the size of the
	 * local subgrid times the number of models!!!!!!!!!! */


	if (WRITE_MODELFILES==1) {
		pwavemod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		swavemod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		qpmod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		qsmod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	}


	if (L){  /*viscoelastic simulation */

		/* vector for maxwellbodies */
		pts=vector(1,L);
		for (l=1;l<=L;l++) {
			pts[l]=1.0/(2.0*PI*FL[l]);
			eta[l]=DT/pts[l];
		}

		//previously : ws=2.0*PI*FL[1];
		/* in the viscoelastic case : reference frequency where no velocity dispersion occurs.
		 * if FREF is not given in input file, the largest center source frequency FC
		 * as specified in input file is used (not that the relation : FC=1/TS) is used here)*/
		
		if (FREF==0.0) ws=2.0*PI/TS;
		else ws=2.0*PI*FREF;

		//fprintf(FP,"MYID=%d \t\t ws=%5.5f \t pts=%5.5f \t FL=%5.5f \n ",MYID,ws,pts[l],FL[l]);
		
		/* loop over global grid */
		for (k=1;k<=NZG;k++){
			for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){
					/*=========================================================
					 * modify below this point for VISCOELASTIC model definition
					 *=========================================================
					 */
					
					/*note that "y" is used for the vertical coordinate*/
					/* calculate vertical coordinate in m */

					y=(float)j*DY;

					/* two layer case */
					if (y<=h){
						Vp=vp1; Vs=vs1; Rhov=rho1; Qp=qp1; Qs=qs1; }


					else{
						Vp=vp2; Vs=vs2; Rhov=rho2; Qp=qp2; Qs=qs2; }

					/*=========================================================
					 * modify up to this point for VSICOELASTIC model definition
					 *=========================================================
					 */

					sumu=0.0;
					sumpi=0.0;
					for (l=1;l<=L;l++){
						sumu=sumu+((ws*ws*pts[l]*pts[l]*(2/Qs))/(1.0+ws*ws*pts[l]*pts[l]));
						sumpi=sumpi+((ws*ws*pts[l]*pts[l]*(2/Qp))/(1.0+ws*ws*pts[l]*pts[l]));
					}

					muv=Vs*Vs*Rhov/(1.0+sumu);
					piv=Vp*Vp*Rhov/(1.0+sumpi);

					/* only the PE which belongs to the current global gridpoint
							is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) &&
							(POS[2]==((j-1)/NY)) &&
							(POS[3]==((k-1)/NZ))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						kk=k-POS[3]*NZ;

						u[jj][ii][kk]=muv;
						rho[jj][ii][kk]=Rhov;
						pi[jj][ii][kk]=piv;

						if (TAU==0.0){
							/*calculation of taus and taup by read-in q-files*/
							taus[jj][ii][kk]=2/Qs;
							taup[jj][ii][kk]=2/Qp;
						}
						else {
							/*constant q (damping) case:*/
							taus[jj][ii][kk]=TAU;
							taup[jj][ii][kk]=TAU;
						}

						if (WRITE_MODELFILES==1) {
							pwavemod[jj][ii][kk]=Vp;
							swavemod[jj][ii][kk]=Vs;
							qsmod[jj][ii][kk]=2/taus[jj][ii][kk];
							qpmod[jj][ii][kk]=2/taup[jj][ii][kk];
						}

					}
				} 
			}
		}
		
	}

	/* each PE writes his model to disk */
	
	/* all models are written to file */
	if (WRITE_MODELFILES==1) {
		sprintf(modfile,"%s.SOFI3D.pi",MFILE);
		writemod(modfile,pi,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.u",MFILE);
		writemod(modfile,u,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.vp",MFILE);
		writemod(modfile,pwavemod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.vs",MFILE);
		writemod(modfile,swavemod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	if ((L) && (WRITE_MODELFILES==1)) {
		sprintf(modfile,"%s.SOFI3D.qp",MFILE);
		writemod(modfile,qpmod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI3D.qs",MFILE);
		writemod(modfile,qsmod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	/* only density is written to file */
	if (WRITE_MODELFILES==2) {
		sprintf(modfile,"%s.SOFI3D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	free_vector(pts,1,L);
	if (WRITE_MODELFILES==1) {
		free_f3tensor(pwavemod,0,NY+1,0,NX+1,0,NZ+1);
		free_f3tensor(swavemod,0,NY+1,0,NX+1,0,NZ+1);
		if ((L) && (TAU==0.0)){
			free_f3tensor(qpmod,0,NY+1,0,NX+1,0,NZ+1);
			free_f3tensor(qsmod,0,NY+1,0,NX+1,0,NZ+1);
		}
	}
	
}
