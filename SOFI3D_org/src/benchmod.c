/* -------------------------------------------------------------
 *   Model for overnight built, fullspace, elastic case
 *
 *   ------------------------------------------------------------- */

    #include "fd.h"

    void model_elastic(float  ***  rho, float ***  pi, float ***  u, 
    float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float DT, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];

	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho, Qp, Qs;
	float *pts=NULL, sumu=0.0, sumpi=0.0;
	float *** pwavemod, *** swavemod;
	float *** qpmod=NULL, *** qsmod=NULL;
	
	int i, j, k, l, ii, jj, kk;
	int writeallmodels=0;
	char filename[STRING_SIZE];
	    
	    /*-----------------material property definition -------------------------*/

	    /* parameters for underground */
	    const float vp1=3500.0, vs1=2000.0, rho1=2000.0,  qp1=20.0, qs1=10.0;
	    
	    /*-----------------------------------------------------------------------*/
	    /*internal switch for writing all models to file (writeallmodels=1)
	    * or just density (writeallmodels=0)
	    * BE AWARE that the output of additional models besides density
	    * cause extra but temporal memory allocation of the size of the
	    * local subgrid times the number of models!!!!!!!!!! */
	    writeallmodels=0;

	    if (writeallmodels) {
		    pwavemod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		    swavemod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		    qpmod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		    qsmod  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	    }

	    /*elastic simulation */
	    if (L==0) {
		    /* loop over global grid */
		    for (k=1;k<=NZG;k++){
			    for (i=1;i<=NXG;i++){
				    for (j=1;j<=NYG;j++){
					    /*for usability reasons, "z" - as commonly used - denotes the depth (vertical direction),
					      however, internally "y" is used for the vertical coordinate,
					      we simply switch the "y" and "z" coordinate as read in the input file,
					      therefore use the variable "j" to, e.g., calculate/address the vertical coordinate:
					      x=(float)i*DX;
					      y=(float)k*DZ;
					      z=(float)j*DY;*/
					    /* calculate vertical coordinate in m */

					    /*=========================================================
					    * modify below to this point for ELASTIC model definition
					    * visco-elastic model is generated BELOW
					    *=========================================================
					    */

					   
					    Vp=vp1; Vs=vs1; Rho=rho1;


					    /*=========================================================
					    * modify up to this point for ELASTIC model definition
					    * visco-elastic model is generated BELOW
					    *=========================================================
					    */

					    muv=Vs*Vs*Rho;
					    piv=Vp*Vp*Rho;

					    /* only the PE which belongs to the current global gridpoint
							    is saving model parameters in his local arrays */
					    if ((POS[1]==((i-1)/NX)) &&
							    (POS[2]==((j-1)/NY)) &&
							    (POS[3]==((k-1)/NZ))){
						    ii=i-POS[1]*NX;
						    jj=j-POS[2]*NY;
						    kk=k-POS[3]*NZ;

						    u[jj][ii][kk]=muv;
						    rho[jj][ii][kk]=Rho;
						    pi[jj][ii][kk]=piv;

						    if (writeallmodels) {
							    pwavemod[jj][ii][kk]=Vp;
							    swavemod[jj][ii][kk]=Vs;
						    }
					    }
				    }
			    }
		    }
	    }

	    if (L){  /*viscoelastic simulation */

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

					    /*for usability reasons, "z" - as commonly used - denotes the depth (vertical direction),
					      however, internally "y" is used for the vertical coordinate,
					      we simply switch the "y" and "z" coordinate as read in the input file,
					      therefore use the variable "j" to, e.g., calculate/address the vertical coordinate:
					      x=(float)i*DX;
					      y=(float)k*DZ;
					      z=(float)j*DY;*/
					    /* calculate vertical coordinate in m */

					    /*=========================================================
					    * modify below to this point for VISOCELASTIC model definition
					    * elastic model is generated above
					    *=========================================================
					    */

					    Vp=vp1; Vs=vs1; Rho=rho1; Qp=qp1; Qs=qs1;
					    

					    /*=========================================================
					    * modify up to this point for ELASTIC model definition
					    * elastic model is generated above
					    *=========================================================
					    */

					    sumu=0.0;
					    sumpi=0.0;
					    for (l=1;l<=L;l++){
						    sumu=sumu+((ws*ws*pts[l]*pts[l]*(2/Qs))/(1.0+ws*ws*pts[l]*pts[l]));
						    sumpi=sumpi+((ws*ws*pts[l]*pts[l]*(2/Qp))/(1.0+ws*ws*pts[l]*pts[l]));
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

						    u[jj][ii][kk]=muv;
						    rho[jj][ii][kk]=Rho;
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

						    if (writeallmodels) {
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

	    if (writeallmodels) {
		    sprintf(filename,"%s.SOFI3D.pi",MFILE);
		    writemod(filename,pi,3);
		    MPI_Barrier(MPI_COMM_WORLD);
		    if (MYID==0) mergemod(filename,3);

		    sprintf(filename,"%s.SOFI3D.u",MFILE);
		    writemod(filename,u,3);
		    MPI_Barrier(MPI_COMM_WORLD);
		    if (MYID==0) mergemod(filename,3);

		    sprintf(filename,"%s.SOFI3D.vp",MFILE);
		    writemod(filename,pwavemod,3);
		    MPI_Barrier(MPI_COMM_WORLD);
		    if (MYID==0) mergemod(filename,3);

		    sprintf(filename,"%s.SOFI3D.vs",MFILE);
		    writemod(filename,swavemod,3);
		    MPI_Barrier(MPI_COMM_WORLD);
		    if (MYID==0) mergemod(filename,3);
	    }

	    sprintf(filename,"%s.SOFI3D.rho",MFILE);
	    writemod(filename,rho,3);
	    MPI_Barrier(MPI_COMM_WORLD);
	    if (MYID==0) mergemod(filename,3);

	    if ((L) && (writeallmodels)) {
		    sprintf(filename,"%s.SOFI3D.qp",MFILE);
		    writemod(filename,qpmod,3);
		    MPI_Barrier(MPI_COMM_WORLD);
		    if (MYID==0) mergemod(filename,3);

		    sprintf(filename,"%s.SOFI3D.qs",MFILE);
		    writemod(filename,qsmod,3);
		    MPI_Barrier(MPI_COMM_WORLD);
		    if (MYID==0) mergemod(filename,3);
	    }


	    free_vector(pts,1,L);
	    if (writeallmodels) {
		    free_f3tensor(pwavemod,0,NY+1,0,NX+1,0,NZ+1);
		    free_f3tensor(swavemod,0,NY+1,0,NX+1,0,NZ+1);
		    if ((L) && (TAU==0.0)){
			    free_f3tensor(qpmod,0,NY+1,0,NX+1,0,NZ+1);
			    free_f3tensor(qsmod,0,NY+1,0,NX+1,0,NZ+1);
		    }
	    }

    }
