/*------------------------------------------------------------------------
 *   Generates elastic model properties (vp,vs,density, Cij coefficients) on the fly
 *
 *   depending on model dimension in vertical direction and local variable "h"
 *   this function can generate a
 *   	-> homogeneneous full space
 *   	-> layer over half space
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void model_elastic(float  ***  rho, float ***  pi, float ***  u,
        float *** C11, float *** C12, float *** C13, float *** C22, float *** C23, float *** C33,
        float *** C44, float *** C55, float *** C66,
        float ***  taus, float ***  taup, float *  eta){

    /*--------------------------------------------------------------------------*/
    /* extern variables */
    extern float DY;
    extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
    extern char  MFILE[STRING_SIZE];
    extern int WRITE_MODELFILES;
    extern FILE *FP;

    /* local variables */
    float muv, piv;
    float Vpv, Vsv, Rho, Poi, Epsx, Delx, Gamx;
    float *** vpv=NULL, *** vsv=NULL, *** epsx=NULL, *** epsy=NULL, *** gamx=NULL;
    float *** delx=NULL, *** dely=NULL, *** delxy=NULL, *** gamy=NULL;
    float y;
    int i, j, k, ii, jj, kk;
    char modfile[STRING_SIZE];

    /*-----------------material property definition -------------------------*/

    /* x=1, y=2 in Tsvankin [1997] (e.g.) epsx=epsion1 & epsy=epsilon2 */

    /* parameters for layer 1 */
   /* const float vpv1=2326.0, poi1=0.25,
     epsx1=0.135, epsy1=-0.082, delx1=-0.166,
     dely1=-0.24, delxy1=-0.089,
     gamx1=0.438, gamy1=0.25,
     rho1=2000.0, h=100000.0; */
   const float vpv1=3500.0, poi1=0.25,
           epsx1=0.0, epsy1=0.0,
           delx1=0.0, dely1=0.0,
           delxy1=0,
           gamx1=0.0, gamy1=0.0, rho1=2000.0, h=750000.0;
    /* parameters for layer 2 */
    //const float vp2=5700.0, vs2=3400.0, rho2=2500.0;
    const float vpv2=3000.0, poi2=0.25, epsx2=0., epsy2=0.0, delx2=0., dely2=0., delxy2=0.,
            gamx2=0., gamy2=0., rho2=1560.0;


    if (WRITE_MODELFILES==1) {
        vpv  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        vsv  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        epsx  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        epsy  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        delx  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        dely  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        delxy  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        gamx  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        gamy  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
    }

    /*elastic simulation */
    if (L==0) {
        /* loop over global grid */
        fprintf(FP,"In HH elastic MYID=%d, POS[1]=%d, POS[2]=%d,POS[3]=%d \n\n",MYID,POS[1],POS[2],POS[3]);
        for (k=1;k<=NZG;k++){
            for (i=1;i<=NXG;i++){
                for (j=1;j<=NYG;j++){

                    /*=========================================================
                     * modify below this point for ELASTIC model definition
                     *=========================================================
                     */
                    /*note that "y" is used for the vertical coordinate*/
                    /* calculate vertical coordinate in m */

                    y=(float)j*DY;



                    /* two layer case */
                    if (y<=h){
                        Vpv=vpv1;
                        Poi=poi1;
                        Rho=rho1;
                    } else{
                        Vpv=vpv2;
                        Poi=poi2;
                        Rho=rho2;
                    }
                    // perturbation in the midle of the model
                    if (((i-(NZG/2))*(i-(NZG/2)) + (j-(NZG/2))*(j-(NZG/2)) + (k-(NZG/2))*(k-(NZG/2))) <= 5*5){
                        //Vpv -= Vpv * 0.1;
                    }


                    /*=========================================================
                     * modify up to this point for ELASTIC model definition
                     *=========================================================
                     */

                    Vsv=Vpv*sqrt((1-2*Poi)/(2-2*Poi));
                    muv=Vsv*Vsv*Rho;
                    piv=Vpv*Vpv*Rho;

                    /* only the PE which belongs to the current global gridpoint
                     * is saving model parameters in his local arrays */

                    if ((POS[1]==((i-1)/NX)) &&
                            (POS[2]==((j-1)/NY)) &&
                            (POS[3]==((k-1)/NZ))){
                        ii=i-POS[1]*NX;
                        jj=j-POS[2]*NY;
                        kk=k-POS[3]*NZ;

                        u[jj][ii][kk]=muv;
                        pi[jj][ii][kk]=piv;
                        //VTI
                        C11[jj][ii][kk] = (1+2*Epsx)*Rho*Vpv*Vpv;
                        C22[jj][ii][kk] = C11[jj][ii][kk];
                        C33[jj][ii][kk] = Rho*Vpv*Vpv;
                        C66[jj][ii][kk] = Rho*Vsv*Vsv;
                        C12[jj][ii][kk] = C11[jj][ii][kk] - 2*C66[jj][ii][kk];
                        C13[jj][ii][kk]=Rho*sqrt((Vpv*Vpv-Vsv*Vsv)*((1+2*Delx)*Vpv*Vpv-Vsv*Vsv))-Rho*Vsv*Vsv;
                        C23[jj][ii][kk]=C13[jj][ii][kk];
                        C44[jj][ii][kk]=Rho*Vsv*Vsv/(1+2*Gamx);
                        C55[jj][ii][kk]=Rho*Vsv*Vsv/(1+2*Gamx);
			rho[jj][ii][kk]=Rho;

                        //ORTHO
                        // C_33 = Rho*Vpv*Vpv;
                        // C_55 = Rho*Vsv*Vsv;
                        // C_66 = (1+2*gamma_1)*C_55;
                        // C_11 = (1+2*eps_2)*C_33;
                        // C_44 = C_66/(1+2*gamma_2);

                        // C33[jj][ii][kk]=C_33;
                        // C55[jj][ii][kk]=C_55;
                        // C11[jj][ii][kk]=C_11;
                        // C22[jj][ii][kk]=(1+2*eps_1)*C_33;
                        // C66[jj][ii][kk]=C_66;
                        // C44[jj][ii][kk]=C_44;
                        // C13[jj][ii][kk]=-C_55+sqrt(2*C_33*(C_33-C_55)*delta_2 + (C_33-C_55)*(C_33-C_55));
                        // C12[jj][ii][kk]=-C_66+sqrt(2*delta_3*C_11*(C_11-C_66) + (C_11-C_66)*(C_11-C_66));
                        // C23[jj][ii][kk]=-C_44+sqrt(2*delta_1*C_33*(C_33-C_44) + (C_33-C_44)*(C_33-C_44));
                        // rho[jj][ii][kk]=Rho;

                        if (WRITE_MODELFILES==1) {
                            vpv[jj][ii][kk]=Vpv;
                            vsv[jj][ii][kk]=Vsv;
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
        writemod(modfile,vpv,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(modfile,3);

        sprintf(modfile,"%s.SOFI3D.vs",MFILE);
        writemod(modfile,vsv,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(modfile,3);

        sprintf(modfile,"%s.SOFI3D.rho",MFILE);
        writemod(modfile,rho,3);
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


    if (WRITE_MODELFILES==1) {
        free_f3tensor(vpv,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(vsv,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(epsx,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(epsy,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(delx,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(dely,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(delxy,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(gamx,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(gamy,0,NY+1,0,NX+1,0,NZ+1);
    }
}
