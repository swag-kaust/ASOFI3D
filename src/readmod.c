/*------------------------------------------------------------------------
 *   Read elastic model properties (vp,vs,density) from files
 *   if L>0 damping model (qp and qs) are from file read, too
 *
 *  ----------------------------------------------------------------------*/
#include <stdbool.h>

#include "fd.h"


void readmod(float ***rho, float ***pi, float ***u,
    float ***C11, float ***C12, float ***C13,
    float ***C22, float ***C23, float ***C33,
    float ***C44, float ***C55, float ***C66,
    float ***taus, float ***taup, float *eta)
{
    // Global variables.
    extern float DT, *FL, TAU, TS, FREF;
    extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
    extern int WRITE_MODELFILES;
    extern char MFILE[STRING_SIZE];
    extern FILE *FP;


    // Local variables.
    float Rho = 0.0, Vp = 0.0, Vs = 0.0, Qp = 0.0, Qs = 0.0;
    float muv = 0.0, piv = 0.0;

    // Stiffness parameters of the model at a given point (i, j, k).
    float C_11, C_22, C_33, C_44, C_55, C_66, C_12, C_13, C_23;

    float *pts = NULL, sumu = 0.0, sumpi = 0.0, ws = 0.0;
    float ***pwavemod = NULL, ***swavemod = NULL;
    float ***qpmod = NULL, ***qsmod = NULL;
    int i, j, k, ii, jj, kk;

    FILE *fp_vs, *fp_vp, *fp_rho, *fp_qp = NULL, *fp_qs = NULL;
    FILE *fp_C11, *fp_C22, *fp_C33, *fp_C44, *fp_C55, *fp_C66;
    FILE *fp_C12, *fp_C13, *fp_C23;

    char filename[STRING_SIZE];

    char fname_C11[STRING_SIZE];
    char fname_C22[STRING_SIZE];
    char fname_C33[STRING_SIZE];
    char fname_C44[STRING_SIZE];
    char fname_C55[STRING_SIZE];
    char fname_C66[STRING_SIZE];
    char fname_C12[STRING_SIZE];
    char fname_C13[STRING_SIZE];
    char fname_C23[STRING_SIZE];

    char fname_epsx[STRING_SIZE];
    char fname_epsy[STRING_SIZE];
    char fname_delx[STRING_SIZE];
    char fname_dely[STRING_SIZE];
    char fname_delxy[STRING_SIZE];
    char fname_gamx[STRING_SIZE];
    char fname_gamy[STRING_SIZE];
    FILE *fp_epsx;
    FILE *fp_epsy;
    FILE *fp_delx;
    FILE *fp_dely;
    FILE *fp_delxy;
    FILE *fp_gamx;
    FILE *fp_gamy;

    bool flag_velocity_files_avail = false;
    bool flag_cij_files_avail = false;


    /* choose data format: ASCII: format=2, BINARY: format=3*/
    const int format = 3;

    /*internal switch for writing all models to file (WRITE_MODELFILES=1)
	 * or just density (WRITE_MODELFILES=0)
	 * BE AWARE that the output of additional models besides density
	 * cause extra but temporal memory allocation of the size of the
	 * local subgrid times the number of models!!!!!!!!!! */


    if (WRITE_MODELFILES) {
        pwavemod = f3tensor(0, NY + 1, 0, NX + 1, 0, NZ + 1);
        swavemod = f3tensor(0, NY + 1, 0, NX + 1, 0, NZ + 1);
        qpmod = f3tensor(0, NY + 1, 0, NX + 1, 0, NZ + 1);
        qsmod = f3tensor(0, NY + 1, 0, NX + 1, 0, NZ + 1);
    }


    fprintf(FP, "\nReading model information from model files...\n");

    // Read density model.
    char fname_rho[STRING_SIZE];
    fprintf(FP, "\t Density:\n\t %s.rho\n\n", MFILE);
    sprintf(fname_rho, "%s.rho", MFILE);
    fp_rho = fopen(fname_rho, "r");
    if (fp_rho == NULL) err(" Could not open model file for densities ! ");

    // Read P- and S-velocity models.
    {
        char fname_vp[STRING_SIZE];
        char fname_vs[STRING_SIZE];
        fprintf(FP, "\t P-wave velocities:\n\t %s.vp\n\n", MFILE);
        sprintf(fname_vp, "%s.vp", MFILE);
        fp_vp = fopen(fname_vp, "r");

        fprintf(FP, "\t Shear wave velocities:\n\t %s.vs\n\n", MFILE);
        sprintf(fname_vs, "%s.vs", MFILE);
        fp_vs = fopen(fname_vs, "r");

        // Checking that either model files for Vp and Vs are both present
        // or both absent, otherwise terminate with error.
        if (fp_vp != NULL && fp_vs != NULL) {
            flag_velocity_files_avail = true;
        } else if (fp_vp == NULL && fp_vs != NULL) {
            err("Could not open model file for P-velocity "
                "but model file for S-velocity is available.");
        } else if (fp_vp != NULL && fp_vs == NULL) {
            err("Could not open model file for S-velocity "
                "but model file for P-velocity is available.");
        } else {
            flag_velocity_files_avail = false;
        }
    }

    // Try opening files for anisotropic parameters.
    {
        fprintf(FP, "\t epsx field:\n\t %s.epsx\n\n", MFILE);
        sprintf(fname_epsx, "%s.epsx", MFILE);
        fp_epsx = fopen(fname_epsx, "r");

        fprintf(FP, "\t epsy field:\n\t %s.epsy\n\n", MFILE);
        sprintf(fname_epsy, "%s.epsy", MFILE);
        fp_epsy = fopen(fname_epsy, "r");

        fprintf(FP, "\t delx field:\n\t %s.delx\n\n", MFILE);
        sprintf(fname_delx, "%s.delx", MFILE);
        fp_delx = fopen(fname_delx, "r");

        fprintf(FP, "\t dely field:\n\t %s.dely\n\n", MFILE);
        sprintf(fname_dely, "%s.dely", MFILE);
        fp_dely = fopen(fname_dely, "r");

        fprintf(FP, "\t delxy field:\n\t %s.delxy\n\n", MFILE);
        sprintf(fname_delxy, "%s.delxy", MFILE);
        fp_delxy = fopen(fname_delxy, "r");

        fprintf(FP, "\t gamx field:\n\t %s.gamx\n\n", MFILE);
        sprintf(fname_gamx, "%s.epsx", MFILE);
        fp_gamx = fopen(fname_gamx, "r");

        fprintf(FP, "\t gamy field:\n\t %s.gamy\n\n", MFILE);
        sprintf(fname_gamy, "%s.gamy", MFILE);
        fp_gamy = fopen(fname_gamy, "r");
    }


    // Try to read model files for anistropic stiffness parameters Cij.
    // Check that either all 9 files for Cij fields are available for reading
    // or all of them are absent.
    // Otherwise, if some Cij files are present and some are absent,
    // terminate the program.
    {
        fprintf(FP, "\t C11:\n\t %s.C11\n\n", MFILE);
        sprintf(fname_C11, "%s.C11", MFILE);
        fp_C11 = fopen(fname_C11, "r");

        fprintf(FP, "\t C22:\n\t %s.C22\n\n", MFILE);
        sprintf(fname_C22, "%s.C22", MFILE);
        fp_C22 = fopen(fname_C22, "r");

        fprintf(FP, "\t C33:\n\t %s.C33\n\n", MFILE);
        sprintf(fname_C33, "%s.C33", MFILE);
        fp_C33 = fopen(fname_C33, "r");

        fprintf(FP, "\t C44:\n\t %s.C44\n\n", MFILE);
        sprintf(fname_C44, "%s.C44", MFILE);
        fp_C44 = fopen(fname_C44, "r");

        fprintf(FP, "\t C55:\n\t %s.C55\n\n", MFILE);
        sprintf(fname_C55, "%s.C55", MFILE);
        fp_C55 = fopen(fname_C55, "r");

        fprintf(FP, "\t C66:\n\t %s.C66\n\n", MFILE);
        sprintf(fname_C66, "%s.C66", MFILE);
        fp_C66 = fopen(fname_C66, "r");

        fprintf(FP, "\t C12:\n\t %s.C12\n\n", MFILE);
        sprintf(fname_C12, "%s.C12", MFILE);
        fp_C12 = fopen(fname_C12, "r");

        fprintf(FP, "\t C13:\n\t %s.C13\n\n", MFILE);
        sprintf(fname_C13, "%s.C13", MFILE);
        fp_C13 = fopen(fname_C13, "r");

        fprintf(FP, "\t C23:\n\t %s.C23\n\n", MFILE);
        sprintf(fname_C23, "%s.C23", MFILE);
        fp_C23 = fopen(fname_C23, "r");

        bool cond_1 = fp_C11 != NULL && fp_C22 != NULL && fp_C33 != NULL;
        bool cond_2 = fp_C44 != NULL && fp_C55 != NULL && fp_C66 != NULL;
        bool cond_3 = fp_C12 != NULL && fp_C13 != NULL && fp_C23 != NULL;

        bool cond_no_1 = fp_C11 == NULL && fp_C22 == NULL && fp_C33 == NULL;
        bool cond_no_2 = fp_C44 == NULL && fp_C55 == NULL && fp_C66 == NULL;
        bool cond_no_3 = fp_C12 == NULL && fp_C13 == NULL && fp_C23 == NULL;

        if (cond_1 && cond_2 && cond_3) {
            flag_cij_files_avail = true;
        } else if (cond_no_1 && cond_no_2 && cond_no_3) {
            flag_cij_files_avail = false;
        } else {
            if (fp_C11 == NULL) err("Could not open model file for C11!");
            if (fp_C22 == NULL) err("Could not open model file for C22!");
            if (fp_C33 == NULL) err("Could not open model file for C33!");
            if (fp_C44 == NULL) err("Could not open model file for C44!");
            if (fp_C55 == NULL) err("Could not open model file for C55!");
            if (fp_C66 == NULL) err("Could not open model file for C66!");
            if (fp_C12 == NULL) err("Could not open model file for C12!");
            if (fp_C13 == NULL) err("Could not open model file for C13!");
            if (fp_C23 == NULL) err("Could not open model file for C23!");
        }
    }

    /*elastic simulation */
    if (L == 0) {
        float tmp;
        float epsx = 0.0f;
        float epsy = 0.0f;
        float delx = 0.0f;
        float dely = 0.0f;
        float delxy = 0.0f;
        float gamx = 0.0f;
        float gamy = 0.0f;
        /* loop over global grid */
        for (k = 1; k <= NZG; k++) {
            for (i = 1; i <= NXG; i++) {
                for (j = 1; j <= NYG; j++) {
                    Rho = readdsk(fp_rho, format);

                    if (flag_velocity_files_avail) {
                        Vp = readdsk(fp_vp, format);
                        Vs = readdsk(fp_vs, format);

                        muv = Vs * Vs * Rho;
                        piv = Vp * Vp * Rho;
                    }

                    if (flag_cij_files_avail) {
                        C_11 = readdsk(fp_C11, format);
                        C_22 = readdsk(fp_C22, format);
                        C_33 = readdsk(fp_C33, format);
                        C_44 = readdsk(fp_C44, format);
                        C_55 = readdsk(fp_C55, format);
                        C_66 = readdsk(fp_C66, format);
                        C_12 = readdsk(fp_C12, format);
                        C_13 = readdsk(fp_C13, format);
                        C_23 = readdsk(fp_C23, format);
                    } else {
                        if (fp_epsx != NULL) {
                            epsx = readdsk(fp_epsx, format);
                        }
                        if (fp_epsy != NULL) {
                            epsy = readdsk(fp_epsy, format);
                        }
                        if (fp_delx != NULL) {
                            delx = readdsk(fp_delx, format);
                        }
                        if (fp_dely != NULL) {
                            dely = readdsk(fp_dely, format);
                        }
                        if (fp_delxy != NULL)
                        {
                            delxy = readdsk(fp_delxy, format);
                        }
                        if (fp_gamx != NULL) {
                            gamx = readdsk(fp_gamx, format);
                        }
                        if (fp_gamy != NULL) {
                            gamy = readdsk(fp_gamy, format);
                        }
                        // clang-format off
                        C_33 = Rho * Vp * Vp;
                        C_55 = Rho * Vs * Vs;
                        C_66 = (1 + 2 * gamx) * C_55;
                        C_11 = (1 + 2 * epsy) * C_33;
                        C_44 = C_66 / (1 + 2 * gamy);
                        C_22 = (1 + 2 * epsx) * C_33;
                        tmp = C_33 - C_55;
                        C_13 = -C_55 + sqrt(2*dely * C_33*tmp + tmp*tmp);
                        tmp = C_11 - C_66;
                        C_12 = -C_66 + sqrt(2*delxy * C_11*tmp + tmp*tmp);
                        tmp = C_33 - C_44;
                        C_23 = -C_44 + sqrt(2*delx * C_33*tmp + tmp*tmp);
                        // clang-format on
                    }

                    /* only the PE which belongs to the current global gridpoint
						is saving model parameters in his local arrays */
                    if ((POS[1] == ((i - 1) / NX)) &&
                        (POS[2] == ((j - 1) / NY)) &&
                        (POS[3] == ((k - 1) / NZ))) {
                        ii = i - POS[1] * NX;
                        jj = j - POS[2] * NY;
                        kk = k - POS[3] * NZ;

                        rho[jj][ii][kk] = Rho;

                        if (flag_velocity_files_avail) {
                            u[jj][ii][kk] = muv;
                            pi[jj][ii][kk] = piv;
                        }


                        C11[jj][ii][kk] = C_11;
                        C33[jj][ii][kk] = C_22;
                        C22[jj][ii][kk] = C_33;

                        C44[jj][ii][kk] = C_44;
                        C66[jj][ii][kk] = C_55;
                        C55[jj][ii][kk] = C_66;

                        C13[jj][ii][kk] = C_12;
                        C12[jj][ii][kk] = C_13;
                        C23[jj][ii][kk] = C_23;


                        if (WRITE_MODELFILES) {
                            pwavemod[jj][ii][kk] = Vp;
                            swavemod[jj][ii][kk] = Vs;
                        }
                    }
                }
            }
        }
    }

    if (L) { /*viscoelastic simulation */


        /* if TAU is specified in input file, q-files will NOT be read-in separately
		thus a constant q-model will be assumed */
		if (TAU==0.0) {
			fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
			sprintf(filename,"%s.qp",MFILE);
			fp_qp=fopen(filename,"r");
			if (fp_qp==NULL) err(" Could not open model file for Qp-values ! ");

			fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
			sprintf(filename,"%s.qs",MFILE);
			fp_qs=fopen(filename,"r");
			if (fp_qs==NULL) err(" Could not open model file for Qs-values ! ");
		}

		/* vector for maxwellbodies */
		pts=vector(1,L);
		for (l=1;l<=L;l++) {
			pts[l]=1.0/(2.0*PI*FL[l]);
			eta[l]=DT/pts[l];
		}
		/* in the viscoelastic case : reference frequency where no velocity dispersion occurs.
		 * if FREF is not given in input file, the largest center source frequency FC
		 * as specified in input file is used (not that the relation : FC=1/TS) is used here)*/
		if (FREF==0.0) ws=2.0*PI/TS;
		else ws=2.0*PI*FREF;

		/* loop over global grid */
		for (k=1;k<=NZG;k++){
			for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){
					Vp=readdsk(fp_vp, format);
					Vs=readdsk(fp_vs, format);
					Rho=readdsk(fp_rho , format);

					/*calculation of taus and taup by read-in q-files*/
					if (TAU==0.0) {
						Qp=readdsk(fp_qp, format);
						Qs=readdsk(fp_qs, format);
					}
					else {
						/*constant q (damping) case:*/
						Qp=2.0/TAU;
						Qs=2.0/TAU;
					}

					sumu=0.0;
					sumpi=0.0;
					for (int l=1;l<=L;l++){
						sumu=sumu+((ws*ws*pts[l]*pts[l]*(2/Qs))/(1.0+ws*ws*pts[l]*pts[l]));
						sumpi=sumpi+((ws*ws*pts[l]*pts[l]*(2/Qp))/(1.0+ws*ws*pts[l]*pts[l]));
					}

					muv=Vs*Vs*Rho/(1.0+sumu);
					piv=Vp*Vp*Rho/(1.0+sumpi);

					/* only the PE which belongs to the current global gridpoint
						is saving model parameters in his local arrays */
                    if ((POS[1] == ((i - 1) / NX)) &&
                        (POS[2] == ((j - 1) / NY)) &&
                        (POS[3] == ((k - 1) / NZ))) {
                        ii = i - POS[1] * NX;
                        jj = j - POS[2] * NY;
                        kk = k - POS[3] * NZ;

                        u[jj][ii][kk] = muv;
                        rho[jj][ii][kk] = Rho;
                        pi[jj][ii][kk] = piv;

                        taus[jj][ii][kk] = 2.0 / Qs;
                        taup[jj][ii][kk] = 2.0 / Qp;

                        if (WRITE_MODELFILES) {
                            pwavemod[jj][ii][kk] = Vp;
                            swavemod[jj][ii][kk] = Vs;
                            qsmod[jj][ii][kk] = 2 / taus[jj][ii][kk];
                            qpmod[jj][ii][kk] = 2 / taup[jj][ii][kk];
                        }
                    }
                }
            }
        }
    }


    fclose(fp_vp);
    fclose(fp_vs);
    fclose(fp_rho);
    if ((L) && (TAU == 0.0)) {
        fclose(fp_qp);
        fclose(fp_qs);
    }

    /* each PE writes his model to disk */

    if (WRITE_MODELFILES) {
        sprintf(filename, "%s.SOFI3D.pi", MFILE);
        writemod(filename, pi, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.u", MFILE);
        writemod(filename, u, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.vp", MFILE);
        writemod(filename, pwavemod, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.vs", MFILE);
        writemod(filename, swavemod, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        // Notice that the stiffness parameters are written
        // to disk in the conventional notation (third axis is vertical);
        // that's why there is a mismatch between filenames and variable names.
        sprintf(filename, "%s.SOFI3D.C11", MFILE);
        writemod(filename, C11, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C22", MFILE);
        writemod(filename, C33, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C33", MFILE);
        writemod(filename, C22, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C44", MFILE);
        writemod(filename, C44, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C55", MFILE);
        writemod(filename, C66, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C66", MFILE);
        writemod(filename, C55, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C12", MFILE);
        writemod(filename, C13, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C13", MFILE);
        writemod(filename, C12, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.C23", MFILE);
        writemod(filename, C23, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);
    }

    sprintf(filename, "%s.SOFI3D.rho", MFILE);
    writemod(filename, rho, 3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYID == 0) mergemod(filename, 3);

    if ((L) && (WRITE_MODELFILES)) {
        sprintf(filename, "%s.SOFI3D.qp", MFILE);
        writemod(filename, qpmod, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);

        sprintf(filename, "%s.SOFI3D.qs", MFILE);
        writemod(filename, qsmod, 3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID == 0) mergemod(filename, 3);
    }


    free_vector(pts, 1, L);
    if (WRITE_MODELFILES) {
        free_f3tensor(pwavemod, 0, NY + 1, 0, NX + 1, 0, NZ + 1);
        free_f3tensor(swavemod, 0, NY + 1, 0, NX + 1, 0, NZ + 1);
        if ((L) && (TAU == 0.0)) {
            free_f3tensor(qpmod, 0, NY + 1, 0, NX + 1, 0, NZ + 1);
            free_f3tensor(qsmod, 0, NY + 1, 0, NX + 1, 0, NZ + 1);
        }
    }
}
