/* ----------------------------------------------------------------------
 * reads checkpoint file for the continuation of a simulation
 ----------------------------------------------------------------------*/
#include "data_structures.h"
#include "fd.h"
#include "globvar.h"


static void read_value(float *dst, FILE *fp);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, Velocity *v,
                     Tensor3d *s,
		     Tensor3d *r,
		     float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y,
		     float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z,
                     float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy,
                     float *** psi_vxz, float *** psi_vyz, float *** psi_vzz) {

	int i,j,k;
	char myid[5];
	FILE *fp;
	char checkptfile[STRING_SIZE];
	extern int MYID;
	extern int ABS_TYPE, NX,NY,NZ,FW,POS[4],NPROCX,NPROCY,NPROCZ,FREE_SURF;
	extern int L;
	extern char  CHECKPTFILE[STRING_SIZE];
        float ***vx = v->x;
        float ***vy = v->y;
        float ***vz = v->z;

        float ***sxx = s->xx;
        float ***syy = s->yy;
        float ***szz = s->zz;
        float ***sxy = s->xy;
        float ***syz = s->yz;
        float ***sxz = s->xz;

        float ***rxx = r->xx;
        float ***ryy = r->yy;
        float ***rzz = r->zz;
        float ***rxy = r->xy;
        float ***ryz = r->yz;
        float ***rxz = r->xz;
	
	sprintf(checkptfile,"%s",CHECKPTFILE);
	sprintf(myid,".%d",MYID);
	strcat(checkptfile,myid);



	fp=fopen(checkptfile,"rb");
	if (fp==NULL) {
		err("CHECKPTFILE can't be opened !");
	}

	for (j=ny1; j<=ny2; j++) {
		for (i=nx1; i<=nx2; i++) {
			for (k=nz1; k<=nz2; k++) {
				read_value(&vx[j][i][k], fp);
				read_value(&vy[j][i][k], fp);
				read_value(&vz[j][i][k], fp);

				read_value(&sxx[j][i][k], fp);
				read_value(&syy[j][i][k], fp);
				read_value(&szz[j][i][k], fp);
				read_value(&sxy[j][i][k], fp);
				read_value(&syz[j][i][k], fp);
				read_value(&sxz[j][i][k], fp);
			}
		}
	}

	
    // Viscoelastic data.
	if (L) {
		for (j=1; j<=NY; j++) {
			for (i=1; i<=NX; i++) {
				for (k=1; k<=NZ; k++) {
					read_value(&rxx[j][i][k], fp);
					read_value(&ryy[j][i][k], fp);
					read_value(&rzz[j][i][k], fp);
					read_value(&rxy[j][i][k], fp);
					read_value(&ryz[j][i][k], fp);
					read_value(&rxz[j][i][k], fp);
				}
			}
		}
	}
	
	
	
	if (ABS_TYPE == 1) {
		if (POS[1]==0) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=FW; i++) {
					for (k=1; k<=NZ; k++) {
						read_value(&psi_sxx_x[j][i][k], fp);
						read_value(&psi_sxy_x[j][i][k], fp);
						read_value(&psi_sxz_x[j][i][k], fp);

						read_value(&psi_vxx[j][i][k], fp);
						read_value(&psi_vyx[j][i][k], fp);
						read_value(&psi_vzx[j][i][k], fp);
					}
				}
			}
		}

		if (POS[1]==NPROCX-1) {
			for (j=1; j<=NY; j++) {
				for (i=FW+1; i<=2*FW; i++) {
					for (k=1; k<=NZ; k++) {
						read_value(&psi_sxx_x[j][i][k], fp);
						read_value(&psi_sxy_x[j][i][k], fp);
						read_value(&psi_sxz_x[j][i][k], fp);

						read_value(&psi_vxx[j][i][k], fp);
						read_value(&psi_vyx[j][i][k], fp);
						read_value(&psi_vzx[j][i][k], fp);
					}
				}
			}
		}

		if (POS[2]==0 && FREE_SURF==0) {
			for (j=1; j<=FW; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=NZ; k++) {
						read_value(&psi_sxy_y[j][i][k], fp);
						read_value(&psi_syy_y[j][i][k], fp);
						read_value(&psi_syz_y[j][i][k], fp);

						read_value(&psi_vxy[j][i][k], fp);
						read_value(&psi_vyy[j][i][k], fp);
						read_value(&psi_vzy[j][i][k], fp);
					}
				}
			}
		}
		if (POS[2]==NPROCY-1) {
			for (j=FW+1; j<=2*FW; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=NZ; k++) {
						read_value(&psi_sxy_y[j][i][k], fp);
						read_value(&psi_syy_y[j][i][k], fp);
						read_value(&psi_syz_y[j][i][k], fp);

						read_value(&psi_vxy[j][i][k], fp);
						read_value(&psi_vyy[j][i][k], fp);
						read_value(&psi_vzy[j][i][k], fp);
					}
				}
			}
		}

		if (POS[3]==0) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=FW; k++) {
						read_value(&psi_sxz_z[j][i][k], fp);
						read_value(&psi_syz_z[j][i][k], fp);
						read_value(&psi_szz_z[j][i][k], fp);

						read_value(&psi_vxz[j][i][k], fp);
						read_value(&psi_vyz[j][i][k], fp);
						read_value(&psi_vzz[j][i][k], fp);
					}
				}
			}
		}
		if (POS[3]==NPROCZ-1) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=NX; i++) {
					for (k=FW+1; k<=2*FW; k++) {
						read_value(&psi_sxz_z[j][i][k], fp);
						read_value(&psi_syz_z[j][i][k], fp);
						read_value(&psi_szz_z[j][i][k], fp);

						read_value(&psi_vxz[j][i][k], fp);
						read_value(&psi_vyz[j][i][k], fp);
						read_value(&psi_vzz[j][i][k], fp);
					}
				}
			}
		}
	} // end if (ABS_TYPE == 1)

	fclose(fp);
}

/**
 * Read a float value from file fp and write it into dst.
 */
static void read_value(float *dst, FILE *fp)
{
    extern char CHECKPTFILE[STRING_SIZE];

    size_t nelems = fread(dst, sizeof(float), 1, fp);

    if (nelems != 1) {
        char msg[STRING_SIZE];
        sprintf(msg,
                "Error occurred while reading from a checkpoint file '%s'",
                CHECKPTFILE);
        err(msg);
    }
}
