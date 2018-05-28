/*----------------------------------------------------------------------
 *  saves checkpoint file for a later continuation of a simulation
 *
 ---------------------------------------------------------------------- */
#include "fd.h"
#include "data_structures.h"

void save_checkpoint(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
        Velocity *v, Tensor3d *s, Tensor3d *r,
        float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y,
        float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z,
        float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy,
        float *** psi_vxz, float *** psi_vyz, float *** psi_vzz) {

	int i,j, k;
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



	fp=fopen(checkptfile,"wb");
	if (fp==NULL) {
		err("CHECKPTFILE can't be opened !");
	}


	for (j=ny1; j<=ny2; j++) {
		for (i=nx1; i<=nx2; i++) {
			for (k=nz1; k<=nz2; k++) {
				fwrite(&vx[j][i][k],sizeof(float),1,fp);
				fwrite(&vy[j][i][k],sizeof(float),1,fp);
				fwrite(&vz[j][i][k],sizeof(float),1,fp);
				fwrite(&sxx[j][i][k],sizeof(float),1,fp);
				fwrite(&syy[j][i][k],sizeof(float),1,fp);
				fwrite(&szz[j][i][k],sizeof(float),1,fp);
				fwrite(&sxy[j][i][k],sizeof(float),1,fp);
				fwrite(&syz[j][i][k],sizeof(float),1,fp);
				fwrite(&sxz[j][i][k],sizeof(float),1,fp);
			}
		}
	}

	if (L) {
		for (j=1; j<=NY; j++) {
			for (i=1; i<=NX; i++) {
				for (k=1; k<=NZ; k++) {
					fwrite(&rxx[j][i][k],sizeof(float),1,fp);
					fwrite(&ryy[j][i][k],sizeof(float),1,fp);
					fwrite(&rzz[j][i][k],sizeof(float),1,fp);
					fwrite(&rxy[j][i][k],sizeof(float),1,fp);
					fwrite(&ryz[j][i][k],sizeof(float),1,fp);
					fwrite(&rxz[j][i][k],sizeof(float),1,fp);
				}
			}
		}
	}

	if (ABS_TYPE == 1) {

		if (POS[1]==0) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=FW; i++) {
					for (k=1; k<=NZ; k++) {
						fwrite(&psi_sxx_x[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_sxy_x[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_sxz_x[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vxx[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vyx[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vzx[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

		if (POS[1]==NPROCX-1) {
			for (j=1; j<=NY; j++) {
				for (i=FW+1; i<=2*FW; i++) {
					for (k=1; k<=NZ; k++) {
						fwrite(&psi_sxx_x[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_sxy_x[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_sxz_x[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vxx[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vyx[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vzx[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

		if (POS[2]==0 && FREE_SURF==0) {
			for (j=1; j<=FW; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=NZ; k++) {
						fwrite(&psi_sxy_y[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_syy_y[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_syz_y[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vxy[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vyy[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vzy[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}
		if (POS[2]==NPROCY-1) {
			for (j=FW+1; j<=2*FW; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=NZ; k++) {
						fwrite(&psi_sxy_y[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_syy_y[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_syz_y[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vxy[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vyy[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vzy[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

		if (POS[3]==0) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=FW; k++) {
						fwrite(&psi_sxz_z[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_syz_z[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_szz_z[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vxz[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vyz[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vzz[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}
		if (POS[3]==NPROCZ-1) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=NX; i++) {
					for (k=FW+1; k<=2*FW; k++) {
						fwrite(&psi_sxz_z[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_syz_z[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_szz_z[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vxz[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vyz[j][i][k],sizeof(float),1,fp);
						fwrite(&psi_vzz[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

	}

	fclose(fp);

}















