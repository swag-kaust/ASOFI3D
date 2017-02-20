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
/* ----------------------------------------------------------------------
 * reads checkpoint file for the continuation of a simulation
 ----------------------------------------------------------------------*/
#include "fd.h"

void read_checkpoint(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *** vx, float *** vy, float *** vz,
                     float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz,
		     float *** rxx, float *** ryy,float *** rzz, float *** rxy, float *** ryz, float *** rxz,
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

				fread(&vx[j][i][k],sizeof(float),1,fp);
				fread(&vy[j][i][k],sizeof(float),1,fp);
				fread(&vz[j][i][k],sizeof(float),1,fp);
				fread(&sxx[j][i][k],sizeof(float),1,fp);
				fread(&syy[j][i][k],sizeof(float),1,fp);
				fread(&szz[j][i][k],sizeof(float),1,fp);
				fread(&sxy[j][i][k],sizeof(float),1,fp);
				fread(&syz[j][i][k],sizeof(float),1,fp);
				fread(&sxz[j][i][k],sizeof(float),1,fp);
			}
		}
	}

	
	if (L) {
		for (j=1; j<=NY; j++) {
			for (i=1; i<=NX; i++) {
				for (k=1; k<=NZ; k++) {
					fread(&rxx[j][i][k],sizeof(float),1,fp);
					fread(&ryy[j][i][k],sizeof(float),1,fp);
					fread(&rzz[j][i][k],sizeof(float),1,fp);
					fread(&rxy[j][i][k],sizeof(float),1,fp);
					fread(&ryz[j][i][k],sizeof(float),1,fp);
					fread(&rxz[j][i][k],sizeof(float),1,fp);
				}
			}
		}
	}
	
	
	
	if (ABS_TYPE == 1) {

		if (POS[1]==0) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=FW; i++) {
					for (k=1; k<=NZ; k++) {
						fread(&psi_sxx_x[j][i][k],sizeof(float),1,fp);
						fread(&psi_sxy_x[j][i][k],sizeof(float),1,fp);
						fread(&psi_sxz_x[j][i][k],sizeof(float),1,fp);
						fread(&psi_vxx[j][i][k],sizeof(float),1,fp);
						fread(&psi_vyx[j][i][k],sizeof(float),1,fp);
						fread(&psi_vzx[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

		if (POS[1]==NPROCX-1) {
			for (j=1; j<=NY; j++) {
				for (i=FW+1; i<=2*FW; i++) {
					for (k=1; k<=NZ; k++) {
						fread(&psi_sxx_x[j][i][k],sizeof(float),1,fp);
						fread(&psi_sxy_x[j][i][k],sizeof(float),1,fp);
						fread(&psi_sxz_x[j][i][k],sizeof(float),1,fp);
						fread(&psi_vxx[j][i][k],sizeof(float),1,fp);
						fread(&psi_vyx[j][i][k],sizeof(float),1,fp);
						fread(&psi_vzx[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

		if (POS[2]==0 && FREE_SURF==0) {
			for (j=1; j<=FW; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=NZ; k++) {
						fread(&psi_sxy_y[j][i][k],sizeof(float),1,fp);
						fread(&psi_syy_y[j][i][k],sizeof(float),1,fp);
						fread(&psi_syz_y[j][i][k],sizeof(float),1,fp);
						fread(&psi_vxy[j][i][k],sizeof(float),1,fp);
						fread(&psi_vyy[j][i][k],sizeof(float),1,fp);
						fread(&psi_vzy[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}
		if (POS[2]==NPROCY-1) {
			for (j=FW+1; j<=2*FW; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=NZ; k++) {
						fread(&psi_sxy_y[j][i][k],sizeof(float),1,fp);
						fread(&psi_syy_y[j][i][k],sizeof(float),1,fp);
						fread(&psi_syz_y[j][i][k],sizeof(float),1,fp);
						fread(&psi_vxy[j][i][k],sizeof(float),1,fp);
						fread(&psi_vyy[j][i][k],sizeof(float),1,fp);
						fread(&psi_vzy[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

		if (POS[3]==0) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=NX; i++) {
					for (k=1; k<=FW; k++) {
						fread(&psi_sxz_z[j][i][k],sizeof(float),1,fp);
						fread(&psi_syz_z[j][i][k],sizeof(float),1,fp);
						fread(&psi_szz_z[j][i][k],sizeof(float),1,fp);
						fread(&psi_vxz[j][i][k],sizeof(float),1,fp);
						fread(&psi_vyz[j][i][k],sizeof(float),1,fp);
						fread(&psi_vzz[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}
		if (POS[3]==NPROCZ-1) {
			for (j=1; j<=NY; j++) {
				for (i=1; i<=NX; i++) {
					for (k=FW+1; k<=2*FW; k++) {
						fread(&psi_sxz_z[j][i][k],sizeof(float),1,fp);
						fread(&psi_syz_z[j][i][k],sizeof(float),1,fp);
						fread(&psi_szz_z[j][i][k],sizeof(float),1,fp);
						fread(&psi_vxz[j][i][k],sizeof(float),1,fp);
						fread(&psi_vyz[j][i][k],sizeof(float),1,fp);
						fread(&psi_vzz[j][i][k],sizeof(float),1,fp);
					}
				}
			}
		}

	}

	fclose(fp);

}
