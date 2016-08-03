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
 *   store amplitudes (particle velocities or pressure) at receiver positions
     in arrays
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo_acoustic(int lsamp, int ntr, int **recpos, float **sectionvx, float **sectionvy, 
float **sectionvz, float **sectiondiv, float **sectioncurl, float **sectionp, 
float ***vx, float ***vy, float ***vz, float ***sxx, float ***pi){ 
		
	extern int SEISMO;	
	int  itr, ins, nxrec, nyrec, nzrec, i, j, k;
	float  dh24x, dh24y, dh24z; 
        float vxx, vyy, vzz;
	/*float amp, vzy, vxz, vyz, vzx, vxy, vyx;*/
	extern float DX, DY, DZ;


        ins=lsamp; /* changed from "ins=lsamp/NDT;" (neccessary after correction of the buggy ns in sofi3D.c) */
	dh24x=1.0/DX;
	dh24y=1.0/DY;
	dh24z=1.0/DZ;
	
	for (itr=1;itr<=ntr;itr++){
		nxrec=recpos[1][itr];
		nyrec=recpos[2][itr];
		nzrec=recpos[3][itr];
		switch (SEISMO){                          
		case 1 : sectionvx[itr][ins]=vx[nyrec][nxrec][nzrec];
			 sectionvy[itr][ins]=vy[nyrec][nxrec][nzrec];
			 sectionvz[itr][ins]=vz[nyrec][nxrec][nzrec]; 
			 break;
		case 2 : sectionp[itr][ins]=-3.0*sxx[nyrec][nxrec][nzrec];
			 break;
		case 3 :				
			i=nxrec; j=nyrec; k=nzrec;
			/*vxy=(-vx[j+2][i][k]+27.0*(vx[j+1][i][k]-vx[j][i][k])+vx[j-1][i][k])*(24.0*DY);
			vxz=(-vx[j][i][k+2]+27.0*(vx[j][i][k+1]-vx[j][i][k])+vx[j][i][k-1])*(24.0*DZ);
			vyx=(-vy[j][i+2][k]+27.0*(vy[j][i+1][k]-vy[j][i][k])+vy[j][i-1][k])*(24.0*DX);
			vyz=(-vy[j][i][k+2]+27.0*(vy[j][i][k+1]-vy[j][i][k])+vy[j][i][k-1])*(24.0*DZ);
			vzx=(-vz[j][i+2][k]+27.0*(vz[j][i+1][k]-vz[j][i][k])+vz[j][i-1][k])*(24.0*DX);
			vzy=(-vz[j+2][i][k]+27.0*(vz[j+1][i][k]-vz[j][i][k])+vz[j-1][i][k])*(24.0*DY);*/
			
			/*vxy=(vx[j+1][i][k]-vx[j][i][k])*(dh24y);
		        vxz=(vx[j][i][k+1]-vx[j][i][k])*(dh24z);
			vyx=(vy[j][i+1][k]-vy[j][i][k])*(dh24x);
			vyz=(vy[j][i][k+1]-vy[j][i][k])*(dh24z);
			vzx=(vz[j][i+1][k]-vz[j][i][k])*(dh24x);
			vzy=(vz[j+1][i][k]-vz[j][i][k])*(dh24y);
			
			amp=u[j][i][k]*((vyz-vzy)*fabs(vyz-vzy)+
					    (vzx-vxz)*fabs(vzx-vxz)+(vxy-vyx)*fabs(vxy-vyx));
			sectioncurl[itr][ins]=fsign(amp)*sqrt(fabs(amp));*/

			/*vxx=(-vx[j][i+1][k]+27.0*(vx[j][i][k]-vx[j][i-1][k])+vx[j][i-2][k])*(24.0*DX);
			vyy=(-vy[j+1][i][k]+27.0*(vy[j][i][k]-vy[j-1][i][k])+vy[j-2][i][k])*(24.0*DY);
			vzz=(-vz[j][i][k+1]+27.0*(vz[j][i][k]-vz[j][i][k-1])+vz[j][i][k-2])*(24.0*DZ);*/
			
			vxx=(vx[j][i][k]-vx[j][i-1][k])*(dh24x);
			vyy=(vy[j][i][k]-vy[j-1][i][k])*(dh24y);
			vzz=(vz[j][i][k]-vz[j][i][k-1])*(dh24z);
			
			sectiondiv[itr][ins]=(vxx+vyy+vzz)*sqrt(pi[j][i][k]);

			break;
		case 4 :				

			sectionvx[itr][ins]=vx[nyrec][nxrec][nzrec];
			sectionvy[itr][ins]=vy[nyrec][nxrec][nzrec];
			sectionvz[itr][ins]=vz[nyrec][nxrec][nzrec]; 
			sectionp[itr][ins]=-3.0*sxx[nyrec][nxrec][nzrec];
			i=nxrec; j=nyrec; k=nzrec;
			
			/*vxy=(-vx[j+2][i][k]+27.0*(vx[j+1][i][k]-vx[j][i][k])+vx[j-1][i][k])*(dh24y);
			vxz=(-vx[j][i][k+2]+27.0*(vx[j][i][k+1]-vx[j][i][k])+vx[j][i][k-1])*(dh24z);
			vyx=(-vy[j][i+2][k]+27.0*(vy[j][i+1][k]-vy[j][i][k])+vy[j][i-1][k])*(dh24x);
			vyz=(-vy[j][i][k+2]+27.0*(vy[j][i][k+1]-vy[j][i][k])+vy[j][i][k-1])*(dh24z);
			vzx=(-vz[j][i+2][k]+27.0*(vz[j][i+1][k]-vz[j][i][k])+vz[j][i-1][k])*(dh24x);
			vzy=(-vz[j+2][i][k]+27.0*(vz[j+1][i][k]-vz[j][i][k])+vz[j-1][i][k])*(dh24y);*/
			
			/*vxy=(vx[j+1][i][k]-vx[j][i][k])*(dh24y);
		        vxz=(vx[j][i][k+1]-vx[j][i][k])*(dh24z);
			vyx=(vy[j][i+1][k]-vy[j][i][k])*(dh24x);
			vyz=(vy[j][i][k+1]-vy[j][i][k])*(dh24z);
			vzx=(vz[j][i+1][k]-vz[j][i][k])*(dh24x);
			vzy=(vz[j+1][i][k]-vz[j][i][k])*(dh24y);
			
			amp=u[j][i][k]*((vyz-vzy)*fabs(vyz-vzy)+
					    (vzx-vxz)*fabs(vzx-vxz)+(vxy-vyx)*fabs(vxy-vyx));
			sectioncurl[itr][ins]=fsign(amp)*sqrt(fabs(amp));*/

			/*vxx=(-vx[j][i+1][k]+27.0*(vx[j][i][k]-vx[j][i-1][k])+vx[j][i-2][k])*(dh24x);
			vyy=(-vy[j+1][i][k]+27.0*(vy[j][i][k]-vy[j-1][i][k])+vy[j-2][i][k])*(dh24y);
			vzz=(-vz[j][i][k+1]+27.0*(vz[j][i][k]-vz[j][i][k-1])+vz[j][i][k-2])*(dh24z);*/
			
			vxx=(vx[j][i][k]-vx[j][i-1][k])*(dh24x);
			vyy=(vy[j][i][k]-vy[j-1][i][k])*(dh24y);
			vzz=(vz[j][i][k]-vz[j][i][k-1])*(dh24z);
			
			sectiondiv[itr][ins]=(vxx+vyy+vzz)*sqrt(pi[j][i][k]);
			break;
 
		}
	}
}
