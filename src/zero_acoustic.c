/*------------------------------------------------------------------------
 *   Initialization of the wave field with zero values (zero wavefield)
 *  
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "data_structures.h"

void zero_acoustic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, Velocity *v,
float *** sxx){



	register int i, j, k;

        float ***vx = v->x;
        float ***vy = v->y;
        float ***vz = v->z;

	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
				vx[j][i][k]=0.0;
				vy[j][i][k]=0.0;
				vz[j][i][k]=0.0;
				sxx[j][i][k]=0.0;
				}
			}
		}
	
}
