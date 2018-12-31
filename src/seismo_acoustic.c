/*------------------------------------------------------------------------
 *   store amplitudes (particle velocities or pressure) at receiver positions
     in arrays
 *  ----------------------------------------------------------------------*/

#include "data_structures.h"
#include "fd.h"
#include "globvar.h"


// Macro ATTR_UNUSED allows to tell a compiler to stop complaining
// about an unused variable.
// Currently, only implemented for gcc compiler (which defines __GNUC__ macro
// internally).
#ifdef __GNUC__
    #define ATTR_UNUSED __attribute__((unused))
#else
    #define ATTR_UNUSED
#endif


void seismo_acoustic(int lsamp, int ntr, int **recpos,
        float **sectionvx, float **sectionvy, float **sectionvz,
        float **sectiondiv,
        float **sectioncurl ATTR_UNUSED,
        float **sectionp, 
        Velocity *v, float ***sxx, float ***pi) {

	extern int SEISMO;	
	int  itr, ins, nxrec, nyrec, nzrec, i, j, k;
	float  dh24x, dh24y, dh24z; 
        float vxx, vyy, vzz;
	/*float amp, vzy, vxz, vyz, vzx, vxy, vyx;*/
	extern float DX, DY, DZ;

    float ***vx = v->x;
    float ***vy = v->y;
    float ***vz = v->z;


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
			
			vxx=(vx[j][i][k]-vx[j][i-1][k])*(dh24x);
			vyy=(vy[j][i][k]-vy[j-1][i][k])*(dh24y);
			vzz=(vz[j][i][k]-vz[j][i][k-1])*(dh24z);
			
			sectiondiv[itr][ins]=(vxx+vyy+vzz)*sqrt(pi[j][i][k]);
			break;
 
		}
	}
}
