/*------------------------------------------------------------------------
 *   generate P-wave source at source nodes
 *
 *  ----------------------------------------------------------------------*/

#include "data_structures.h"
#include "fd.h"
#include "globvar.h"


void psource(int nt, Tensor3d *s, float **srcpos_loc, float **signals, int nsrc, int *stype)
{
    extern float DX, DY, DZ;
    extern float DT;

    int i, j, k, l;
    float amp = 0.0;

    float ***sxx = s->xx;
    float ***syy = s->yy;
    float ***szz = s->zz;

    /* adding source wavelet to stress components 
	   (explosive source) at source points */

    for (l = 1; l <= nsrc; l++)
    {
        if (stype[l] == 1)
        {
            i = (int)srcpos_loc[1][l];
            j = (int)srcpos_loc[2][l];
            k = (int)srcpos_loc[3][l];

            // scaled explosive source with respect to spatial
            // and temporal discretization, seismic Moment = 1 Nm
            amp = DT * (signals[l][nt]) / (DX * DY * DZ);

            sxx[j][i][k] -= amp;
            syy[j][i][k] -= amp;
            szz[j][i][k] -= amp;
        }
    }
}
