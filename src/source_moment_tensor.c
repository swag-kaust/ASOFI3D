#include "data_structures.h"
#include "enum.h"
#include "fd.h"

/**
 * Add general moment tensor source for simulation of earthquakes.
 *
 * Source term is added to the stress components.
 * This source type is more general than the one computed in `eqsource`
 * function as here the moment tensor components can be specified directly
 * instead of being computed using other parameters.
 * Components of the moment tensor are specified via parameter file
 * and passed to this function via global variables.
 *
 * nt          Time step number
 * s           Stress tensor
 * srcpos_loc  Source positions (local to the MPI process)
 * nsrc        Number of sources
 * stype       Types of sources
 */
void source_moment_tensor(
        int nt, Tensor3d *s, float **srcpos_loc,
        float **signals, int nsrc, int *stype)
{
    extern float DT, DX, DY, DZ;
    extern float AMON;
    extern float M11, M12, M13;
    extern float M22, M23, M33;
    int i, j, k, l;
    float amp, scale_amp;
    ;
    float m11, m12, m13, m22, m23, m33;

    float ***sxx = s->xx;
    float ***syy = s->yy;
    float ***szz = s->zz;
    float ***sxy = s->xy;
    float ***syz = s->yz;
    float ***sxz = s->xz;

    /* adding source wavelet to stress components
           (moment tensor source) at source points */
    for (l = 1; l <= nsrc; l++)
    {
        if (stype[l] == SOURCE_TYPE_MOMENT_TENSOR)
        {
            i = (int)srcpos_loc[1][l];
            j = (int)srcpos_loc[2][l];
            k = (int)srcpos_loc[3][l];

            amp = signals[l][nt];

            scale_amp = DT / (DX * DY * DZ);

            amp = AMON * amp * scale_amp;

            m33 = M33;
            m13 = M13;
            m23 = M23;
            m11 = M11;
            m12 = M12;
            m22 = M22;

            sxx[j][i][k] += amp * m11;
            syy[j][i][k] += amp * m22;
            szz[j][i][k] += amp * m33;

            sxy[j][i][k] += 0.25 * amp * m12;
            sxy[j][i - 1][k] += 0.25 * amp * m12;
            sxy[j - 1][i][k] += 0.25 * amp * m12;
            sxy[j - 1][i - 1][k] += 0.25 * amp * m12;

            syz[j][i][k] += 0.25 * amp * m23;
            syz[j][i][k - 1] += 0.25 * amp * m23;
            syz[j][i - 1][k] += 0.25 * amp * m23;
            syz[j][i - 1][k - 1] += 0.25 * amp * m23;

            sxz[j][i][k] += 0.25 * amp * m13;
            sxz[j - 1][i][k] += 0.25 * amp * m13;
            sxz[j][i][k - 1] += 0.25 * amp * m13;
            sxz[j - 1][i][k - 1] += 0.25 * amp * m13;
        }
    }
}
