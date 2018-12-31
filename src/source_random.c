#include "enum.h"
#include "data_structures.h"
#include "fd.h"
#include "globvar.h"


/*
 * Generate random P-wave source at all grid nodes.
 * This special source term is used for checking self-adjointness
 * of the code through dot-product test.
 * See, e.g.,
 * http://sepwww.stanford.edu/sep/prof/pvi/conj/paper_html/node9.html
 * for details.
 *
 * Parameters
 * ----------
 * nt :
 *     Time step number (starting with 1).
 * s :
 *     Stress tensor on the local grid.
 * source_field :
 *     Source field on the local grid.
 */
void source_random(int nt, Tensor3d *s, float ***source_field)
{
    extern int NX, NY, NZ;
    extern int IDX, IDY, IDZ;
    extern float DX, DY, DZ;
    extern float DT;
    extern int MYID;
    extern int SOURCE_TYPE;

    float amp = 0.0;

    float random_signal;

    float ***sxx = s->xx;
    float ***syy = s->yy;
    float ***szz = s->zz;

    // If user specifies SOURCE_TYPE = 7 in the parameter file,
    // then generate random source value at the each point of the local grid
    // at the initial time step.
    srand((unsigned int) time(NULL));
    if (SOURCE_TYPE == SOURCE_TYPE_RANDOM) {
        if (nt == 1) {
            for (int k = 1; k <= NZ; k += IDZ) {
                for (int i = 1; i <= NX; i += IDX) {
                    for (int j = 1; j <= NY; j += IDY) {
                        random_signal = 1e-2 * ((float) rand()) / RAND_MAX;
                        // scaled explosive source with respect to spatial
                        // and temporal discretization, seismic Moment = 1 Nm
                        amp = DT * random_signal / (DX * DY * DZ);

                        source_field[j][i][k] = amp;

                        sxx[j][i][k] -= amp;
                        syy[j][i][k] -= amp;
                        szz[j][i][k] -= amp;
                    }
                }
            }

            // Write source data for the current time step to file.
            // STRING_SIZE is #define
            char source_field_file[STRING_SIZE];
            sprintf(source_field_file, "source_field/%d.bin", nt);
            write_source_field(source_field_file, source_field, 3);
            MPI_Barrier(MPI_COMM_WORLD);
            if (MYID == 0) merge_source_field(source_field_file, 3);
        }
    }
}
