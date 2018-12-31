#include "fd.h"
#include "globvar.h"


/*
 * Write local part of the source field to the file `source_field_file`.
 * 'Local part' means 'belonging to the current MPI process'.
 *
 * Parameters
 * ----------
 * source_field_file :
 *     Prefix of the name of the source-field file, e.g., 'source_field/test_'.
 * s :
 *     Quantity being written to file (source term in the equations).
 * format :
 *     File format.
 *
 * See also
 * --------
 * writedsk
 *     To see possible file formats.
 *
 */
void write_source_field(
        char source_field_file[STRING_SIZE],
        float ***s,
        int format) {
    // External (global) variables.
    extern int NX, NY, NZ, POS[4], IDX, IDY, IDZ;

    int i, j, k;
    FILE *fp;
    char file[STRING_SIZE];

    sprintf(file, "%s.%d.%d.%d", source_field_file, POS[1], POS[2], POS[3]);
    fp = fopen(file, "w");
    if (fp == NULL) {
        err2("Cannot open file %s for writing\n", source_field_file);
    }

    for (k = 1; k <= NZ; k += IDZ) {
        for (i = 1; i <= NX; i += IDX) {
            for (j = 1; j <= NY; j += IDY) {
                writedsk(fp, s[j][i][k], format);
            }
        }
    }

    fclose(fp);
}
