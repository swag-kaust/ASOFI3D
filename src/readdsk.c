/*------------------------------------------------------------------------
 *   Read one single amplitude from file
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/

float readdsk(FILE *fp_in, int format)
{
    float amp = 0.0;

    size_t nelems;

    switch (format) {
        case 1: /* SU*/
            err(" Sorry, SU-format for snapshots not implemented yet. \n");
            break;
        case 2: /*ASCII*/
            if (fscanf(fp_in, "%e\n", &amp) != 1) {
                err("[readdsk] Could not read an amplitude "
                    "from a file in ASCII format\n");
            }
            break;
        case 3: /* BINARY */
            nelems = fread(&amp, sizeof(float), 1, fp_in);
            if (nelems != 1) {
                err("[readdsk] Could not read an amplitude "
                    "from a file in binary format\n");
            }
            break;

        default:
            printf(" Don't know the format for the snapshot-data !\n");
            err(" No output was written. ");
    }

    return amp;
}
