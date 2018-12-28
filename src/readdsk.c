#include "enum.h"
#include "fd.h"


/*
 * Read one single float value from a given file.
 *
 * The given file can be in one of the given formats:
 * format = 1 : SU (IEEE)
 * format = 2 : ASCII
 * format = 3 : BINARY (IEEE)
 *
 * fp_in   Pointer to the opened file
 * format  File format
 */
float readdsk(FILE *fp_in, int format)
{
    float amp = 0.0;

    size_t nelems;

    switch (format) {
        case FILE_FORMAT_SU:
            err(" Sorry, SU-format for snapshots not implemented yet. \n");
            break;
        case FILE_FORMAT_ASCII:
            if (fscanf(fp_in, "%e\n", &amp) != 1) {
                err("[%s] Could not read an amplitude "
                    "from a file in ASCII format\n", __func__);
            }
            break;
        case FILE_FORMAT_BINARY:
            nelems = fread(&amp, sizeof(float), 1, fp_in);
            if (nelems != 1) {
                err("[%s] Could not read an amplitude "
                    "from a file in binary format\n", __func__);
            }
            break;

        default:
            err("[%s] Unsupported file format. "
                "Supported formats are SU, ASCII, and BINARY", __func__);
    }

    return amp;
}
