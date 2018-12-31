/*------------------------------------------------------------------------
 *   merge snapshots files written by the different processes to 
 *   a single file                                 
  *
 *  ----------------------------------------------------------------------*/

#include "enum.h"
#include "fd.h"
#include "globvar.h"


void merge(int nsnap, int type)
{
    extern char SNAP_FILE[STRING_SIZE];
    extern int NXG, NYG, SNAP_FORMAT, NPROCX, NPROCY, NPROCZ;
    extern int NX, NY, NZ, IDX, IDY, IDZ;
    extern FILE *FP;

    char file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[10];
    FILE *fp[NPROCY_MAX][NPROCX_MAX][NPROCZ_MAX], *fpout;
    int i, j, k, ip, jp, kp, n;
    float a;

    if ((NPROCX > NPROCX_MAX) || (NPROCY > NPROCY_MAX) || (NPROCZ > NPROCZ_MAX))
        err(" merge.c: constant expression NPROC?_MAX < NPROC? ");

    switch (SNAP_FORMAT) {
        case 1:
            sprintf(ext, ".su");
            break;
        case 2:
            sprintf(ext, ".asc");
            break;
        case 3:
            sprintf(ext, ".bin");
            break;
    }

    switch (type) {
        case 1:
            fprintf(FP, "x-component of particle velocity");
            strcat(ext, ".vx");
            break;
        case 2:
            fprintf(FP, "y-component of particle velocity");
            strcat(ext, ".vy");
            break;
        case 3:
            fprintf(FP, "z-component of particle velocity");
            strcat(ext, ".vz");
            break;
        case 4:
            fprintf(FP, "P-wave energyfield");
            strcat(ext, ".div");
            break;
        case 5:
            fprintf(FP, "S-wave energyfield");
            strcat(ext, ".curl");
            break;
        case 6:
            fprintf(FP, "pressure");
            strcat(ext, ".p");
            break;
        default:
            err("[%s] Incorrect type of the snapshot quantity"
                "Correct types are from 1 to 6", __func__);
            break;
    }

    sprintf(mfile, "%s%s", SNAP_FILE, ext);
    fprintf(FP, " (files: %s.x.y.z).\n", mfile);

    sprintf(outfile, "%s%s", SNAP_FILE, ext);
    fprintf(FP, "Writing merged snapshot file to %s\n", outfile);
    fpout = fopen(outfile, "w");

    fprintf(FP, "Opening snapshot files with template %s.x.y.z\n", mfile);
    for (kp = 0; kp <= NPROCZ - 1; kp++)
        for (ip = 0; ip <= NPROCX - 1; ip++)
            for (jp = 0; jp <= NPROCY - 1; jp++)
            {
                sprintf(file, "%s.%i.%i.%i", mfile, ip, jp, kp);
                fprintf(FP, "%s\n", file);
                fp[jp][ip][kp] = fopen(file, "r");
                if (fp[jp][ip][kp] == NULL)
                {
                    err(file);
                }
            }
    fprintf(FP, " ... finished. \n");

    fprintf(FP, " Copying...");
    for (n = 0; n < nsnap; n++)
        for (kp = 0; kp <= NPROCZ - 1; kp++)
            for (k = 1; k <= NZ; k += IDZ)
                for (ip = 0; ip <= NPROCX - 1; ip++)
                    for (i = 1; i <= NX; i += IDX)
                        for (jp = 0; jp <= NPROCY - 1; jp++)
                            for (j = 1; j <= NY; j += IDY)
                            {
                                a = readdsk(fp[jp][ip][kp], SNAP_FORMAT);
                                writedsk(fpout, a, SNAP_FORMAT);
                            }
    fprintf(FP, " ... finished. \n");

    for (kp = 0; kp <= NPROCZ - 1; kp++)
        for (ip = 0; ip <= NPROCX - 1; ip++)
            for (jp = 0; jp <= NPROCY - 1; jp++) {
                fclose(fp[jp][ip][kp]);
            }

    fclose(fpout);

    if (SNAP_FORMAT == FILE_FORMAT_BINARY) {
        fprintf(FP, "Use command:\n");
        fprintf(FP, "xmovie n1=%d n2=%d < %s loop=1 label1=Y label2=X\n",
                ((NYG - 1) / IDY) + 1, ((NXG - 1) / IDY) + 1, outfile);
        fprintf(FP, "to play movie. \n");
        fprintf(FP, "=====================================================\n");
    }
}
