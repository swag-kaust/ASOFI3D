#include "fd.h"


/**
 * Merge together source-field data scattered over many files written by
 * separate MPI processes into a single file.
 */
void merge_source_field(char source_field_file[STRING_SIZE], int format) {
	extern int NXG, NYG, MYID, NPROCX, NPROCY, NPROCZ;
	extern int NX, NY, NZ, NPROC, IDX, IDY, IDZ;
	extern FILE *FP;


	char file[STRING_SIZE];
	FILE *fp[NPROCY_MAX][NPROCX_MAX][NPROCZ_MAX], *fp_out;
	int i, j, k, ip, jp, kp;
	float a;


	if ((NPROCX>NPROCX_MAX)||(NPROCY>NPROCY_MAX)||(NPROCZ>NPROCZ_MAX))
		err("merge_source_field.c: constant expression NPROC?_MAX < NPROC? ");

	printf("\n **Message from merge_source_field (printed by PE %d):\n", MYID);
	printf(" PE %d starts merge of %d source field files\n", MYID, NPROC);	

	fprintf(FP, " writing merged source-field file '%s'\n", source_field_file);
	fp_out = fopen(source_field_file, "w");

    for (kp = 0; kp <= NPROCZ-1; kp++) {
        for (ip = 0; ip <= NPROCX-1; ip++) {
            for (jp = 0; jp <= NPROCY-1; jp++) {
				sprintf(file, "%s.%i.%i.%i", source_field_file, ip, jp, kp);
				fp[jp][ip][kp] = fopen(file, "r");
				if (fp[jp][ip][kp]==NULL) {
                    err("merge: can't read source_field_file !");
                }
			}
        }
    }

	fprintf(FP, " ... finished. \n");
	fprintf(FP, " Copying...");

	for (kp = 0; kp <= NPROCZ-1; kp++)
		for (k = 1; k <= NZ; k += IDZ)
			for (ip = 0; ip <= NPROCX-1; ip++)
				for (i = 1; i <= NX; i += IDX)
					for (jp = 0; jp <= NPROCY-1; jp++)
						for (j = 1; j <= NY; j += IDY){
							a = readdsk(fp[jp][ip][kp],format);
							writedsk(fp_out, a, format);
						}

	fprintf(FP," finished. \n");

	for (kp=0;kp<=NPROCZ-1; kp++)
		for (ip=0;ip<=NPROCX-1; ip++)
			for (jp=0;jp<=NPROCY-1; jp++){
				fclose(fp[jp][ip][kp]);
			}


	fclose(fp_out);

	fprintf(FP," Use \n");
	fprintf(FP," xmovie n1=%d n2=%d < %s loop=1 label1=\"Y (vertical)\" label2=\"X (horizontal)\" title=%%g \n",
			((NYG-1)/IDY)+1,((NXG-1)/IDX)+1,source_field_file);
	fprintf(FP," to visualize source field. \n");
}
