/*------------------------------------------------------------------------
 *   Read extern source wavelet
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float *rd_sour(int *nts,FILE* fp_source){

	/* local variables */
	float *psource;
	int i, c;
	extern int NT;

	if (fp_source==NULL) err(" SIGNAL_FILE could no be opened as a required by SOURCE_TYPE is set to 4 !");
	/* fscanf(fp_source,"%i", nts); */
        *nts=0;
        while ((c=fgetc(fp_source)) != EOF)
         if (c=='\n') ++(*nts);
        rewind(fp_source);
	psource=vector(1,NT);
    size_t ret;
	for (i=1;i<=*nts;i++) {
        ret = fscanf(fp_source,"%e",&psource[i]);
        if (ret != 1) {
            fclose(fp_source);
            err("Error occurred while reading a signal file");
        }
    }
	fclose(fp_source);
	return psource;
}
