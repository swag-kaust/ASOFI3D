/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters to FILE *fp */
void info(FILE *fp){

	fprintf(fp," \n\n\n***********************************************************\n");
	fprintf(fp," This is program ASOFI3D.          \n");
	fprintf(fp," Parallel 3-D anisotropic elastic/viscoelastic Finite Difference Modelling      \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," Distributed by									    \n");
	fprintf(fp," SWAG,         \n");
	fprintf(fp," KAUST         \n");
	fprintf(fp," isotropic version from: http://www.gpi.kit.edu \n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");
}
