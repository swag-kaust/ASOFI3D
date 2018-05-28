/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters to FILE *fp */
void info(FILE *fp){

	fprintf(fp," \n\n\n***********************************************************\n");
	fprintf(fp," This is program SOFI3D.          \n");
	fprintf(fp," Parallel 3-D acoustic/elastic/viscoelastic Finite Difference Modelling      \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," Distributed by									    \n");
	fprintf(fp," Geophysical Institute, Department of Physics,         \n");
	fprintf(fp," Institute of Technology, Karlsruhe, Germany         \n");
	fprintf(fp," http://www.gpi.kit.edu \n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");
}
