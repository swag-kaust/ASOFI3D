/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI3D.
 * 
 * SOFI3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI3D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
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
