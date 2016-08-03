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
 *   Write one single amplitude on disk                                   
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/


void writedsk(FILE *fp_out, float amp, int format){



	switch(format){
                case 1 : /* SU*/ 
                        err(" Sorry, SU-format for snapshots not implemented yet. \n");
                        break;
		case 2 :  /*ASCII*/
                        fprintf(fp_out,"%e\n", amp); 
                        break;
		case 3 :   /* BINARY */

			fwrite(&amp, sizeof(float), 1, fp_out);
              		break;
	                
		default :
			printf(" Don't know the format for the snapshot-data !\n");
			err(" No output was written. ");
	}
}
