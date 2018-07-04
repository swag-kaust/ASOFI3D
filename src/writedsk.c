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
