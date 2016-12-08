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
 *   output source signal e.g. for cross-correlation, deconvolution 
 *   or comparison with analytical solutions                                  
 *
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form){

	/* declaration of extern variables */
	extern float DT;
	extern int NDT;
	
	/* declaration of local variables */
	int i,j, ntr=1;
	segy tr;
	int tracl ;
	
	
	
	switch(seis_form){
	case 1 :
		for(tracl=1;tracl<=ntr;tracl++){        /*SEGY (without file-header)*/
         		
			tr.tracl=tracl;     
			tr.trid=(short)1;           /* trace identification code: 1=seismic*/
			
			tr.ns=(unsigned short)iround(ns/NDT); /* number of samples in this trace */
			tr.dt=(unsigned short)iround(((float)NDT*DT)*1.0e6); /* sample interval in micro-seconds */
			tr.d1=(float)NDT*DT;        /* sample spacing for non-seismic data */

			for(j=1;j<=tr.ns;j++) tr.data[j]=signals[tracl][j*NDT];

			fwrite(&tr,240,1,fp);
			fwrite(&tr.data[1],4,tr.ns,fp);
		}
		break;


	case 2 :
		for(i=1;i<=ntr;i++){         /*ASCII ONE COLUMN*/
			for(j=1;j<=ns;j+=NDT) fprintf(fp,"%e\n", signals[i][j]);
		}
		break;

	case 3 :                             /*BINARY */

		for(i=1;i<=ntr;i++)
			for(j=1;j<=ns;j+=NDT){
				fwrite(&signals[i][j],sizeof(float),1,fp); }
		break;

	default :
		fprintf(stdout," Message from output_source_signal: Don't know data format for seismograms !\n");
		fprintf(stdout," No output written. ");
	}

	fclose(fp);
}
