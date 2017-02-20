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
 *   prepares to write seismograms from each PE individually to disk
 *------------------------------------------------------------------------*/

#include "fd.h"

void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionvz,
		float **sectionp, float **sectioncurl, float **sectiondiv,
		int  **recpos, int  **recpos_loc, int ntr, float ** srcpos, int ishot,int ns){

	extern int SEISMO, SEIS_FORMAT[6], MYID, RUN_MULTIPLE_SHOTS;	
	extern char  SEIS_FILE[STRING_SIZE];

	char vxf[STRING_SIZE], vyf[STRING_SIZE], vzf[STRING_SIZE], curlf[STRING_SIZE], divf[STRING_SIZE], pf[STRING_SIZE],file_ext[5];
	int nsrc=1;

	switch (SEIS_FORMAT[0]){
	case 0: sprintf(file_ext,"sgy"); break;
	case 1: sprintf(file_ext,"su");  break;
	case 2: sprintf(file_ext,"txt"); break;
	case 3: sprintf(file_ext,"bin"); break;
	case 4: sprintf(file_ext,"sgy"); break;
	case 5: sprintf(file_ext,"sgy"); break;
	}

	/*note that "y" denotes the vertical coordinate*/

	if (RUN_MULTIPLE_SHOTS){
		sprintf(vxf,"%s_vx.%s.shot%d.%d",SEIS_FILE,file_ext,ishot,MYID);
		sprintf(vyf,"%s_vy.%s.shot%d.%d",SEIS_FILE,file_ext,ishot,MYID);
		sprintf(vzf,"%s_vz.%s.shot%d.%d",SEIS_FILE,file_ext,ishot,MYID);
		sprintf(curlf,"%s_curl.%s.shot%d.%d",SEIS_FILE,file_ext,ishot,MYID);
		sprintf(divf,"%s_div.%s.shot%d.%d",SEIS_FILE,file_ext,ishot,MYID);
		sprintf(pf,"%s_p.%s.shot%d.%d",SEIS_FILE,file_ext,ishot,MYID);
	}
	else{
		sprintf(vxf,"%s_vx.%s.%d",SEIS_FILE,file_ext,MYID);
		sprintf(vyf,"%s_vy.%s.%d",SEIS_FILE,file_ext,MYID);
		sprintf(vzf,"%s_vz.%s.%d",SEIS_FILE,file_ext,MYID);
		sprintf(curlf,"%s_curl.%s.%d",SEIS_FILE,file_ext,MYID);
		sprintf(divf,"%s_div.%s.%d",SEIS_FILE,file_ext,MYID);
		sprintf(pf,"%s_p.%s.%d",SEIS_FILE,file_ext,MYID);

	}


	switch (SEISMO){
	case 1 : /* particle velocities only */

		fprintf(fp,"\n PE %d is writing %d seismogram traces (vx)   to  %s \n",MYID,ntr,vxf);
		outseis(fp,fopen(vxf,"w"),sectionvx,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 1);
		fprintf(fp," PE %d is writing %d seismogram traces (vy)   to  %s \n",MYID,ntr,vyf);
		outseis(fp,fopen(vyf,"w"),sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 2);
		fprintf(fp," PE %d is writing %d seismogram traces (vz)   to  %s \n",MYID,ntr,vzf);
		outseis(fp,fopen(vzf,"w"),sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 3);

		break;
	case 2 : /* pressure only */

		fprintf(fp," PE %d is writing %d seismogram traces (p)    to  %s \n",MYID,ntr,pf);
		outseis(fp,fopen(pf,"w"),sectionp, recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 0);

		break;
	case 3 : /* curl and div only */

		fprintf(fp," PE %d is writing %d seismogram traces (div)  to  %s \n",MYID,ntr,divf);
		outseis(fp,fopen(divf,"w"),sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 0);
		fprintf(fp," PE %d is writing %d seismogram traces (curl) to  %s \n\n",MYID,ntr,curlf);
		outseis(fp,fopen(curlf,"w"),sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 0);

		break;	
	case 4 : /* everything */

		fprintf(fp,"\n PE %d is writing %d seismogram traces (vx)   to  %s \n",MYID,ntr,vxf);
		outseis(fp,fopen(vxf,"w"),sectionvx,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 1);
		fprintf(fp," PE %d is writing %d seismogram traces (vy)   to  %s \n",MYID,ntr,vyf);
		outseis(fp,fopen(vyf,"w"),sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 2);
		fprintf(fp," PE %d is writing %d seismogram traces (vz)   to  %s \n",MYID,ntr,vzf);
		outseis(fp,fopen(vzf,"w"),sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 3);

		fprintf(fp," PE %d is writing %d seismogram traces (p)    to  %s \n",MYID,ntr,pf);
		outseis(fp,fopen(pf,"w"),sectionp, recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 0);

		fprintf(fp," PE %d is writing %d seismogram traces (div)  to  %s \n",MYID,ntr,divf);
		outseis(fp,fopen(divf,"w"),sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 0);
		fprintf(fp," PE %d is writing %d seismogram traces (curl) to  %s \n\n",MYID,ntr,curlf);
		outseis(fp,fopen(curlf,"w"),sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT, ishot, 0);
		break;

	}
}
