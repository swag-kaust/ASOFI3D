/*------------------------------------------------------------------------
 *   prepares to write seismograms from each PE individually to disk
 *------------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"

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
