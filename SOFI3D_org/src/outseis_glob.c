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
 *   Write seismograms merged from each PE collectively to disk
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "segy.h"

/* ****************************  UNDER CONSTRUCTION !!!  ***************************** */


void  outseis_glob(FILE *fp, FILE *fpdata, float **section,
		int **recpos, int **recpos_loc, int ntr, float ** srcpos,
		int nsrc, int ns, int seis_form[6], int ishot, int comp){

	/* declaration of extern variables */
	extern int NDT,NDTSHIFT, SOFI3DVERS;
	extern float  DX, DY, DZ, TIME, DT, REFREC[4];
	extern char * FILEINP[STRING_SIZE];

	extern int ASCIIEBCDIC,LITTLEBIG,IEEEIBM;

	/* declaration of extern functions from rwsegy.c */

	extern int DDN_wsegytxth(FILE * outstream, int asciiebcdic, char * kindofdata, char * infilename, int ns, float dt, int ndt);
	extern int DDN_wsegybinh(FILE * outstream, int lbendian, int ieeeibm, int meterfeet, int ns, 
			float dt, int ndt, int ntrpr, int nart);
	extern int DDN_wbindata(FILE * outstream, int outlen, float * indata, int inlen, int first, 
			int step, int padding, int lbendian, int ieeeibm, int meterfeet);
	extern int DDN_wtxtdata(FILE * outstream, int outlen, float * indata, int inlen, int first, int step, int padding, int asciiebcdic, int meterfeet, char * dataformat, char * dataseperator, char * dataendmark);
	extern int DDN_wtraceh(FILE * outstream, int lbendian, int ieeeibm, int meterfeet, int susegy, 
			int ns, float dt, float dtshift, float srctime, int traceno, int globrecno, int shotno,
			int recno, int comp, int trid, int ntraces, int **recpos, int **recpos_loc, float ** srcpos,
			int nsrc, float xcoo, float ycoo);

	/* declaration of local variables */
	int i,j, * pint;
	segy tr;
	int tracl ;
	float xr, yr, zr, y, z, scalefac, tfloat;
	//float x;
	float XS=0.0, YS=0.0, ZS=0.0;
	const float scale=3.0;
	char kindofdata[STRING_SIZE];
	char * infilename;
	float indatap[ns];

	infilename=(char *)FILEINP;

	strncpy(kindofdata,"synthetic seismograms modeled by SOFI3D",STRING_SIZE);
	if (SOFI3DVERS==33)
		strncat(kindofdata," (3D isotropic elastic)",STRING_SIZE-sizeof("synthetic seismograms modeled by SOFI3D"));
	else if (SOFI3DVERS==32)
		strncat(kindofdata," (3D isotropic acoustic)",STRING_SIZE-sizeof("synthetic seismograms modeled by SOFI3D"));
	kindofdata[STRING_SIZE-1]='\0';


	/* if there are more than one source position/coordinate specified
	 * in SOURCE_FILE, only the first position is written into trace header fields,
	 * in case of RUN_MULTIPLE_SHOTS is activated the individual shot position is used */
	XS=srcpos[1][ishot];
	YS=srcpos[2][ishot];
	ZS=srcpos[3][ishot];

	scalefac=pow(10.0,scale);
	switch(seis_form[0]){
	case 0 :
	case 4 :
	case 5 :
		/* fprintf(stderr,"BEEN HERE !!!\n")  */
		DDN_wsegytxth(fpdata, seis_form[1], kindofdata, infilename,ns,DT,NDT);
		DDN_wsegybinh(fpdata, seis_form[2], seis_form[3], seis_form[4], ns, DT, NDT, ntr, 0);
		for(tracl=1;tracl<=ntr;tracl++){
			DDN_wtraceh(fpdata, seis_form[2], seis_form[3], seis_form[4], 1, ns, DT*NDT, DT*NDTSHIFT, 0.0, tracl, recpos_loc[4][tracl], nsrc, 
					tracl, comp, 1, ntr, recpos, recpos_loc, srcpos, nsrc, 0, 0);
			for(j=0;j<ns;j++) indatap[j]=section[tracl][j+1];
			DDN_wbindata(fpdata, ns, indatap, ns, 0, 1, 0, seis_form[2], seis_form[3], seis_form[4]);
		}
		break;
	case 1 : /* SU caution: allows IBM floats!*/ 
		for(tracl=1;tracl<=ntr;tracl++){
			DDN_wtraceh(fpdata, seis_form[2], seis_form[3], seis_form[4], 0, ns, DT*NDT, DT*NDTSHIFT, 0.0, tracl, recpos_loc[4][tracl], nsrc, 
					tracl, comp, 1, ntr, recpos, recpos_loc, srcpos, nsrc, 0, 0);
			for(j=0;j<ns;j++) indatap[j]=section[tracl][j+1];
			DDN_wbindata(fpdata, ns, indatap, ns, 0, 1, 0, seis_form[2], seis_form[3], seis_form[4]);
		}
		break;
	case 7 : /* SU ~ (IEEE) SEGY without file-headers and slightly modified trace headers (original version) */
		for(tracl=1;tracl<=ntr;tracl++){ 
			xr=recpos[1][recpos_loc[4][tracl]]*DX;
			yr=recpos[2][recpos_loc[4][tracl]]*DY;
			zr=recpos[3][recpos_loc[4][tracl]]*DZ;
			//x=xr-REFREC[1]; /* not in use*/
			y=yr-REFREC[2];
			z=zr-REFREC[3];
			tr.tracl=recpos_loc[4][tracl];      /* trace sequence number within line */
			tr.ep=comp;
			tr.cdp=recpos_loc[4][tracl];
			tr.trid=1;           /* trace identification code: 1=seismic*/
			tr.offset=iround(sqrt((XS-xr)*(XS-xr)
					+(YS-yr)*(YS-yr)
					+(ZS-zr)*(ZS-zr))*scalefac);
			tr.gelev=iround(yr*scalefac);
			tr.sdepth=iround(YS*scalefac);   /* source depth (positive) */

			/* angle between receiver position and reference point
                           (sperical coordinate system: used for tunnel geometry) */
			tr.gdel=iround(atan2(-y,z)*180.0*scalefac/PI);   
			tr.gwdep=iround(sqrt(z*z+y*y)*scalefac);

			tr.scalel=(short)(-scale);
			tr.scalco=(short)(-scale);
			tr.sx=iround(XS*scalefac);  /* X source coordinate */
			tr.sy=iround(ZS*scalefac);  /* Z source coordinate */

			/* group coordinates */
			tr.gx=iround(xr*scalefac);
			tr.gy=iround(zr*scalefac); 


			tr.ns=(unsigned short)ns; /* number of samples in this trace */
			tr.dt=(unsigned short)iround(((float)NDT*DT)*1.0e6); /* sample interval in micro-seconds */
			tr.d1=(float)(TIME/ns);        /* sample spacing for non-seismic data */
			tr.f1=0.0;              /* first sample location for non-seismic data */
			tr.d2=0.0;        /* sample spacing between traces */
			tr.f2=0.0;

			fwrite(&tr,240,1,fpdata);
			if (seis_form[3]) for(j=0;j<ns;j++) tr.data[j]=section[tracl][j+1]/0.3048; /* FEET */
			else for(j=0;j<ns;j++) tr.data[j]=section[tracl][j+1];
			fwrite(tr.data,4,ns,fpdata);
		}
		break;
	case 6 : /* pseudo-SU (segy-traces headers & traces) rwsegy.c */
		for(tracl=1;tracl<=ntr;tracl++){
			DDN_wtraceh(fpdata, seis_form[2], seis_form[3], seis_form[4], 1, ns, DT*NDT, DT*NDTSHIFT, 0.0, tracl, recpos_loc[4][tracl], nsrc, 
					tracl, comp, 1, ntr, recpos, recpos_loc, srcpos, nsrc, 0, 0);
			for(j=0;j<ns;j++) indatap[j]=section[tracl][j+1];
			DDN_wbindata(fpdata, ns, indatap, ns, 0, 1, 0, seis_form[2], seis_form[3], seis_form[4]);
		}
		break;
	case 2 :if (ASCIIEBCDIC==seis_form[1]){
		if (seis_form[4]==1) /*OUTPUT IN FEET*/ switch(seis_form[3]){
		case 1:	for(j=1;j<=ns;j++){         /*ASCII ONE COLUMN PER TRACE */
			for(i=1;i<=ntr;i++) fprintf(fpdata,"%e\t", section[i][j]/0.3048);
			fprintf(fpdata,"\n");
		}
		break;
		case 2: for(i=1;i<=ntr;i++) /*ASCII ONE LINE*/ for(j=1;j<=ns;j++) fprintf(fpdata,"%e ",section[i][j]/0.3048);
		break;
		case 3: for(i=1;i<=ntr;i++){         /*ASCII ONE LINE PER RECEIVER */
			for(j=1;j<=ns;j++) fprintf(fpdata,"%e\t", section[i][j]/0.3048);
			fprintf(fpdata,"\n");
		}
		break;
		default:for(i=1;i<=ntr;i++) /*ASCII ONE COLUMN*/ for(j=1;j<=ns;j++) fprintf(fpdata,"%e\n",section[i][j]/0.3048);
		}
		else /*OUTPUT IN METER*/ switch(seis_form[3]){
		case 1:	for(j=1;j<=ns;j++){         /*ASCII ONE COLUMN PER TRACE */
			for(i=1;i<=ntr;i++) fprintf(fpdata,"%e\t", section[i][j]);
			fprintf(fpdata,"\n");
		}
		break;
		case 2: for(i=1;i<=ntr;i++) /*ASCII ONE LINE*/ for(j=1;j<=ns;j++) fprintf(fpdata,"%e ",section[i][j]);
		break;
		case 3: for(i=1;i<=ntr;i++){         /*ASCII ONE LINE PER RECEIVER */
			for(j=1;j<=ns;j++) fprintf(fpdata,"%e\t", section[i][j]);
			fprintf(fpdata,"\n");
		}
		break;
		default:for(i=1;i<=ntr;i++) /*ASCII ONE COLUMN*/ for(j=1;j<=ns;j++) fprintf(fpdata,"%e\n",section[i][j]);
		}
	}
	else{
		switch(seis_form[3]){
		case 1:for(j=1;j<=ns;j++){         /*ASCII ONE COLUMN PER TRACE */
			for(i=1;i<=ntr;i++) fprintf(fpdata,"%e\t", section[i][j]);
			fprintf(fpdata,"\n");
		}
		break;
		case 2:
		case 3:
		default:

			/* dseperator[]="\t", ddataendmark[]="\n", ddataformat[]="%e"

		DDN_wtxtdata(fpdata,outlen,&indata,inlen,0,1,0,seis_form[1],seis_form[4], char * dataformat, char * dataseperator, char * dataendmark); */

			;
		}
		break;

		case 3 :                             /*BINARY */
			if ((seis_form[3]==IEEEIBM)&&(seis_form[2]==LITTLEBIG)){ /* OUTPUT NATIVE FLOATS */
				if (seis_form[4]==1){ /*OUTPUT IN FEET*/
					for(i=1;i<=ntr;i++) for(j=1;j<=ns;j++){
						tfloat=section[i][j]/0.3048;
						fwrite(&tfloat,sizeof(float),1,fpdata);
					}
				}
				else for(i=1;i<=ntr;i++) for(j=1;j<=ns;j++) fwrite(&section[i][j],sizeof(float),1,fpdata);
			}
			else {
				if ((seis_form[3]==IEEEIBM)&&(seis_form[2]!=LITTLEBIG)) /* SWAP FLOATS */ {
					if (seis_form[4]==1) /*OUTPUT IN FEET*/ for(i=1;i<=ntr;i++) for(j=1;j<=ns;j++){
						tfloat=section[i][j]/0.3048;
						pint=(int *) &tfloat;
						*pint=((*pint>>24)&0xff)|((*pint&0xff)<<24)|((*pint>>8)&0xff00)|((*pint&0xff00)<<8);
						fwrite(&tfloat,sizeof(float),1,fpdata);
					}
					else for(i=1;i<=ntr;i++) for(j=1;j<=ns;j++){
						tfloat=section[i][j];
						pint=(int *) &tfloat;
						*pint=((*pint>>24)&0xff)|((*pint&0xff)<<24)|((*pint>>8)&0xff00)|((*pint&0xff00)<<8);
						fwrite(&tfloat,sizeof(float),1,fpdata);
					}
				}
				else
					for(i=1;i<=ntr;i++) for(j=1;j<=ns;j++)		DDN_wbindata(fpdata,1,&section[i][j],1,0,1,0,seis_form[2],seis_form[3],seis_form[4]); /* probably extremely slow !!! */
			}	
	}
	break;

	default :
		fprintf(fp," Don't know data format for seismograms !\n");
		fprintf(fp," No output written. ");
	}

	fclose(fpdata);
}
