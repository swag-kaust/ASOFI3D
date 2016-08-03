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
/* ------------------------------------------------------------------------
 * This is function receiver.
 * Purpose: Find global grid positions for the receivers.
 *
------------------------------------------------------------------------*/

#include "fd.h"
#include <stdbool.h>

int **receiver(FILE *fp, int *ntr){

	/* declaration of extern variables */
	extern char REC_FILE[STRING_SIZE];
	extern float XREC1, YREC1, ZREC1, XREC2, YREC2, ZREC2;
	extern float DX, DY, DZ, REFREC[4], REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern int READREC, NGEOPH, NXG, NZG, REC_ARRAY, BOUNDARY;
	extern int MYID, DRX, DRZ, FW;

	int **recpos1, **recpos=NULL, nxrec=0, nyrec=0, nzrec=0;
	int itr=1, itr1=0, itr2=0, recflag=0, i, j, k, ifw, n;
	int nxrec1, nxrec2, nyrec1, nyrec2, nzrec1, nzrec2;
	float xrec, yrec, zrec;
	/*char rec_file_sub[STRING_SIZE]; */ /* variable not in use*/
	char bufferstring[10], buffer[STRING_SIZE];
	bool testbuff1, testbuff2, testbuff3;
	//bool testbuff4;
	FILE *fpr;


	if (MYID==0){
		if (READREC){ /* read receiver positions from file */
			fprintf(fp,"\n Reading receiver positions from file: \n\t%s\n",REC_FILE);
			
			fpr=fopen(REC_FILE,"r");
			
			if (fpr==NULL) err(" Receiver file could not be opened !");
			*ntr=0;
	 	 	
	 	 	/* counts the number of receivers in the receiver file */
			while(fgets(buffer, STRING_SIZE, fpr)){
					testbuff1=strchr(buffer,'#');
				testbuff2=strchr(buffer,'%');
				testbuff3=sscanf(buffer,"%s",bufferstring)==1;

				/*the following output is for debugging*/
				//testbuff4=(testbuff1==1 || testbuff2==1);
				//fprintf(fp," buffer: _%s_with testbuff1=_%i_ testbuff2=_%i_testbuff3=_%i_ testbuff4=_%i_\n",buffer,testbuff1, testbuff2, testbuff3,testbuff4);
				/* checks if the line contains a '%' or '#' character which indicates a 
					comment line, and if the reading of a string was successful, 
					which is not the case for an empty line*/
				if (((testbuff1==1 || testbuff2==1)==0) && testbuff3==1) ++(*ntr);
			}
			rewind(fpr);
		
			recpos1=imatrix(1,4,1,*ntr);
			for (itr=1;itr<=*ntr;itr++){
				fscanf(fpr,"%f%f%f\n",&xrec, &yrec, &zrec);
				/*note that "y" denotes the vertical coordinate*/
				recpos1[1][itr]=iround((xrec+REFREC[1])/DX);
				recpos1[2][itr]=iround((yrec+REFREC[2])/DY);
				recpos1[3][itr]=iround((zrec+REFREC[3])/DZ);
				recpos1[4][itr]=itr;
			}
			fclose(fpr);
			fprintf(fp," Message from function receiver (written by PE %d):\n",MYID);
			fprintf(fp," Number of receiver positions found: %i\n",*ntr);

			/* check if more than one receiver is located
								         at the same gridpoint */
			for (itr=1;itr<=(*ntr-1);itr++)
				for (itr1=itr+1;itr1<=*ntr;itr1++)
					if ((recpos1[1][itr]==recpos1[1][itr1])
					    && (recpos1[2][itr]==recpos1[2][itr1])
					    && (recpos1[3][itr]==recpos1[3][itr1]))
						recpos1[1][itr1]=-(++recflag);
			recpos=imatrix(1,4,1,*ntr-recflag);
			for (itr=1;itr<=*ntr;itr++)
				if (recpos1[1][itr]>0){
					recpos[1][++itr2]=recpos1[1][itr];
					recpos[2][itr2]=recpos1[2][itr];
					recpos[3][itr2]=recpos1[3][itr];
					recpos[4][itr2]=itr2;
				}
			*ntr=itr2;
			if ((recflag>0)||(itr2<(itr-1))){
				fprintf(fp,"\n\n");
				fprintf(fp," Warning:\n");
				fprintf(fp," Several receivers located at the same gridpoint !\n");
				fprintf(fp," Number of receivers reduced to %i\n", *ntr);
				fprintf(fp,"\n\n");
			}
			free_imatrix(recpos1,1,4,1,*ntr);

		}
		else if (REC_ARRAY>0){
			/* receiver array in the horizonzal X-Z plane */
			ifw=FW;  /* frame width in gridpoints */
			if (BOUNDARY==1) ifw=0;
			*ntr=((1+(NZG-2*ifw)/DRZ)*(1+(NXG-2*ifw)/DRX))*REC_ARRAY;
			recpos=imatrix(1,4,1,*ntr);
			itr=0;
			for (n=0;n<=REC_ARRAY-1;n++){
				j=iround((REC_ARRAY_DEPTH+REC_ARRAY_DIST*(float)n)/DY);
				for (k=ifw;k<=NZG-ifw;k+=DRZ)
				for (i=ifw;i<=NXG-ifw;i+=DRX){
						itr++;
						recpos[1][itr]=i;
						recpos[2][itr]=j;
						recpos[3][itr]=k;
						recpos[4][itr]=itr;
				}
			}
		}
		else{         /* straight horizontal or vertical line of receivers */
			if ((XREC1>XREC2) || ((YREC1>YREC2) ||(ZREC1>ZREC2))){
				fprintf(fp," Coordinates of first receiver specified in input file :\n");				  
				fprintf(fp,"    %5.2f (x) , %5.2f (y) , %5.2f (z) :\n", XREC1,YREC1,ZREC1);
				fprintf(fp," Coordinates of last receiver specified in input file :\n");
				fprintf(fp,"    %5.2f (x) , %5.2f (y) , %5.2f (z) :\n", XREC2,YREC2,ZREC2);
				err("\n\n Receiver coordinates of first receiver should be equal/smaller than last receiver coordinates!");
			}
			
			nxrec1=iround(XREC1/DX); /* (nxrec1,nyrec1,nzrec1) and (nxrec2,nyrec2,nzrec2) */
			nxrec2=iround(XREC2/DX); /* are the positions of the first and last receiver*/
			nyrec1=iround(YREC1/DY); /* in gridpoints */
			nyrec2=iround(YREC2/DY);
			nzrec1=iround(ZREC1/DZ);
			nzrec2=iround(ZREC2/DZ);
			if ((abs(nyrec2-nyrec1)<=abs(nxrec2-nxrec1))||
			    (abs(nyrec2-nyrec1)<=abs(nzrec2-nzrec1))){
				if (abs(nzrec2-nzrec1)<=abs(nxrec2-nxrec1)){
					/* geophone-array horizontal x-dirextion */
					*ntr=iround((nxrec2-nxrec1)/NGEOPH)+1;
					recpos=imatrix(1,4,1,*ntr);
					for (nxrec=nxrec1;nxrec<=nxrec2;nxrec+=NGEOPH){
						nyrec=nyrec1+((nyrec2-nyrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
						nzrec=nzrec1+((nzrec2-nzrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
						itr=iround((nxrec-nxrec1)/NGEOPH)+1;
						recpos[1][itr]=nxrec;
						recpos[2][itr]=nyrec;
						recpos[3][itr]=nzrec;
						recpos[4][itr]=itr;
					}
				}
				else{   /* geophone-array horizontal z-direction */
					*ntr=iround((nzrec2-nzrec1)/NGEOPH)+1;
					recpos=imatrix(1,4,1,*ntr);
					for (nzrec=nzrec1;nzrec<=nzrec2;nzrec+=NGEOPH){
						nyrec=nyrec1+((nyrec2-nyrec1)/(nzrec2-nzrec1)*(nzrec-nzrec1));
						nxrec=nxrec1+((nxrec2-nxrec1)/(nzrec2-nzrec1)*(nzrec-nzrec1));
						itr=iround((nzrec-nzrec1)/NGEOPH)+1;
						recpos[1][itr]=nxrec;
						recpos[2][itr]=nyrec;
						recpos[3][itr]=nzrec;
						recpos[4][itr]=itr;
					}
				}
			}
			else{         /* receiver-line vertical y-direction*/
				*ntr=iround((nyrec2-nyrec1)/NGEOPH)+1;
				recpos=imatrix(1,4,1,*ntr);
				for (nyrec=nyrec1;nyrec<=nyrec2;nyrec+=NGEOPH){
					nxrec=nxrec1+((nxrec2-nxrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
					nzrec=nzrec1+((nzrec2-nzrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
					itr=iround((nyrec-nyrec1)/NGEOPH)+1;
					recpos[1][itr]=nxrec;
					recpos[2][itr]=nyrec;
					recpos[3][itr]=nzrec;
					recpos[4][itr]=itr;
				}
			}

		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(ntr,1,MPI_INT,0,MPI_COMM_WORLD);
	if (MYID!=0) recpos=imatrix(1,4,1,*ntr);
	MPI_Bcast(&recpos[1][1],(*ntr)*4,MPI_INT,0,MPI_COMM_WORLD);

	if (MYID==0){
		fprintf(fp,"\n **Message from function receiver (written by PE %d):\n",MYID);
		fprintf(fp," Number of receiver positions found: %i\n",*ntr);
		if (*ntr>25) fprintf(fp," List of receiver positions is truncated to the first 25 entries only! \n");
		fprintf(fp," If receiver positions are not exactly placed at a grid point, they are shifted to the nearest grid point.\n");
		fprintf(fp," Receiver positions (in m) in the global model-system:\n");
		fprintf(fp," x  \t\ty \t\tz \n");
		fprintf(fp," -  \t\t- \t\t- \n");
		if (*ntr>25) {
			for (k=1;k<=25;k++)
				fprintf(fp," %5.2f \t %5.2f \t %5.2f\n",recpos[1][k]*DX,recpos[2][k]*DY,recpos[3][k]*DZ);
			fprintf(fp,"\n\n");
		}
		else {
			for (k=1;k<=*ntr;k++)
				fprintf(fp," %5.2f \t %5.2f \t %5.2f\n",recpos[1][k]*DX,recpos[2][k]*DY,recpos[3][k]*DZ);
			fprintf(fp,"\n\n");
		}

	}


	return recpos;
}
