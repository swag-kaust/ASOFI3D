#include <string.h>

#include "fd.h"
#include "globvar.h"


void removespace(char *str) {

    char *p1 = str, *p2 = str;
    do
	while (*p2 == ' ' || *p2 == '\t')
	    p2++;
    while ((*p1++ = *p2++));
}



void madinput(char header[], float *** DEN ){


    // -------------  header file Reading ---------------------
    extern float DX, DY, DZ;//, OX, OY, OZ;
    extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], MYID;
    extern FILE *FP;
    fprintf(FP,"\n\n\n--------------------------------------------------------- \n");
    fprintf(FP," \n \n *********** Madagascar Input Start *************** \n \n");
    fprintf(FP,"--------------------------------------------------------- \n");
    fprintf(FP,"\n \n \t \t Madagascar file : \t%s \n",header);
    fflush( FP);
    // local variables
    int i, ii, j, jj, k, kk;

    float tempRho=0.0;

    // -------------  header file Reading ---------------------
    char *pch;
    FILE *ioh_file;
    char filename[STRING_SIZE];
    char *s="\"";
    char *token;
    char binary_file[STRING_SIZE];

    char tempbuff[STRING_SIZE];


    // Read the header file
    sprintf(filename,"%s",header);
    ioh_file=fopen(filename,"r");
    if (ioh_file==NULL) err("\t \t \t :( Could not open Header :( ");




    while(!feof(ioh_file)){
	if(fgets(tempbuff,STRING_SIZE,ioh_file)!=NULL ){
	    pch = strtok (tempbuff," =");
	    removespace(tempbuff);

	    if (strcmp(pch,"in")==0)    {

		pch = strtok (NULL, " =");
		pch[strcspn(pch,"\r\n")]=0; // to remove the newline character
		token = strtok(pch,s); // remove the ""
		strcpy(binary_file,token);
	    }
	    else if (strcmp(pch,"n1")==0)       {
		pch = strtok (NULL, " = ");
		NZ=atoi(pch);
	    }
	    else if (strcmp(pch,"n2")==0)       {
		pch = strtok (NULL, " = ");
		NX=atoi(pch);
	    }
	    else if (strcmp(pch,"n3")==0)       {
		pch = strtok (NULL, " = ");
		NY=atoi(pch);
	    }
	    else if (strcmp(pch,"d1")==0)       {
		pch = strtok (NULL, " = ");
		DZ=atof(pch);
	    }

	    else if (strcmp(pch,"d2")==0)       {
		pch = strtok (NULL, " = ");
		DX=atof(pch);
	    }

	    else if (strcmp(pch,"d3")==0)       {
		pch = strtok (NULL, " = ");
		DY=atof(pch);
	    }
	    else if (strcmp(pch,"o1")==0)       {
		pch = strtok (NULL, " = ");
	    }
	    else if (strcmp(pch,"o2")==0)       {
		pch = strtok (NULL, " = ");
	    }

	    else if (strcmp(pch,"o3")==0)       {
		pch = strtok (NULL, " = ");
	    }

	}
    }

    fclose(ioh_file);

    /*
    // over-write the json parameters
    NXG=Nx;
    NZG=Nz;
    NYG=Ny;
    DZ=Dz;
    DX=Dx;
    DY=Dy;
    // ----------------------- Header informaition Reading Completed ---------------
    */

    // --------- Binary allocation to Memory ---------------------------

    fprintf(FP,"i\n \n \t MYID = \t%d \n\n",MYID);

    const int format=3;

    ioh_file=fopen(binary_file,"r");
    fprintf(FP,"\n\n \t \t Binary : \t%s\n\n",binary_file);
    if (ioh_file==NULL) err("\t \t \t :( No binaries  present :( ");


    for (k=1;k<=NZG;k++){
	for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){

		tempRho=readdsk(ioh_file, format);

		if ((POS[1]==((i-1)/NX)) 
			&& (POS[2]==((j-1)/NY)) 
			&& (POS[3]==((k-1)/NZ)))
		{
		    if (tempRho!=5000) fprintf(FP,"\n New in %g Nx %d Ny %d Nz %d Nxg %d Nyg %d Nzg %d",tempRho,NX,NY,NZ,POS[1], POS[2], POS[3]);
		    ii=i-POS[1]*NX;
		    jj=j-POS[2]*NY;
		    kk=k-POS[3]*NZ;
		    DEN[jj][ii][kk]=tempRho;
		}
	    }
	}
    }

    fclose(ioh_file);
    fprintf(FP,"\n\n\n--------------------------------------------------------- \n");
    fprintf(FP," \n \n *********** Madagascar Input Finish *************** \n \n");
    fprintf(FP,"--------------------------------------------------------- \n");
}
