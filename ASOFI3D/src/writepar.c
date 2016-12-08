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
/* ----------------------------------------------------------------------
 * This is function write_par.
   Purpose: Writing FD-Parameters to stdout or log-file                           
   
   ----------------------------------------------------------------------*/
#include <limits.h>
#include "fd.h"

/* printing all important parameters to FILE *fp */
void writepar(FILE *fp, int ns){

	/* declaration of extern variables */
	extern int   NX, NY, NZ, NT, SOURCE_SHAPE, SOURCE_TYPE, FDORDER,FDORDER_TIME, RUN_MULTIPLE_SHOTS, NDT;
	extern int  SNAP, SNAP_FORMAT, REC_ARRAY, L, SNAP_PLANE,FW;
	extern float DX, DY, DZ, TIME, DT, TS, *FL, TAU, PLANE_WAVE_DEPTH;
	extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2;
	extern float SOURCE_ALPHA, SOURCE_BETA, AMON, STR, DIP, RAKE;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern int SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT[6], FREE_SURF;
	extern int  READMOD, READREC, DRX, DRZ, BOUNDARY, SRCREC, IDX, IDY, IDZ;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4], DAMPING;
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], REC_FILE[STRING_SIZE], SEIS_FILE[STRING_SIZE];
	extern char  MFILE[STRING_SIZE];
	extern int NP, NPROCX, NPROCY, NPROCZ, MYID;
	
	/* definition of local variables */
	char th1[3], file_ext[8];
	char th2[3];	
	int l;

	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ********* PARAMETERS AS SPECIFIED IN INPUT FILE **********");
	fprintf(fp,"\n **********************************************************\n\n");
	
	/*note that "y" denotes the vertical coordinate*/
	fprintf(fp,"\n\n **Message from write_par (printed by PE %d):\n\n",MYID);
	fprintf(fp,"------------------------- Processors ------------------------\n");
	fprintf(fp," Number of PEs in horizontal x-direction (NPROCX): %d\n",NPROCX);
	fprintf(fp," Number of PEs in horizontal y-direction (NPROCY): %d\n",NPROCY);
	fprintf(fp," Number of PEs in vertical   z-direction (NPROCZ): %d\n",NPROCZ);
	fprintf(fp," Total number of PEs in use: %d\n",NP);
	fprintf(fp,"\n");
	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", NX);
	fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", NY);
	fprintf(fp," Number of gridpoints in z-direction (NZ): %i\n", NZ);
	fprintf(fp," Grid-spacing in x-direction (DX): %5.4f meter\n", DX);
	fprintf(fp," Grid-spacing in y-direction (DY): %5.4f meter\n", DY);
	fprintf(fp," Grid-spacing in z-direction (DZ): %5.4f meter\n", DZ);
	fprintf(fp," Time of wave propagation (T): %5.4f seconds\n",TIME);
	fprintf(fp," Sampling rate of the seismogram output (NDT): %i \n",NDT);
	fprintf(fp," Timestep (DT): %5.4e seconds\n", DT);
	
	fprintf(fp," Number of timesteps: %i \n",NT);
	fprintf(fp,"\n");


	fprintf(fp," ------------------------- ORDER OF FD OPERATORS --------------\n");
	fprintf(fp," Order of spatial FD operators: %i \n",FDORDER);
	if ((FDORDER<0)||(FDORDER%2!=0)||(FDORDER>12))
		err(" Incorrect FDORDER (must be 2, 4, 8, or 12) ! ");
	fprintf(fp," Order of temporal FD operators: %i \n",FDORDER_TIME);
	
	
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");

	if ((SRCREC) && (!PLANE_WAVE_DEPTH)){
		fprintf(fp," Reading source positions, time delay, centre frequency \n");
		fprintf(fp," and initial amplitude from ASCII-file \n");
		fprintf(fp,"\t%s\n\n",SOURCE_FILE);
		if (RUN_MULTIPLE_SHOTS)	fprintf(fp,"\n SOFI3D will run (independent) simulations for each source defined in %s\n\n", SOURCE_FILE);
 
	} else {
		fprintf(fp," Plane wave excitation: depth= %5.2f meter \n",PLANE_WAVE_DEPTH);
 		fprintf(fp," duration of source signal: %e seconds\n",TS);
 		fprintf(fp," (centre frequency is approximately %e Hz)\n",1.0/TS);
	}


	
	fprintf(fp," Wavelet of source:");

	switch (SOURCE_SHAPE){
	case 1 :
		fprintf(fp," Ricker\n");
		break;
	case 2 :
		fprintf(fp," Fuchs-Mueller\n");
		break;
	case 3 :
		fprintf(fp," reading from \n\t %s\n",SIGNAL_FILE);
		break;
	case 4 :
		fprintf(fp," sinus raised to the power of 3.0 \n");
		break;
	
	default :
		err(" Sorry, incorrect specification of source wavelet ! ");
	}

	fprintf(fp," Default type of source:");
	switch (SOURCE_TYPE){
	case 1 :
		fprintf(fp," explosive point source (concentrated at a single gridpoint)\n");
		break;
	case 2 :
		fprintf(fp," point source with directive force in x-direction\n");
		break;
	case 3 :
		fprintf(fp," point source with directive force in (vertical) y-direction\n");
		break;
	case 4 :
		fprintf(fp," point source with directive force in  z-direction\n");
		break;
	case 5 :
		fprintf(fp," point source with directive force in  custom-direction\n");
		fprintf(fp," Angle between x and y (vertical) direticon (SOURCE_ALPHA): %f\n", SOURCE_ALPHA);
		fprintf(fp," Angle between x and z direticon (SOURCE_BETA): %f\n", SOURCE_BETA);
		break;
	case 6 :
		fprintf(fp," Seismic moment tensor point source with double couple mechanism \n");
		fprintf(fp," Seismic moment of source (AMON): %f\n", AMON);
		fprintf(fp," Strike angle of fault plane (STR): %f\n", STR);
		fprintf(fp," Dip angle of fault plane (DIP): %f\n", DIP);
		fprintf(fp," Rake angle of slip vector (RAKE): %f\n", RAKE);
		
		break;
	default :
		fprintf(fp," WARNING: Default type of source ('%d') not available -> changed to explosive! ", SOURCE_TYPE);
		SOURCE_TYPE=1;
		break;
	}
	
	fprintf(fp,"\n");

	if (SEISMO){
		fprintf(fp," ------------------------- RECEIVER  ------- -------------------\n");
		if (READREC){
			fprintf(fp," Reading receiver positions from file \n");
			fprintf(fp,"\t%s\n\n",REC_FILE);
			fprintf(fp," reference_point_for_receiver_coordinate_system:\n");
			fprintf(fp," x=%f \ty=%f\t z=%f\n",REFREC[1], REFREC[2], REFREC[3]);

		} else if (REC_ARRAY>0){
				fprintf(fp," Horitontal plane of receivers.\n");
				fprintf(fp," Number of planes: %d \n",REC_ARRAY);
				fprintf(fp," Depth of upper plane: %e m \n",REC_ARRAY_DEPTH);
				fprintf(fp," Vertical increment between planes: %e m \n",REC_ARRAY_DIST);
				fprintf(fp," Distance between receivers in x-direction within plane: %i \n", DRX);		
				fprintf(fp," Distance between receivers in y-direction within plane: %i \n", DRZ);
			}			
			else{
			fprintf(fp," Receiver line: \n");			
			fprintf(fp," First receiver position (XREC1,YREC1,ZREC1) = (%5.3f, %5.3f, %5.3f) m\n",
			    XREC1,YREC1,ZREC1);
			fprintf(fp," Last receiver position (XREC2,YREC2,ZREC2)  = (%5.3f, %5.3f, %5.3f) m\n",
			    XREC2,YREC2,ZREC2);
			fprintf(fp,"\n");
		}
	}
	fprintf(fp," ------------------------- FREE SURFACE ------------------------\n");
	if (FREE_SURF) fprintf(fp," There is a free surface at the top of the model ! \n");
	else fprintf(fp," There is no free surface at the top of the model ! \n");
	fprintf(fp,"\n");

	fprintf(fp," ------------------------- ABSORBING FRAME ---------------------\n");
	if (FW>0){
		fprintf(fp," Width of absorbing frame is %i grid points.\n",FW);
		fprintf(fp," The percentage of amplitude decay at the edge is set to %f \n",DAMPING);
	}
	else {
		fprintf(fp," Absorbing frame not installed ! \n");
		fprintf(fp," Be aware of artificial reflections from the edges of the numerical mesh ! \n");
	}

	switch (BOUNDARY){
		case 0 :
			fprintf(fp," No periodic boundary condition.\n");
			break;
		case 1 :
			fprintf(fp," Periodic boundary condition at left/right and front/back sides of global grid.\n");
			break;
		default :
			warning(" Wrong integer value for BOUNDARY specified in parameter file! ");
			warning(" No periodic boundary condition will be applied ");
			BOUNDARY=0;
			break;
	}


	
	if (READMOD){
		fprintf(fp," ------------------------- MODEL-FILES -------------------------\n");
		fprintf(fp," Names of model-files: \n");
		fprintf(fp,"\t shear wave velocities:\n\t %s.vs\n",MFILE);
		fprintf(fp,"\t tau for shear waves:\n\t %s.ts\n",MFILE);
		fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
		fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
		fprintf(fp,"\t tau for P-waves:\n\t %s.tp\n",MFILE);
		for (l=1;l<=L;l++) fprintf(fp,"\t %1i. relaxation frequencies: %s.f%1i\n",l,MFILE,l);
	}

	fprintf(fp,"\n");
	fprintf(fp," ------------------------- Q-APROXIMATION --------------------\n");
	fprintf(fp," Number of relaxation mechanisms (L): %i\n",L);
	fprintf(fp," The L relaxation frequencies are at:  \n");
	for (l=1;l<=L;l++) fprintf(fp,"\t%f",FL[l]);
	fprintf(fp," Hz\n");
	fprintf(fp," Value for tau is : %f\n",TAU);


	if (SNAP){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of");
		switch(SNAP){
		case 1:
			fprintf(fp," particle velocity.\n");
			break;
		case 2:
			fprintf(fp," pressure field.\n");
			break;
		case 3:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			break;
		case 4:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			fprintf(fp," and particle velocity.\n");
			break;
		default:
			err(" sorry, incorrect value for SNAP ! \n");
		}

		fprintf(fp," \t first (TSNAP1)= %8.5f s\n", TSNAP1);
		fprintf(fp," \t last (TSNAP2)=%8.5f s\n",TSNAP2);
		fprintf(fp," \t increment (TSNAPINC) =%8.5f s\n\n",TSNAPINC);
		fprintf(fp," \t spacing in x-direction (IDX*DX) =%8.5f m\n",IDX*DX);
		fprintf(fp," \t spacing in y-direction (IDY*DY) =%8.5f m\n",IDY*DY);
		fprintf(fp," \t spacing in z-direction (IDZ*DZ) =%8.5f m\n",IDZ*DZ);
		fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);
		switch (SNAP_FORMAT){
		case 1 :
			err(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
			break;
		default:
			err(" Don't know the format for the Snapshot-data ! \n");
		}
		switch (SNAP_PLANE){
		case 1 :
			fprintf(fp," \nDiv and curl output will be as Energy without sign. \n");
			break;
		case 2 :
			fprintf(fp," \nDiv and curl output will be as Energy with sign true for xz-plane. \n");
			break;
		case 3 :
			fprintf(fp," \nDiv and curl output will be as Energy with sign true for xy-plane. \n");
			break;
		case 4 :
			fprintf(fp," \nDiv and curl output will be as Energy with sign true for yz-plane. \n");
			break;
		}

		fprintf(fp,"\n\n");
	}
	if (SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
		switch (SEIS_FORMAT[0]){
			case 0: sprintf(file_ext,"sgy"); break;
			case 1: sprintf(file_ext,"su");  break;
			case 2: sprintf(file_ext,"txt"); break;
			case 3: sprintf(file_ext,"bin"); break;
			case 4: sprintf(file_ext,"sgy"); break;
			case 5: sprintf(file_ext,"sgy"); break;
		}
		
		if ((SEISMO==1) || (SEISMO==4)){
			fprintf(fp," Seismograms of ");
			fprintf(fp," x-, y-, and z-component");
			fprintf(fp," of particle velocity.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s_vx.%s\n\t%s_vy.%s\n\t%s_vz.%s\n",SEIS_FILE,file_ext,SEIS_FILE,file_ext,SEIS_FILE,file_ext);
		}
		if ((SEISMO==2) || (SEISMO==4)){
			fprintf(fp," Seismograms of pressure field (hydrophones).\n");
			fprintf(fp," output-file: \n ");
			fprintf(fp,"\t%s_p.%s\n",SEIS_FILE,file_ext);
		}
		if ((SEISMO==3) || (SEISMO==4)){
			fprintf(fp," Seismograms of curl (S-wave component) and div (P-wave component of wavefield).\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s_rot.%s \n\t%s_div.%s\n",SEIS_FILE,file_ext,SEIS_FILE,file_ext);
			
		}		
		
		if (NDT==0) {NDT=1; fprintf(fp," NDT set to %d.\n",NDT);}
		else if (NDT<0) {
			NDT=-NDT;
			fprintf(fp," Negative NDT set to its absolute value %d.\n",NDT);
		}
		if (NDTSHIFT<0){
			NDTSHIFT=-NDTSHIFT;
			fprintf(fp," Negative NDTSHIFT set to its absolute value %d.\n",NDTSHIFT);
		}
		if (ns) {
			switch (NDT) {
				case 1 : strcpy(th1,"st"); break;
				case 2 : strcpy(th1,"nd"); break;
				case 3 : strcpy(th1,"rd"); break;
				default: strcpy(th1,"th"); break;
			}
			switch (NDTSHIFT) {
				case 1 : strcpy(th2,"st"); break;
				case 2 : strcpy(th2,"nd"); break;
				case 3 : strcpy(th2,"rd"); break;
				default: strcpy(th2,"th"); break;
			}			
			fprintf(fp," Amplitudes will be written every %d%s time-step, starting at the %d%s.\n",NDT,th1,NDTSHIFT,th2);
		}
		else{
			if ((SEIS_FORMAT[0]==2)||(SEIS_FORMAT[0]==3)) 
				fprintf(fp," Warning: seismogram files will be empty! \n");
		 	else if ((SEIS_FORMAT[0]==0)||(SEIS_FORMAT[0]==1)||(SEIS_FORMAT[0]==4)) 
				fprintf(fp," Warning: seismogram files will contain only headers! \n");
		}

		switch (SEIS_FORMAT[0]){
		case 0 :
		case 5 :
			fprintf(fp," Seismograms are written in SEG-Y format. \n");
			if (!SEIS_FORMAT[1]) fprintf(fp," \t textual header: ASCII \n");
			else if (SEIS_FORMAT[1]==1) fprintf(fp," \t textual header: EBCDIC \n");
			if (!SEIS_FORMAT[2]) fprintf(fp," \t byte order: little endian \n");
			else if (SEIS_FORMAT[2]==1) fprintf(fp," \t byte order: big endian \n");
			if (!SEIS_FORMAT[3]) fprintf(fp," \t data type: IEEE 4-byte floats \n");
			else if (SEIS_FORMAT[3]==1) fprintf(fp," \t data type: IBM 4-byte floats \n");
			if (!SEIS_FORMAT[4]) fprintf(fp," \t coordinate unit: meter \n");
			else if (SEIS_FORMAT[4]==1) fprintf(fp," \t coordinate unit: feet \n");
			break;
		
		case 1 :
			fprintf(fp," Seismograms are written in SU-format. \n");
			if (!SEIS_FORMAT[2]) fprintf(fp," \t byte order: little endian \n");
			else if (SEIS_FORMAT[2]==1) fprintf(fp," \t byte order: big endian \n");
			if (!SEIS_FORMAT[3]) fprintf(fp," \t data type: IEEE 4-byte floats \n");
			else if (SEIS_FORMAT[3]==1) fprintf(fp," \t CAUTION: data type: IBM 4-byte floats \n");
			if (!SEIS_FORMAT[4]) fprintf(fp," \t coordinate unit: meter \n");
			else if (SEIS_FORMAT[4]==1) fprintf(fp," \t coordinate unit: feet \n");
			break;
		case 2 :
			if (!SEIS_FORMAT[1]) fprintf(fp," Seismograms are written in ASCII. \n");
			else if (!SEIS_FORMAT[1]==1) fprintf(fp," Seismograms are written in EBCDIC. \n");
			if (!SEIS_FORMAT[4]) fprintf(fp," \t coordinate unit: meter \n");
			else if (SEIS_FORMAT[4]==1) fprintf(fp," \t coordinate unit: feet \n");			
			break;
		case 3 :
			fprintf(fp," Seismograms are written in binary format.");
			if (!SEIS_FORMAT[2]) fprintf(fp," \t byte order: little endian \n");
			else if (SEIS_FORMAT[2]==1) fprintf(fp," \t byte order: big endian \n");
			if (!SEIS_FORMAT[3]) fprintf(fp," \t data type: IEEE 4-byte floats \n");
			else if (SEIS_FORMAT[3]==1) fprintf(fp," \t data type: IBM 4-byte floats \n");
			if (!SEIS_FORMAT[4]) fprintf(fp," \t coordinate unit: meter \n");
			else if (SEIS_FORMAT[4]==1) fprintf(fp," \t coordinate unit: feet \n");			
			break;
		case 6 :fprintf(fp," Seismograms are written in pseudo SU-format (SEG-Y with trace headers only). \n");
			if (!SEIS_FORMAT[2]) fprintf(fp," \t byte order: little endian \n");
			else if (SEIS_FORMAT[2]==1) fprintf(fp," \t byte order: big endian \n");
			if (!SEIS_FORMAT[3]) fprintf(fp," \t data type: IEEE 4-byte floats \n");
			else if (SEIS_FORMAT[3]==1) fprintf(fp," \t CAUTION: data type: IBM 4-byte floats \n");
			if (!SEIS_FORMAT[4]) fprintf(fp," \t coordinate unit: meter \n");
			else if (SEIS_FORMAT[4]==1) fprintf(fp," \t coordinate unit: feet \n");
			break;
		case 7 :fprintf(fp," Seismograms are written in SU-format (output in meter, native endian and floats). \n");
			break;	
		default:
			err(" Sorry. Unknown format for seismic data! \n");
		}
		fprintf(fp," samplingrate of seismic data:                %e s\n",NDT*DT);
		if (!READREC) fprintf(fp," Trace-spacing: %e m\n", NGEOPH*DX);
		fprintf(fp," Number of samples per receiver:              %i \n", ns);
		if (!SEIS_FORMAT[5]) SEIS_FORMAT[5]= USHRT_MAX; /* default */
		if (ns<=abs(SEIS_FORMAT[5])) fprintf(fp," Number of samples per trace:                 %i \n", ns);
		else if (SEIS_FORMAT[5]<0){
			fprintf(fp," Number of samples per trace:                 %i \n", -SEIS_FORMAT[5]);
			fprintf(fp," Number of traces per receiver:               %i \n", ns / -SEIS_FORMAT[5]);
			fprintf(fp," Number of significant samples in last trace: %i \n", ns % -SEIS_FORMAT[5]);
		}
		else {
			fprintf(fp," Maximum allowed number of samples per trace: %i \n", SEIS_FORMAT[5]);
			err(" Sorry. Too many samples per receiver! \n");
		} 
		fprintf(fp," ----------------------------------------------------------\n");
		fprintf(fp,"\n");
		fprintf(fp,"\n");
	}

	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ******* PARAMETERS READ or PROCESSED within SOFI3D *******");
	fprintf(fp,"\n **********************************************************\n\n");


}
