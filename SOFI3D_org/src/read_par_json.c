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
 * This is function read_par_json.
 *  Purpose: Reading FD-Parameters from parameter-file (formatted according to json standard)
 *
 ------------------------------------------------------------------------*/

#include "fd.h"

void read_par_json(FILE *fp, char *fileinp){

	/* declaration of extern variables */
	extern int   NX, NY, NZ, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, SNAP_PLANE;
	extern int DRX, DRZ, L, SRCREC, FDORDER, FW, FDCOEFF, FDORDER_TIME;
	extern float DX, DY, DZ, TIME, DT, TS, *FL, TAU, FREF, PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE;
	extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2, SOURCE_ALPHA, SOURCE_BETA, AMON, STR, DIP, RAKE;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern int SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT[6], FREE_SURF, READMOD, MOD_FORMAT[6], READREC, RUN_MULTIPLE_SHOTS;
	extern int BOUNDARY, REC_ARRAY, LOG, IDX, IDY, IDZ, ABS_TYPE, WRITE_MODELFILES;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4], DAMPING, FPML, VPPML,NPOWER, K_MAX_CPML;
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE[STRING_SIZE];
	extern int NPROCX,NPROCY,NPROCZ, CHECKPTREAD, CHECKPTWRITE, OUTNTIMESTEPINFO, OUTSOURCEWAVELET;
	extern int ASCIIEBCDIC,LITTLEBIG,IEEEIBM;

	extern float FC,AMP, REFSRC[3], SRC_DT, SRCTSHIFT;
	extern int SRC_MF, SIGNAL_FORMAT[6], SRCOUT_PAR[6], FSRC, JSRC, LSRC;
	extern char SRCOUT_FILE[STRING_SIZE];

	/* definition of local variables */
	int number_readobjects=0;
	int number_defaultobjects=0;
	char tempstring[STRING_SIZE];
	char varname_tmp1[STRING_SIZE], value_tmp1[STRING_SIZE];

	char ** varname_list, ** value_list;
	char ** varnamedefault_list, ** valuedefault_list;

	//allocate first object in list
	varname_list = malloc(STRING_SIZE*sizeof(char*));
	value_list = malloc(STRING_SIZE*sizeof(char*));
	varnamedefault_list = malloc(STRING_SIZE*sizeof(char*));
	valuedefault_list = malloc(STRING_SIZE*sizeof(char*));

	//read in objects from file
	number_readobjects=read_objects_from_intputfile(fp, fileinp, varname_list, value_list);
	fprintf(fp,"\nFrom input file %s, %i objects have been read in. \n",fileinp, number_readobjects);

	//print objects to screen
	fprintf(fp, "\n===========================================================");
	fprintf(fp, "\n=   List of Parameters read by the built in Json Parser   =");
	print_objectlist_screen(fp, number_readobjects, varname_list, value_list);

	//extract variables form object list

	//if a NON-CRITICAL or NOT-IN-USE value is not in the input file, it will be set to default
	//and a list of all these objects (varnamedefault_list )will be displayed to screen
	//note that most default values are set in globvar.h, here values are displayed only!

	/*=================================
	  section general grid and discretization parameters
	  =================================*/

	if (get_int_from_objectlist("NPROCX",number_readobjects,&NPROCX,varname_list, value_list))
		err("Variable NPROCX could not be retrieved from the json input file!");
	if (get_int_from_objectlist("NPROCY",number_readobjects,&NPROCY,varname_list, value_list))
		err("Variable NPROCY could not be retrieved from the json input file!");
	if (get_int_from_objectlist("NPROCZ",number_readobjects,&NPROCZ,varname_list, value_list))
		err("Variable NPROCY could not be retrieved from the json input file!");
	/*note that "y" is used for the vertical coordinate */
	if (get_int_from_objectlist("FDORDER",number_readobjects,&FDORDER,varname_list, value_list))
		err("Variable FDORDER could not be retrieved from the json input file!");
    if (get_int_from_objectlist("FDORDER_TIME",number_readobjects,&FDORDER_TIME,varname_list, value_list))
        err("Variable FDORDER_TIME could not be retrieved from the json input file!");
	if (get_int_from_objectlist("FDCOEFF",number_readobjects,&FDCOEFF,varname_list, value_list))
		err("Variable FDCOEFF could not be retrieved from the json input file!");
	if (get_int_from_objectlist("NX",number_readobjects,&NX,varname_list, value_list))
		err("Variable NX could not be retrieved from the json input file!");
	if (get_int_from_objectlist("NY",number_readobjects,&NY,varname_list, value_list))
		err("Variable NY could not be retrieved from the json input file!");
	if (get_int_from_objectlist("NZ",number_readobjects,&NZ,varname_list, value_list))
		err("Variable NZ could not be retrieved from the json input file!");
	/*note that "y" is used for the vertical coordinate */
	if (get_float_from_objectlist("DX",number_readobjects,&DX,varname_list, value_list))
		err("Variable DX could not be retrieved from the json input file!");
	if (get_float_from_objectlist("DY",number_readobjects,&DY,varname_list, value_list))
		err("Variable DY could not be retrieved from the json input file!");
	if (get_float_from_objectlist("DZ",number_readobjects,&DZ,varname_list, value_list))
		err("Variable DZ could not be retrieved from the json input file!");
	/*note that "y" is used for the vertical coordinate */
	if (get_float_from_objectlist("TIME",number_readobjects,&TIME,varname_list, value_list))
		err("Variable TIME could not be retrieved from the json input file!");
	if (get_float_from_objectlist("DT",number_readobjects,&DT,varname_list, value_list))
		err("Variable DT could not be retrieved from the json input file!");

	/*=================================
		  section snapshot parameters
	  =================================*/
	if (get_int_from_objectlist("SNAP",number_readobjects,&SNAP,varname_list, value_list))
		err("Variable SNAP could not be retrieved from the json input file!");
	else{
		if (SNAP>0) {
			if (get_int_from_objectlist("SNAP_PLANE",number_readobjects,&SNAP_PLANE,varname_list, value_list))
				err("Variable SNAP_PLANE could not be retrieved from the json input file!");
			if (get_int_from_objectlist("SNAP_FORMAT",number_readobjects,&SNAP_FORMAT,varname_list, value_list))
				err("Variable SNAP_FORMAT could not be retrieved from the json input file!");
			if (get_float_from_objectlist("TSNAP1",number_readobjects,&TSNAP1,varname_list, value_list))
				err("Variable TSNAP1 could not be retrieved from the json input file!");
			if (get_float_from_objectlist("TSNAP2",number_readobjects,&TSNAP2,varname_list, value_list))
				err("Variable TSNAP2 could not be retrieved from the json input file!");
			if (get_float_from_objectlist("TSNAPINC",number_readobjects,&TSNAPINC,varname_list, value_list))
				err("Variable TSNAPINC could not be retrieved from the json input file!");
			if (get_string_from_objectlist("SNAP_FILE",number_readobjects,SNAP_FILE,varname_list, value_list))
				err("Variable SNAP_FILE could not be retrieved from the json input file!");
		}
	}
	/* increments are read in any case, because they will be also used as increment for model output */
	if (get_int_from_objectlist("IDX",number_readobjects,&IDX,varname_list, value_list))
		err("Variable IDX could not be retrieved from the json input file!");
	if (get_int_from_objectlist("IDY",number_readobjects,&IDY,varname_list, value_list))
		err("Variable IDY could not be retrieved from the json input file!");
	if (get_int_from_objectlist("IDZ",number_readobjects,&IDZ,varname_list, value_list))
		err("Variable IDZ could not be retrieved from the json input file!");
	/*note that "y" is used for the vertical coordinate */

	/*=================================
		  section seismogramm parameters
	  =================================*/

	if (get_int_from_objectlist("SEISMO",number_readobjects,&SEISMO,varname_list, value_list))
		err("Variable SEISMO could not be retrieved from the json input file!");
	else {
		if (SEISMO>0) {
			if (get_int_from_objectlist("READREC",number_readobjects,&READREC,varname_list, value_list))
				err("Variable READREC could not be retrieved from the json input file!");
			else {
				if (READREC==0) {
					if (get_int_from_objectlist("NGEOPH",number_readobjects,&NGEOPH,varname_list, value_list))
						err("Variable NGEOPH could not be retrieved from the json input file!");
					if (get_float_from_objectlist("XREC1",number_readobjects,&XREC1,varname_list, value_list))
						err("Variable XREC1 could not be retrieved from the json input file!");
					if (get_float_from_objectlist("XREC2",number_readobjects,&XREC2,varname_list, value_list))
						err("Variable XREC2T could not be retrieved from the json input file!");
					if (get_float_from_objectlist("YREC1",number_readobjects,&YREC1,varname_list, value_list))
						err("Variable YREC1 could not be retrieved from the json input file!");
					if (get_float_from_objectlist("YREC2",number_readobjects,&YREC2,varname_list, value_list))
						err("Variable YREC2 could not be retrieved from the json input file!");
					if (get_float_from_objectlist("ZREC1",number_readobjects,&ZREC1,varname_list, value_list))
						err("Variable ZREC1 could not be retrieved from the json input file!");
					if (get_float_from_objectlist("ZREC2",number_readobjects,&ZREC2,varname_list, value_list))
						err("Variable ZREC2 could not be retrieved from the json input file!");
					/*note that "y" is used for the vertical coordinate */
				}
				else
					if (get_string_from_objectlist("REC_FILE",number_readobjects,REC_FILE,varname_list, value_list))
						err("Variable REC_FILE could not be retrieved from the json input file!");
			}

			if (get_string_from_objectlist("SEIS_FILE",number_readobjects,SEIS_FILE,varname_list, value_list))
				err("Variable SEIS_FILE could not be retrieved from the json input file!");
			if (get_float_from_objectlist("REFRECX",number_readobjects,&REFREC[1],varname_list, value_list))
				err("Variable REFRECX could not be retrieved from the json input file!");
			if (get_float_from_objectlist("REFRECY",number_readobjects,&REFREC[2],varname_list, value_list))
				err("Variable REFRECY could not be retrieved from the json input file!");
			if (get_float_from_objectlist("REFRECZ",number_readobjects,&REFREC[3],varname_list, value_list))
				err("Variable REFRECZ could not be retrieved from the json input file!");
			/*note that "y" is used for the vertical coordinate */

			if (get_int_from_objectlist("REC_ARRAY",number_readobjects,&REC_ARRAY,varname_list, value_list))
				err("Variable REC_ARRAY could not be retrieved from the json input file!");
			else {
				if (REC_ARRAY>0) {
					if (get_float_from_objectlist("REC_ARRAY_DEPTH",number_readobjects,&REC_ARRAY_DEPTH,varname_list, value_list))
						err("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
					if (get_float_from_objectlist("REC_ARRAY_DIST",number_readobjects,&REC_ARRAY_DIST,varname_list, value_list))
						err("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
					if (get_int_from_objectlist("DRX",number_readobjects,&DRX,varname_list, value_list))
						err("Variable DRX could not be retrieved from the json input file!");
					if (get_int_from_objectlist("DRZ",number_readobjects,&DRZ,varname_list, value_list))
						err("Variable DRZ could not be retrieved from the json input file!");
					/*note that "y" is used for the vertical coordinate */
				}
			}
			if (get_int_from_objectlist("NDT",number_readobjects,&NDT,varname_list, value_list))
				err("Variable NDT could not be retrieved from the json input file!");
			if (get_int_from_objectlist("NDTSHIFT",number_readobjects,&NDTSHIFT,varname_list, value_list))
				err("Variable NDTSHIFT could not be retrieved from the json input file!");
			if (get_int_from_objectlist("SEIS_FORMAT",number_readobjects,&SEIS_FORMAT[0],varname_list, value_list))
				err("Variable SEIS_FORMAT could not be retrieved from the json input file!");
			else {
				if (get_int_from_objectlist("SEIS_FORMAT1",number_readobjects,&SEIS_FORMAT[1],varname_list, value_list)) {
					strcpy(varname_tmp1,"SEIS_FORMAT1");
					strcpy(value_tmp1,"0");
					add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
				}
				if (get_int_from_objectlist("SEIS_FORMAT2",number_readobjects,&SEIS_FORMAT[2],varname_list, value_list)) {
					strcpy(varname_tmp1,"SEIS_FORMAT2");
					strcpy(value_tmp1,"0");
					add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
				}
				if (get_int_from_objectlist("SEIS_FORMAT3",number_readobjects,&SEIS_FORMAT[3],varname_list, value_list)) {
					strcpy(varname_tmp1,"SEIS_FORMAT3");
					strcpy(value_tmp1,"0");
					add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
				}
				if (get_int_from_objectlist("SEIS_FORMAT4",number_readobjects,&SEIS_FORMAT[4],varname_list, value_list)) {
					strcpy(varname_tmp1,"SEIS_FORMAT4");
					strcpy(value_tmp1,"0");
					add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
				}
				if (get_int_from_objectlist("SEIS_FORMAT5",number_readobjects,&SEIS_FORMAT[6],varname_list, value_list)) {
					strcpy(varname_tmp1,"SEIS_FORMAT5");
					strcpy(value_tmp1,"0");
					add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);

				}
				if (get_int_from_objectlist("SEIS_FORMAT6",number_readobjects,&SEIS_FORMAT[6],varname_list, value_list)) {
					strcpy(varname_tmp1,"SEIS_FORMAT6");
					strcpy(value_tmp1,"0");
					add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
				}

				if (SEIS_FORMAT[0]==4) {
					SEIS_FORMAT[0]=0;
					/* fprintf(stderr," SEIS_FORMAT[1] set to %d.\n",SEIS_FORMAT[1]); */
				}
				if (SEIS_FORMAT[0]==5) {
					SEIS_FORMAT[0]=0;
					SEIS_FORMAT[1]=1;
					/* fprintf(stderr," SEIS_FORMAT[1] set to %d.\n",SEIS_FORMAT[1]); */
					SEIS_FORMAT[2]=1;
					/* fprintf(stderr," SEIS_FORMAT[2] set to %d.\n",SEIS_FORMAT[2]); */
					SEIS_FORMAT[3]=1;
					/* fprintf(stderr," SEIS_FORMAT[3] set to %d.\n",SEIS_FORMAT[3]); */
					SEIS_FORMAT[4]=0;
					SEIS_FORMAT[5]=0;
				}
			}

		}
	}

	/*=================================
		  section boundary parameters
	  =================================*/

	if (get_int_from_objectlist("FREE_SURF",number_readobjects,&FREE_SURF,varname_list, value_list))
		err("Variable FREE_SURF could not be retrieved from the json input file!");
	if (get_int_from_objectlist("BOUNDARY",number_readobjects,&BOUNDARY,varname_list, value_list))
		err("Variable BOUNDARY could not be retrieved from the json input file!");
	if (get_int_from_objectlist("ABS_TYPE",number_readobjects,&ABS_TYPE,varname_list, value_list))
		err("Variable ABS_TYPE could not be retrieved from the json input file!");

	if (ABS_TYPE==1) {
		if (get_float_from_objectlist("NPOWER",number_readobjects,&NPOWER,varname_list, value_list)) {
			strcpy(varname_tmp1,"NPOWER");
			strcpy(value_tmp1,"2");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);

		}
		if (get_float_from_objectlist("K_MAX_CPML",number_readobjects,&K_MAX_CPML,varname_list, value_list)) {
			strcpy(varname_tmp1,"K_MAX_CPML");
			strcpy(value_tmp1,"10");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}

		if (get_float_from_objectlist("FPML",number_readobjects,&FPML,varname_list, value_list))
			err("Variable FPML could not be retrieved from the json input file!");
		if (get_float_from_objectlist("VPPML",number_readobjects,&VPPML,varname_list, value_list))
			err("Variable VPPML could not be retrieved from the json input file!");
	}

	if (get_int_from_objectlist("FW",number_readobjects,&FW,varname_list, value_list))
		err("Variable FW could not be retrieved from the json input file!");
	if (get_float_from_objectlist("DAMPING",number_readobjects,&DAMPING,varname_list, value_list))
		err("Variable DAMPING could not be retrieved from the json input file!");

	/*=================================
		  section source parameters
	  =================================*/

	if (get_int_from_objectlist("SOURCE_SHAPE",number_readobjects,&SOURCE_SHAPE,varname_list, value_list))
		err("Variable SOURCE_SHAPE could not be retrieved from the json input file!");
	else {
		if (get_float_from_objectlist("FC",number_readobjects,&FC,varname_list, value_list)) {
			strcpy(varname_tmp1,"FC");
			strcpy(value_tmp1,"0.0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_float_from_objectlist("AMP",number_readobjects,&AMP,varname_list, value_list)) {
			strcpy(varname_tmp1,"AMP");
			strcpy(value_tmp1,"1.0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
	}
	if (SOURCE_SHAPE==3) {
		if (get_string_from_objectlist("SIGNAL_FILE",number_readobjects,SIGNAL_FILE,varname_list, value_list))
			err("Variable SIGNAL_FILE could not be retrieved from the json input file!");
		else {
			if (get_int_from_objectlist("SIGNAL_FORMAT0",number_readobjects,&SIGNAL_FORMAT[0],varname_list, value_list)) {
				strcpy(varname_tmp1,"SIGNAL_FORMAT0");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SIGNAL_FORMAT1",number_readobjects,&SIGNAL_FORMAT[1],varname_list, value_list)) {
				strcpy(varname_tmp1,"SIGNAL_FORMAT1");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SIGNAL_FORMAT2",number_readobjects,&SIGNAL_FORMAT[2],varname_list, value_list)) {
				strcpy(varname_tmp1,"SIGNAL_FORMAT2");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SIGNAL_FORMAT3",number_readobjects,&SIGNAL_FORMAT[3],varname_list, value_list)) {
				strcpy(varname_tmp1,"SIGNAL_FORMAT3");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SIGNAL_FORMAT4",number_readobjects,&SIGNAL_FORMAT[4],varname_list, value_list)) {
				strcpy(varname_tmp1,"SIGNAL_FORMAT0");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SIGNAL_FORMAT5",number_readobjects,&SIGNAL_FORMAT[5],varname_list, value_list)) {
				strcpy(varname_tmp1,"SIGNAL_FORMAT0");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
		}
	}
	if (get_int_from_objectlist("SRCREC",number_readobjects,&SRCREC,varname_list, value_list))
		err("Variable SRCREC could not be retrieved from the json input file!");
	else {
		if (get_int_from_objectlist("SRC_MF",number_readobjects,&SRC_MF,varname_list, value_list)) {
			strcpy(varname_tmp1,"SRC_MF");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_float_from_objectlist("REFSRC0",number_readobjects,&REFSRC[0],varname_list, value_list)) {
			strcpy(varname_tmp1,"REFSRC0");
			strcpy(value_tmp1,"0.0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_float_from_objectlist("REFSRC1",number_readobjects,&REFSRC[1],varname_list, value_list)) {
			strcpy(varname_tmp1,"REFSRC1");
			strcpy(value_tmp1,"0.0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_float_from_objectlist("REFSRC2",number_readobjects,&REFSRC[2],varname_list, value_list)) {
			strcpy(varname_tmp1,"REFSRC2");
			strcpy(value_tmp1,"0.0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}

	}
	if (SRCREC==1) {
		if (get_string_from_objectlist("SOURCE_FILE",number_readobjects,SOURCE_FILE,varname_list, value_list))
			err("Variable SOURCE_FILE could not be retrieved from the json input file!");
		else {
			if (get_string_from_objectlist("SRCOUT_FILE",number_readobjects,SRCOUT_FILE,varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_FILE");
				strcpy(value_tmp1,"");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SRCOUT_PAR0",number_readobjects,&SRCOUT_PAR[0],varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_PAR0");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SRCOUT_PAR1",number_readobjects,&SRCOUT_PAR[1],varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_PAR1");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SRCOUT_PAR2",number_readobjects,&SRCOUT_PAR[2],varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_PAR2");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SRCOUT_PAR3",number_readobjects,&SRCOUT_PAR[3],varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_PAR3");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SRCOUT_PAR4",number_readobjects,&SRCOUT_PAR[4],varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_PAR4");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("SRCOUT_PAR5",number_readobjects,&SRCOUT_PAR[5],varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCOUT_PAR5");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
		}
		if (get_int_from_objectlist("RUN_MULTIPLE_SHOTS",number_readobjects,&RUN_MULTIPLE_SHOTS,varname_list, value_list))
			err("Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");
		else {
			if (get_float_from_objectlist("SRCTSHIFT",number_readobjects,&SRCTSHIFT,varname_list, value_list)) {
				strcpy(varname_tmp1,"SRCTSHIFT");
				strcpy(value_tmp1,"0.0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("FSRC",number_readobjects,&FSRC,varname_list, value_list)) {
				strcpy(varname_tmp1,"FSRC");
				strcpy(value_tmp1,"0");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("JSRC",number_readobjects,&JSRC,varname_list, value_list)) {
				strcpy(varname_tmp1,"JSRC");
				strcpy(value_tmp1,"1");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
			if (get_int_from_objectlist("LSRC",number_readobjects,&LSRC,varname_list, value_list)) {
				strcpy(varname_tmp1,"LSRC");
				strcpy(value_tmp1,"2147483647");
				add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
			}
		}
	}

	if (SRCREC==2) {
		if (get_float_from_objectlist("PLANE_WAVE_DEPTH",number_readobjects,&PLANE_WAVE_DEPTH,varname_list, value_list))
			err("Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");
		else {
			if (PLANE_WAVE_DEPTH>0) {
				if (get_float_from_objectlist("PLANE_WAVE_ANGLE",number_readobjects,&PLANE_WAVE_ANGLE,varname_list, value_list))
					err("Variable PLANE_WAVE_ANGLE could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TS",number_readobjects,&TS,varname_list, value_list))
					err("Variable TS could not be retrieved from the json input file!");
				else {
					if (get_float_from_objectlist("SRC_DT",number_readobjects,&SRC_DT,varname_list, value_list)) {
						strcpy(varname_tmp1,"SRC_DT");
						get_string_from_objectlist("DT",number_readobjects,tempstring,varname_list, value_list);
						strcpy(value_tmp1,tempstring);
						add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
						SRC_DT=DT;
					}
				}
			}
		}
	}

	if (get_int_from_objectlist("SOURCE_TYPE",number_readobjects,&SOURCE_TYPE,varname_list, value_list))
		err("Variable SOURCE_TYPE could not be retrieved from the json input file!");
	if (SOURCE_TYPE==5) {
		if (get_float_from_objectlist("SOURCE_ALPHA",number_readobjects,&SOURCE_ALPHA,varname_list, value_list))
			err("Variable SOURCE_ALPHA could not be retrieved from the json input file!");
		if (get_float_from_objectlist("SOURCE_BETA",number_readobjects,&SOURCE_BETA,varname_list, value_list))
			err("Variable SOURCE_BETA could not be retrieved from the json input file!");
	}
	if (SOURCE_TYPE==6) {
		if (get_float_from_objectlist("AMON",number_readobjects,&AMON,varname_list, value_list))
			err("Variable AMON could not be retrieved from the json input file!");
		if (get_float_from_objectlist("STR",number_readobjects,&STR,varname_list, value_list))
			err("Variable STR could not be retrieved from the json input file!");
		if (get_float_from_objectlist("DIP",number_readobjects,&DIP,varname_list, value_list))
			err("Variable DIP could not be retrieved from the json input file!");
		if (get_float_from_objectlist("RAKE",number_readobjects,&RAKE,varname_list, value_list))
			err("Variable RAKE could not be retrieved from the json input file!");
	}


	/*=================================
	section model and log file parameters
	  =================================*/

	if (get_string_from_objectlist("MFILE",number_readobjects,MFILE,varname_list, value_list))
		err("Variable MFILE could not be retrieved from the json input file!");
	if (get_int_from_objectlist("WRITE_MODELFILES",number_readobjects,&WRITE_MODELFILES,varname_list, value_list)) {
		strcpy(varname_tmp1,"WRITE_MODELFILES");
		strcpy(value_tmp1,"2");
		add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
	}
	if (get_string_from_objectlist("LOG_FILE",number_readobjects,LOG_FILE,varname_list, value_list))
		err("Variable LOG_FILE could not be retrieved from the json input file!");
	if (get_string_from_objectlist("CHECKPT_FILE",number_readobjects,CHECKPTFILE,varname_list, value_list))
		err("Variable CHECKPT_FILE could not be retrieved from the json input file!");
	if (get_int_from_objectlist("LOG",number_readobjects,&LOG,varname_list, value_list))
		err("Variable LOG could not be retrieved from the json input file!");
	if (get_int_from_objectlist("CHECKPTREAD",number_readobjects,&CHECKPTREAD,varname_list, value_list))
		err("Variable CHECKPTREAD could not be retrieved from the json input file!");
	if (get_int_from_objectlist("CHECKPTWRITE",number_readobjects,&CHECKPTWRITE,varname_list, value_list))
		err("Variable CHECKPTWRITE could not be retrieved from the json input file!");
	if (get_int_from_objectlist("OUT_TIMESTEP_INFO",number_readobjects,&OUTNTIMESTEPINFO,varname_list, value_list)) {
		strcpy(varname_tmp1,"OUT_TIMESTEP_INFO");
		strcpy(value_tmp1,"1");
		add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
	}
	if (get_int_from_objectlist("OUT_SOURCE_WAVELET",number_readobjects,&OUTSOURCEWAVELET,varname_list, value_list)) {
		strcpy(varname_tmp1,"OUT_SOURCE_WAVELET");
		strcpy(value_tmp1,"0");
		add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
	}

	if (get_int_from_objectlist("L",number_readobjects,&L,varname_list, value_list))
		err("Variable L could not be retrieved from the json input file!");
	else {
		FL=vector(1,L);
		switch(L) {
		case 0:
			break;
		case 1:
			if (get_float_from_objectlist("FL1",number_readobjects,&FL[1],varname_list, value_list))
				err("Variable FL1 could not be retrieved from the json input file!");
			break;
		default:
			err("More than one relaxation Parameter (L>1) is not implemented yet!");
			break;
		}
		(get_float_from_objectlist("FREF",number_readobjects,&FREF,varname_list, value_list));
		//	err("Variable FREF could not be retrieved from the json input file!");

		/* in the viscoelastic case : reference frequency where no velocity dispersion occurs */
	}
	if (L) {
		if (get_float_from_objectlist("TAU",number_readobjects,&TAU,varname_list, value_list))
			err("Variable TAU could not be retrieved from the json input file!");
	}

	if (get_int_from_objectlist("ASCIIEBCDIC",number_readobjects,&ASCIIEBCDIC,varname_list, value_list)) {
		strcpy(varname_tmp1,"ASCIIEBCDIC");
		strcpy(value_tmp1,"0");
		add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
	}
	if (get_int_from_objectlist("LITTLEBIG",number_readobjects,&LITTLEBIG,varname_list, value_list)) {
		strcpy(varname_tmp1,"LITTLEBIG");
		strcpy(value_tmp1,"0");
		add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
	}
	if (get_int_from_objectlist("IEEEIBM",number_readobjects,&IEEEIBM,varname_list, value_list)) {
		strcpy(varname_tmp1,"IEEEIBM");
		strcpy(value_tmp1,"0");
		add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
	}

	if (get_int_from_objectlist("READMOD",number_readobjects,&READMOD,varname_list, value_list))
		err("Variable READMOD could not be retrieved from the json input file!");
	else {
		if (get_int_from_objectlist("MOD_FORMAT0",number_readobjects,&MOD_FORMAT[0],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT0");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_int_from_objectlist("MOD_FORMAT1",number_readobjects,&MOD_FORMAT[1],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT1");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_int_from_objectlist("MOD_FORMAT2",number_readobjects,&MOD_FORMAT[2],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT2");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_int_from_objectlist("MOD_FORMAT3",number_readobjects,&MOD_FORMAT[3],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT3");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_int_from_objectlist("MOD_FORMAT4",number_readobjects,&MOD_FORMAT[4],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT4");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_int_from_objectlist("MOD_FORMAT5",number_readobjects,&MOD_FORMAT[5],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT5");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
		if (get_int_from_objectlist("MOD_FORMAT6",number_readobjects,&MOD_FORMAT[6],varname_list, value_list)) {
			strcpy(varname_tmp1,"MOD_FORMAT6");
			strcpy(value_tmp1,"0");
			add_object_tolist(varname_tmp1,value_tmp1,&number_defaultobjects, varnamedefault_list, valuedefault_list);
		}
	}
	fprintf(fp, "\n===========================================================");
	fprintf(fp, "\n= List of Parameters NOT read by the built in Json Parser =");
	fprintf(fp, "\n=      Values might not be in used, set to default        =");

	print_objectlist_screen(fp, number_defaultobjects, varnamedefault_list, valuedefault_list);

}

