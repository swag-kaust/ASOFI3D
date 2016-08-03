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
/******************************************************************************/
/* RWSEGY: reading and writing SEG-Y, SU, BIN, TXT and UKOOA P190             */
/*                                                                            */
/*                                                                            */
/* 	CAUTION: tested for LITTLEBIG=ASCIIEBCDIC=IEEEIBM=0 only              */
/*	         (no other suitable computer for further test available)      */
/*                                                                            */
/*               Some routines are not optimized for speed!                   */
/******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#ifndef STRING_SIZE
#define STRING_SIZE 74
#endif

#ifndef PI
#define PI (3.141592653589793)
#endif


#define mymin(x,y) ((x<y) ? x:y)
#define mymax(x,y) ((y<x) ? x:y)
inline int myisign(int x){ if (x>0) return 1; else if (x<0) return -1; else return 0;}
inline float myfsign(float x){ if (x>0.0) return 1; else if (x<0.0) return -1; else return 0;}


	float * m2ft(float * iovec, int len);
	float * ft2m(float * iovec, int len);
	short * swap2(short * ip);
	int * swap4(int * ip);
	float * swap4fv(float * fp, int len);
	int * ieee2ibmv(int * ieee, int n);
	int * ibm2ieeev(int * ibm, int n);
	unsigned char a2e(unsigned char ain);
	unsigned char e2a(unsigned char ein);
	char * a2ev(char * inout, int len);
	char * e2av(char * inout, int len);
	float * float2native(float * indata, int len, int lbendian, int ieeeibm, int meterfeet);
	float * native2float(float * outdata, int len, int lbendian, int ieeeibm, int meterfeet);
	/* const int a2etab[128]; */	
	
	extern FILE * FP;
	FILE * curerr=NULL;   /* stream for error messages etc. */
	
	extern int LITTLEBIG, ASCIIEBCDIC, IEEEIBM, MYID;
	/* LITTLEBIG	: 0 for little endian machines, 1 for big endian machines */
	/* ASCIIEBCDIC	: 0 native format for reading and writing characters shall be ascii or equivanlent (not ebcdic) */
	/* IEEEIBM	: native float format: 0=IEEE, 1=IBM */	
	

/* reading and writing UKOOA files --- caution: metric grid easting, northing and elevation are expected, only S and R entries, but no H-entries are evaluated */

/* UNDER CONSTRUCTION !!! */

	
/* reading and writing SU/SEG-Y trace headers */	
	
int DDN_rtraceh(FILE * instream, int lbendian, int ieeeibm, int meterfeet, int susegy, int * parvec[27]){
	/* structure from segy.h is not used in order to access variable data portion independently */
	/* however, using segy.h probably makes the code faster */

	/* internal variables */
	short dummyshort, scalel, scalco, trid;
	int * pint, dummyint, doswap=0, m, n=0;
	unsigned short /*uh=0,*/ udt, uns;

	/* to be returned */
	int ns, dt, traceid;
	int tracl, traceno, shotno, globrecno, comp, cdp, gcoo[3], scoo[3], scoo2, offset, gdel, sdel, swdep, gwdep;
	int xcoo, ycoo, inlineno, crosslineno, sno; /* segy only */
	float  d1, f1, d2, f2; /* su only */ 
	int ntraces; /* su only */
	
	if (lbendian!=LITTLEBIG) doswap=1;	
	/* READING TRACE HEADER - not optimized for speed! */
	
	/*   1-  4 int tracl;    */ n+=fread(&tracl,1,4,instream); if (doswap) swap4(&tracl); /* trace sequence number within line --- FOR SEISMOGRAMS: RECEIVER NO. FOR THIS SHOT */
	/*   5-  8 int tracr;    */ n+=fread(&traceno,1,4,instream); if (doswap) swap4(&traceno); /* trace sequence number within reel (jth trace in file) */
	/*   9- 12 int fldr;     */ n+=fread(&shotno,1,4,instream); if (doswap) swap4(&shotno); /* field record number --- FOR SEISMOGRAMS: CURENT SHOT NO. */
	/*  13- 16 int tracf;    */ n+=fread(&globrecno,1,4,instream); if (doswap) swap4(&globrecno); /* trace number within field record --- RECEIVER NO. */
	/*  17- 20 int ep;       */ n+=fread(&comp,1,4,instream); if (doswap) swap4(&comp); /* energy source point no. --- FOR SEISMOGRAMS: component, WARNING: ORIGINAL MEANING IS SOURCE POINT! */
	/*  21- 24 int cdp;      */ n+=fread(&cdp,1,4,instream); if (doswap) swap4(&cdp); /* CDP ensemble number --- FOR SEISMOGRAMS: RECEIVER NO. FOR THIS SHOT */
	/*  25- 28 int cdpt;     */ n+=fread(&dummyint,1,4,instream);  /* trace number within CDP ensemble */
	/*  29- 30 short trid;   */ n+=fread(&trid,1,2,instream); if (doswap) swap2(&trid); /* trace identification code:
			1 = seismic data (SEISMOGRAM)
			2 = dead (RECEIVER OUTSIDE MODEL)
			3 = dummy (SOURCE and possibly receiver OUTSIDE MODEL)
			4 = time break (SOURCE TIMES)
			6 = sweep (SOURCE SIGNATURE)
		FOR SU DATA: CWP id flags:
			43 = Seismic Data, Vertical Component 
			44 = Seismic Data, Horizontal Component 1 
			45 = Seismic Data, Horizontal Component 2 
			46 = Seismic Data, Radial Component
			47 = Seismic Data, Transverse Component  
		FOR SEG-Y DATA:
			11 = Seismic pressure sensor  
			12 = Multicomponent seismic sensor - Vertical component
			13 = Multicomponent seismic sensor - Cross-line component
			14 = Multicomponent seismic sensor - In-line component
			15 = Rotated multicomponent seismic sensor - Vertical component
			16 = Rotated multicomponent seismic sensor - Transverse component
			17 = Rotated multicomponent seismic sensor - Radial component */
	/*  31- 32 short nvs;    */ n+=fread(&dummyshort,1,2,instream); /* number of vertically summed traces (see vscode in bhed structure) */
	/*  33- 34 short nhs;    */ n+=fread(&dummyshort,1,2,instream); /* number of horizontally summed traces (see vscode in bhed structure) */
	/*  35- 36 short duse;   */ n+=fread(&dummyshort,1,2,instream); /* data use: 1 = production, 2 = test */
	/*  37- 40 int offset;   */ n+=fread(&offset,1,4,instream); if (doswap) swap4(&offset); /* distance from source point to receiver group 
	                                 (CAUTION: should be negative if opposite to direction in which the line was shot) */
	/*  41- 44 int gelev;    */ n+=fread(&gcoo[2],1,4,instream); if (doswap) swap4(&gcoo[2]); /* receiver group elevation from sea level (above sea level is positive) --- CAUTION: SEISMOGRAMS: Y rec. coordinate */
	/*  45- 48 int selev;    */ n+=fread(&scoo[2],1,4,instream); if (doswap) swap4(&scoo[2]); /* source elevation from sea level (above sea level is positive) --- CAUTION: SEISMOGRAMS: Y source coordinate  */
	/*  49- 52 int sdepth;   */ n+=fread(&scoo2,1,4,instream); if (doswap) swap4(&scoo2); /* source depth (positive) --- SEISMOGRAMS: Y source coordinate  */
	/*  53- 56 int gdel;     */ n+=fread(&gdel,1,4,instream); if (doswap) swap4(&gdel); /* datum elevation at receiver group --- WARNING: SEISMOGRAMS: DIRECTIVITY NON-STANDARD */
	/*  57- 60 int sdel;     */ n+=fread(&sdel,1,4,instream); if (doswap) swap4(&sdel); /* datum elevation at source */
	/*  61- 64 int swdep;    */ n+=fread(&swdep,1,4,instream); if (doswap) swap4(&sdel); /* water depth at source */
	/*  65- 68 int gwdep;    */ n+=fread(&gwdep,1,4,instream); if (doswap) swap4(&sdel); /* water depth at receiver group 
									--- WARNING: SEISMOGRAMS: DIRECTIVITY NON-STANDARD */
	/*  69- 70 short scalel; */ n+=fread(&scalel,1,2,instream); if (doswap) swap2(&scalel); /* scale factor for previous 7 entries with value plus or minus 10 to the power 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) */
	/*  71- 72 short scalco; */ n+=fread(&scalco,1,2,instream); if (doswap) swap2(&scalco); /* scale factor for next 4 entries with value plus or minus 10 to the power 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) */
	/*  73- 76 int  sx;      */ n+=fread(&scoo[0],1,4,instream); if (doswap) swap4(&scoo[0]); /* X source coordinate */
	/*  77- 80 int  sy;      */ n+=fread(&scoo[1],1,4,instream); if (doswap) swap4(&scoo[1]); /* Y source coordinate */
	/*  81- 84 int  gx;      */ n+=fread(&gcoo[0],1,4,instream); if (doswap) swap4(&gcoo[0]); /* X group coordinate  */
	/*  85- 88 int  gy;      */ n+=fread(&gcoo[1],1,4,instream); if (doswap) swap4(&gcoo[1]); /* Y group coordinate  */
	/*  89- 90 short counit; */ n+=fread(&dummyshort,1,2,instream); if (doswap) swap2(&dummyshort); /* coordinate units code: for previous four entries */
			if ((dummyshort!=0)&&(dummyshort!=1)) fprintf(curerr,"Warning [rtraceh]: co-ordinate unit flag is '%d'\n",dummyshort);
			/* 1 = length (meters or feet),	2 = seconds of arc (X = longitude; Y = latitude,
			   a positive value designates the number of seconds east of Greenwich	or north of the equator) */
	/*  91- 92 short wevel;  */ n+=fread(&dummyshort,1,2,instream); /* weathering velocity */
	/*  93- 94 short swevel; */ n+=fread(&dummyshort,1,2,instream); /* subweathering velocity */
	/*  95- 96 short sut;    */ n+=fread(&dummyshort,1,2,instream); /* uphole time at source */
	/*  97- 98 short gut;    */ n+=fread(&dummyshort,1,2,instream); /* uphole time at receiver group */
	/*  99-100 short sstat;	 */ n+=fread(&dummyshort,1,2,instream); /* source static correction */
	/* 101-102 short gstat;	 */ n+=fread(&dummyshort,1,2,instream); /* group static correction */
	/* 103-104 short tstat;	 */ n+=fread(&dummyshort,1,2,instream); /* total static applied */
		/* sh=(int)-floor(dtshift*1e6+0.5); if (doswap) swap2(&sh);*/ /* time break = start of the modeling = 0 */
	/* 105-106 short laga;   */ n+=fread(&dummyshort,1,2,instream);
			/* lag time A, time in ms between end of 240-
			   byte trace identification header and time
			   break, positive if time break occurs after
			   end of header, time break is defined as
			   the initiation pulse which maybe recorded
			   on an auxiliary trace or as otherwise
			   specified by the recording system */
		/* sh=(int)floor(srctime*1e6+0.5); if (doswap) swap2(&sh); */		   
	/* 107-108 short lagb;   */ n+=fread(&dummyshort,1,2,instream);
			/* lag time B, time in ms between the time break
			   and the initiation time of the energy source,
			   may be positive or negative */
		/* sh=(int)floor((dtshift-srctime)*1e6+0.5); if (doswap) swap2(&sh); */		   
	/* 109-110 short delrt;  */ n+=fread(&dummyshort,1,2,instream);
			/* delay recording time, time in ms between
			   initiation time of energy source and time
			   when recording of data samples begins
			   (for deep water work if recording does not
			   start at zero time) */
	/* 111-112 short muts;   */ n+=fread(&dummyshort,1,2,instream); /* mute time--start */
	/* 113-114 short mute;   */ n+=fread(&dummyshort,1,2,instream); /* mute time--end */
	/* 115-116 unsigned short ns; */ n+=fread(&uns,1,2,instream); swap2((short *) &uns); /* number of samples in this trace */
	/* 117-118 unsigned short dt; */ n+=fread(&udt,1,2,instream); swap2((short *) &udt); /* sample interval; in micro-seconds */
	/* 119-120 short gain;	 */ n+=fread(&dummyshort,1,2,instream); /* gain type of field instruments code:
				1 = fixed,2 = binary, 3 = floating point, 4 ---- N = optional use */
	/* 121-122 short igc;	 */ n+=fread(&dummyshort,1,2,instream); /* instrument gain constant */
	/* 123-124 short igi;	 */ n+=fread(&dummyshort,1,2,instream); /* instrument early or initial gain */
	/* 125-126 short corr;	 */ n+=fread(&dummyshort,1,2,instream); /* correlated: 1 = no, 2 = yes */
	/* 127-128 short sfs;	 */ n+=fread(&dummyshort,1,2,instream); /* sweep frequency at start */
	/* 129-130 short sfe;	 */ n+=fread(&dummyshort,1,2,instream); /* sweep frequency at end */
	/* 131-132 short slen;	 */ n+=fread(&dummyshort,1,2,instream); /* sweep length in ms */
	/* 133-134 short styp;	 */ n+=fread(&dummyshort,1,2,instream); /* sweep type code:
				1 = linear, 2 = parapolic,   3 = exponential, 4 = other */
	/* 135-136 short stas;   */ n+=fread(&dummyshort,1,2,instream); /* sweep trace taper length at start in ms */
	/* 137-138 short stae;   */ n+=fread(&dummyshort,1,2,instream); /* sweep trace taper length at end in ms */
	/* 139-140 short tatyp;  */ n+=fread(&dummyshort,1,2,instream); /* sweep taper type code : 1=linear, 2=cos^2, 3=other */
	/* 141-142 short afilf;  */ n+=fread(&dummyshort,1,2,instream); /* alias filter frequency if used */
	/* 143-144 short afils;  */ n+=fread(&dummyshort,1,2,instream); /* alias filter slope */
	/* 145-146 short nofilf; */ n+=fread(&dummyshort,1,2,instream); /* notch filter frequency if used */
	/* 147-148 short nofils; */ n+=fread(&dummyshort,1,2,instream); /* notch filter slope */
	/* 149-150 short lcf;	 */ n+=fread(&dummyshort,1,2,instream); /* low cut frequency if used */
	/* 151-152 short hcf;	 */ n+=fread(&dummyshort,1,2,instream); /* high cut frequncy if used */
	/* 153-154 short lcs;	 */ n+=fread(&dummyshort,1,2,instream); /* low cut slope */
	/* 155-156 short hcs;	 */ n+=fread(&dummyshort,1,2,instream); /* high cut slope */
	/* 157-158 short year;   */ n+=fread(&dummyshort,1,2,instream); /* year data recorded */
	/* 159-160 short day;	 */ n+=fread(&dummyshort,1,2,instream); /* day of year */
	/* 161-162 short hour;   */ n+=fread(&dummyshort,1,2,instream); /* hour of day (24 hour clock) */
	/* 163-164 short minute; */ n+=fread(&dummyshort,1,2,instream); /* minute of hour */
	/* 165-166 short sec;	 */ n+=fread(&dummyshort,1,2,instream); /* second of minute */
	/* 167-168 short timbas; */ n+=fread(&dummyshort,1,2,instream); /* time basis code: 1 = local, 2 = GMT, 3 = other */
	/* 169-170 short trwf;	 */ n+=fread(&dummyshort,1,2,instream); /* trace weighting factor, i.e. 1/2^N  volts for the least sigificant bit */
	/* 171-172 short grnors; */ n+=fread(&dummyshort,1,2,instream); /* geophone group number of roll switch  position one */
	/* 173-174 short grnofr; */ n+=fread(&dummyshort,1,2,instream); /* geophone group number of trace one within original field record */
	/* 175-176 short grnlof; */ n+=fread(&dummyshort,1,2,instream); /* geophone group number of last trace within original field record */
	/* 177-178 short gaps;	 */ n+=fread(&dummyshort,1,2,instream); /* gap size (total number of groups dropped) */
	/* 179-180 short otrav;	 */ n+=fread(&dummyshort,1,2,instream); /* overtravel taper code: 1 = down (or behind), 2 = up (or ahead) */
	if (susegy==1){
	 /* 181-184 int ???;      */ n+=fread(&xcoo,1,4,instream); if (doswap) swap4(&xcoo); /* 1st horizontal coordinate for models (in SOFI3D ususally X) */
	 							      /* X coordinate of ensemble (CDP) position of this trace (scalar in Trace
									 Header bytes 71-72 applies). The coordinate reference system should be
									 identified through an extended header Location Data stanza. */
	 /* 185-188 int ???;      */ n+=fread(&ycoo,1,4,instream); if (doswap) swap4(&ycoo); /* 2nd horizontal coordinate for models (in SOFI3D ususally Z) */   
	 							      /* Y coordinate of ensemble (CDP) position of this trace (scalar in Trace
									 Header bytes 71-72 applies). The coordinate reference system should be
									 identified through an extended header Location Data stanza. */
	 /* 189-192 int ???;      */ n+=fread(&inlineno,1,4,instream);   /* For 3-D poststack data, this field should be used for the in-line
	 								 number. If one in-line per SEG Y file is being recorded, this value 
									 should be the same for all traces in the file and the same value will 
									 be recorded in bytes 3205-3208 of the Binary File Header. */
	 /* 193-196 int ???;      */ n+=fread(&crosslineno,1,4,instream);   /* For 3-D poststack data, this field should be used for the cross-line
									 number.  This will typically be the same value as the ensemble (CDP)
									 number in Trace Header bytes 21-24, but this does not have to be the 
									 case. */
	 /* 197-200 int ???;      */ n+=fread(&sno,1,4,instream); if (doswap) swap4(&sno);
	 							      /* Shotpoint number - This is probably only applicable to 2-D poststack 
									 data. - Note that it is assumed that the shotpoint number refers to the
									 source location nearest to the ensemble (CDP) location for a particular
									 trace.  If this is not the case, there should be a comment in the Textual
									 File Header explaining what the  shotpoint number actually refers to. */
	 /* 201-202 short ???;    */ n+=fread(&dummyshort,1,2,instream); /* Scalar to be applied to the shotpoint number in Trace Header bytes 
									 197-200 to  give the real value. If positive, scalar is used as 
									 multiplier; if negative as a divisor; if zero the shotpoint number is not
									 scaled (i.e. it is an integer. A typical value will be -10, allowing
									 shotpoint numbers with one decimal digit to the right of the decimal
									 point). */
									 
	 /* 203-204 short ???;    */ n+=fread(&dummyshort,1,2,instream); /* Trace value measurement unit:  
									 -1 = Other (should be described in Data Sample Measurement Units Stanza) 
								          0 = Unknown    
									  1 = Pascal (Pa)    
									  2 = Volts (v)  
									  3 = Millivolts (mV) 
  									  4 = Amperes (A) 
 									  5 = Meters (m) 
									  6 = Meters per second (m/s) 
									  7 = Meters per second squared (m/s2) 
									  8 = Newton (N) 
									  9 = Watt (W) 	*/
									  
	 /* 205-208 int ???;    */ n+=fread(&dummyint,1,4,instream);  /* see bytes 209-210! */
	 /* 209-210 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Transduction Constant - The multiplicative constant used to convert the 
									 Data  Trace samples to the Transduction Units (specified in Trace Header
									 bytes 211- 212).  The constant is encoded as a four-byte, two's 
									 complement integer (bytes 205-208) which is the mantissa and a two-byte,
									 two's complement integer (bytes 209-210) which is the power of ten
									 exponent (i.e. Bytes 205-208 * 10**Bytes  209-210).  */
									 
	 /* 211-212 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Transduction Units - The unit of measurement of the Data Trace samples
									 after  they have been multiplied by the Transduction Constant specified 
									 in Trace  Header bytes 205-210.
									 0 = Unknown    
									 1 = Pascal (Pa)
									 2 = Volts (v)
									 3 = Millivolts (mV) 
									 4 = Amperes (A) 
									 5 = Meters (m) 
									 6 = Meters per second (m/s) 
									 7 = Meters per second squared (m/s2)
									 8 = Newton (N) 
									 9 = Watt (W)  */
	 /* 213-214 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Device/Trace Identifier ? The unit number or id number of the device
									 associated  with the Data Trace (i.e. 4368 for vibrator serial number 
									 4368 or 20316 for gun 16 on string 3 on vessel 2).  This field allows
									 traces to be associated across trace ensembles independently of the 
									 trace number (Trace Header bytes 25-28). */
	 /* 215-216 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Scalar to be applied to times specified in Trace Header bytes 95-114 to
									 give the  true time value in milliseconds.  Scalar = 1, +10, +100, +1000,
									 or +10,000.  If positive, scalar is used as a multiplier; if negative,
									 scalar is used as divisor.  A  value of zero is assumed to be a scalar
									 value of 1. */
									 
	 /* 217-218 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Source Type/Orientation ? Defines the type and the orientation of the
									 energy  source.  The terms vertical, cross-line and in-line refer to the
									 three axes of an  orthogonal coordinate system.  The absolute azimuthal
									 orientation of the  coordinate system axes can be defined in the Bin Grid
									 Definition Stanza.  
									 -1 to -n = Other (should be described in Source Type/Orientation stanza) 
									 0 = Unknown 
									 1 = Vibratory - Vertical orientation 	
									 2 = Vibratory - Cross-line orientation
									 3 = Vibratory - In-line orientation
									 4 = Impulsive - Vertical orientation 
									 5 = Impulsive - Cross-line orientation 
									 6 = Impulsive - In-line orientation 
									 7 = Distributed Impulsive - Vertical orientation 
									 8 = Distributed Impulsive - Cross-line orientation 
									 9 = Distributed Impulsive - In-line orientation */
 	 /* 219-224 ?!?!? ???;  */ n+=fread(&dummyshort,1,2,instream);
		    		   n+=fread(&dummyint,1,4,instream);  /* Source Energy Direction with respect to the source orientation  - The
									 positive  orientation direction is defined in Bytes 217-218 of the Trace
									 Header.  The energy direction is encoded in tenths of degrees 
									 (i.e. 347.8� is encoded as 3478). */ 							 
	 /* 225-228 int ???;    */ n+=fread(&dummyint,1,4,instream);  /* see bytes 229-230! */
	 /* 229-230 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Measurement - Describes the source effort used to generate the
									 trace.   The measurement can be simple, qualitative measurements such as
									 the total  weight of explosive used or the peak air gun pressure or the
									 number of vibrators  	times the sweep duration.  Although these simple
									 measurements are acceptable,  it is preferable to use true measurement
									 units of energy or work.  The constant is encoded as a four-byte, two's
									 complement integer (bytes 225-228) which is the mantissa and a two-byte,
									 two's complement integer (bytes 209-230) which is the power of ten
									 exponent (i.e. Bytes 225-228 * 10**Bytes 229- 230). */
									 
	 /* 231-232 short ???;  */ n+=fread(&dummyshort,1,2,instream); /* Source Measurement Unit - The unit used for the Source Measurement, 
	 								   Trace header bytes 225-230.  
									 -1 = Other (should be described in Source Measurement Unit stanza)   
									  0 = Unknown    
									  1 = Joule (J)
									  2 = Kilowatt (kW)  
									  3 = Pascal (Pa) 
									  4 = Bar (Bar) 
									  4 = Bar-meter (Bar-m) 
									  5 = Newton (N) 
									  6 = Kilograms (kg) */
  
	 /* 233-240 int[2] ???;    */ n+=fread(&dummyint,1,4,instream); 
	 			      n+=fread(&dummyint,1,4,instream); /* Unassigned ? For optional information. */
	}
	else { /* local SU assignments */
	   /* 181-184 float d1;     */ n+=fread(&d1,1,4,instream); float2native(&d1,1,lbendian,ieeeibm,0); /* sample spacing for non-seismic data */
	   /* 185-188 float f1;     */ n+=fread(&f1,1,4,instream); float2native(&f1,1,lbendian,ieeeibm,0); /* first sample location for non-seismic data */
	   /* 189-192 float d2;     */ n+=fread(&d2,1,4,instream); float2native(&d2,1,lbendian,ieeeibm,0); /* sample spacing between traces */
	   /* 193-296 float f2;     */ n+=fread(&f2,1,4,instream); float2native(&f2,1,lbendian,ieeeibm,0); /* first trace location */
	   /* 197-200 float ungpow; */ n+=fread(&dummyint,1,4,instream); /* negative of power used for dynamic range compression */
	   /* 201-204 float unscale;*/ n+=fread(&dummyint,1,4,instream); /* reciprocal of scaling factor to normalize range */
	   /* 205-208 int ntr;      */ n+=fread(&ntraces,1,4,instream); if (doswap) swap4(&ntraces); /* number of traces --- NUMBER OF TRACE IN FILE */
	   /* 209-210 short mark;   */ n+=fread(&dummyshort,1,2,instream); /* mark selected traces */
       	   /* 211-212 short shortpad; */ n+=fread(&dummyshort,1,2,instream); /* alignment padding */
	   /* 213-240 short unass[14] */ for (m=0;m<14;m++) n+=fread(&dummyshort,1,2,instream); /* unassigned */
	}
	if (n!=240) fprintf(curerr,"Warning [rtraceh]: segy-header consists of %d bytes instead of 240?!\n",n);

	dt=udt*1.0e-6;
	ns=(int)uns;
	traceid=(int)trid;
	if (meterfeet==1){
	    offset=offset*0.3048;
	    d1=d1*0.3048;
	    f1=f1*0.3048;
	    d2=d2*0.3048;
	    f2=f2*0.3048;	
	    if (susegy==1){ /* SEG-Y */
		if (scalel>0){
			/* NOT TO BE SCALED: offset */
			gcoo[2]=scalel*gcoo[2]*0.3048; /* gelev */
			scoo[2]=scalel*scoo[2]*0.3048; /* selev */
			scoo2=scalel*scoo2*0.3048;     /* sdepth */
			gdel=scalel*gdel*0.3048;
			sdel=scalel*sdel*0.3048;
		} else if (scalel<0){
			/* NOT TO BE SCALED: offset */
			gcoo[2]=gcoo[2]*0.3048/-scalel; /* gelev */
			scoo[2]=scoo[2]*0.3048/-scalel; /* selev */
			scoo2=scoo2*0.3048/-scalel;     /* sdepth */
			gdel=gdel*0.3048/-scalel;
			sdel=sdel*0.3048/-scalel;		
		} else {
			gcoo[2]=gcoo[2]*0.3048; /* gelev */
			scoo[2]=scoo[2]*0.3048; /* selev */
			scoo2=scoo2*0.3048;     /* sdepth */
			gdel=gdel*0.3048;
			sdel=sdel*0.3048;
		}		
		if (scalco>0){
		 	scoo[0]=scalco*scoo[0]*0.3048; /* sx */
			scoo[1]=scalco*scoo[1]*0.3048; /* sy */
		 	gcoo[0]=scalco*gcoo[0]*0.3048; /* gx */
			gcoo[1]=scalco*gcoo[1]*0.3048; /* gy */
			xcoo=scalco*xcoo*0.3048;
			ycoo=scalco*ycoo*0.3048;			
		} else if (scalco<0){
		 	scoo[0]=scoo[0]*0.3048/-scalco; /* sx */
			scoo[1]=scoo[1]*0.3048/-scalco; /* sy */
		 	gcoo[0]=gcoo[0]*0.3048/-scalco; /* gx */
			gcoo[1]=gcoo[1]*0.3048/-scalco; /* gy */
			xcoo=xcoo*0.3048/-scalco;
			ycoo=ycoo*0.3048/-scalco;		
		} else {
		 	scoo[0]=scoo[0]*0.3048; /* sx */
			scoo[1]=scoo[1]*0.3048; /* sy */
		 	gcoo[0]=gcoo[0]*0.3048; /* gx */
			gcoo[1]=gcoo[1]*0.3048; /* gy */
			xcoo=xcoo*0.3048;
			ycoo=ycoo*0.3048;		
		}
	    }	
	    else { /* SU */
			gcoo[2]=pow(10.0,scalel)*gcoo[2]*0.3048; /* gelev */
			scoo[2]=pow(10.0,scalel)*scoo[2]*0.3048; /* selev */
			scoo2=pow(10.0,scalel)*scoo2*0.3048;     /* sdepth */
			gdel=pow(10.0,scalel)*gdel*0.3048;
			sdel=pow(10.0,scalel)*sdel*0.3048;
		 	scoo[0]=pow(10.0,scalco)*scoo[0]*0.3048; /* sx */
			scoo[1]=pow(10.0,scalco)*scoo[1]*0.3048; /* sy */
		 	gcoo[0]=pow(10.0,scalco)*gcoo[0]*0.3048; /* gx */
			gcoo[1]=pow(10.0,scalco)*gcoo[1]*0.3048; /* gy */			
	    }	
	}
	else{
	    if (susegy==1) { /* SEG-Y */
		if (scalel>0){
			/* NOT TO BE SCALED: offset */
			gcoo[2]=scalel*gcoo[2]; /* gelev */
			scoo[2]=scalel*scoo[2]; /* selev */
			scoo2=scalel*scoo2;     /* sdepth */
			gdel=scalel*gdel;
			sdel=scalel*sdel;
		} else if (scalel<0){
			/* NOT TO BE SCALED: offset */
			gcoo[2]=gcoo[2]/-scalel; /* gelev */
			scoo[2]=scoo[2]/-scalel; /* selev */
			scoo2=scoo2/-scalel;     /* sdepth */
			gdel=gdel/-scalel;
			sdel=sdel/-scalel;		
		}
		if (scalco>0){
		 	scoo[0]=scalco*scoo[0]; /* sx */
			scoo[1]=scalco*scoo[1]; /* sy */
		 	gcoo[0]=scalco*gcoo[0]; /* gx */
			gcoo[1]=scalco*gcoo[1]; /* gy */
			xcoo=scalco*xcoo;
			ycoo=scalco*ycoo;			
		} else if (scalco<0){
		 	scoo[0]=scoo[0]/-scalco; /* sx */
			scoo[1]=scoo[1]/-scalco;; /* sy */
		 	gcoo[0]=gcoo[0]/-scalco;; /* gx */
			gcoo[1]=gcoo[1]/-scalco;; /* gy */
			xcoo=xcoo/-scalco;
			ycoo=ycoo/-scalco;		
		}
	    }		
	    else /* SU */ {
		if (scalel!=0){
			/* NOT TO BE SCALED: offset */
			gcoo[2]=pow(10.0,scalel)*gcoo[2]; /* gelev */
			scoo[2]=pow(10.0,scalel)*scoo[2]; /* selev */
			scoo2=pow(10.0,scalel)*scoo2;     /* sdepth */
			gdel=pow(10.0,scalel)*gdel;
			sdel=pow(10.0,scalel)*sdel;
		}
		if (scalco!=0){
		 	scoo[0]=pow(10.0,scalco)*scoo[0]; /* sx */
			scoo[1]=pow(10.0,scalco)*scoo[1]; /* sy */
		 	gcoo[0]=pow(10.0,scalco)*gcoo[0]; /* gx */
			gcoo[1]=pow(10.0,scalco)*gcoo[1]; /* gy */			
		}
	    }
	}
	
	/* this s a suboptimum solution ...  */
	*parvec[0]=ns;
	*parvec[1]=dt;
	*parvec[2]=traceid;
	*parvec[3]=tracl;
	*parvec[4]=traceno;
	*parvec[5]=shotno;
	*parvec[6]=globrecno;
	*parvec[7]=comp;
	*parvec[8]=cdp;
	*parvec[9]=gcoo[0];
	*parvec[10]=gcoo[1];
	*parvec[11]=gcoo[2];
	*parvec[12]=scoo[0];
	*parvec[13]=scoo[1];
	*parvec[14]=scoo[2];
	*parvec[15]=scoo2;
	*parvec[16]=offset;
	*parvec[17]=gdel;
	*parvec[18]=sdel;
	*parvec[19]=swdep;
	*parvec[20]=gwdep;
	if (susegy==1){
		*parvec[21]=xcoo;
		*parvec[22]=ycoo;
		*parvec[23]=inlineno;
		*parvec[24]=crosslineno;
		*parvec[25]=sno;
	}
	else{
		pint=(int *)(&d1);
		*parvec[21]=*pint;
		pint=(int *)(&f1);
		*parvec[22]=*pint;
		pint=(int *)(&d2);
		*parvec[23]=*pint;
		pint=(int *)(&f2);
		*parvec[24]=*pint;		
		*parvec[25]=ntraces;	
	}	
	return n;
}

int DDN_wtraceh(FILE * outstream, int lbendian, int ieeeibm, int meterfeet, int susegy, int ns, float dt, float dtshift, float srctime, int traceno, int globrecno, int shotno, int recno, int comp, int trid, int ntraces, int **recpos, int **recpos_loc, float ** srcpos, int nsrc, float xcoo, float ycoo){
	/* sructure from segy.h is not used in order to access variable data portion independently */
	/* however, using segy.h probably makes the code faster */

	extern float DX, DY, DZ, REFREC[4];
	
	float fnull=0.0, fh=0.0, xshot=0.0, yshot=0.0, zshot=0.0, y, z, scalefac;
	//float x;
	short snull=0, sone=1, sh=0, scalef;
	unsigned short /*uh=0,*/ udt, uns;
	int inull=0, ione=1, ih=0, n=0, m, scoo[3], gcoo[3], gwdep, gdel, tracl, offset, scale, doswap=0;
	
	if (lbendian!=LITTLEBIG) doswap=1;
	
	if (nsrc==1){ /* source coordinates are written into trace header fields in case of a single source point per shot */
		xshot=srcpos[1][1];	
		yshot=srcpos[2][1];	
		zshot=srcpos[3][1];
	}
	if (trid==-1) { /* extended definitions of trid in SEG-Y rev. 1 */
		if (susegy==1){
			switch (comp){
				case 0: trid=11; break; /* P (pressure) */
				case 1: trid=12; break; /* Y (vertical) */
				case 2: trid=14; break; /* X (in-line    = horizontal 1 for su) */
				case 3: trid=13; break; /* Z (cross-line = horizontal 2 for su) */
				default: trid=1;
			}		
		}
		else { /* extended definitions of trid in SU */
			switch (comp){
				case 0: trid=1; break;
				case 1: trid=43; break; /* Y (vertical) */
				case 2: trid=44; break; /* X (horizontal 1) */
				case 3: trid=45; break; /* Z (horizontal 2) */
				default: trid=1;
			}
		}
	}
	/* if (trid!=1) fprintf(curerr,"Caution [DDN_wtraceh]: trace ident. code.\n",trid); */
	tracl=recpos_loc[4][recno];
        //x=recpos[1][tracl]*DX-REFREC[1]; /* not in use*/
        y=recpos[2][tracl]*DY-REFREC[2];
        z=recpos[3][tracl]*DZ-REFREC[3];
	scale=3;
	scalefac=pow(10.0,scale);
	if (susegy==1) scalef=-myisign(scale)*(int)pow(10.0,abs(scale)); else scalef=-scale;
	if (meterfeet==1){
		gcoo[0]=(int)floor((recpos[1][tracl]*DX*scalefac)/0.3048+0.5);
		gcoo[2]=(int)floor((recpos[2][tracl]*DY*scalefac)/0.3048+0.5); /* SOFI3D swaps Y and Z */
		gcoo[1]=(int)floor((recpos[3][tracl]*DZ*scalefac)/0.3048+0.5);
		scoo[0]=(int)floor((xshot*scalefac)/0.3048+0.5);
		scoo[2]=(int)floor((yshot*scalefac)/0.3048+0.5); /* SOFI3D swaps Y and Z */
		scoo[1]=(int)floor((zshot*scalefac)/0.3048+0.5);
		gwdep=(int)floor((sqrt(z*z+y*y)*scalefac)/0.3048+0.5);  /* misused parameter for tunnel geometry */
		if (nsrc==1) offset=(int)floor((sqrt(
	                 (recpos[1][recpos_loc[4][traceno]]*DX-xshot)*(recpos[1][recpos_loc[4][traceno]]*DX-xshot)
	                +(recpos[2][recpos_loc[4][traceno]]*DY-yshot)*(recpos[2][recpos_loc[4][traceno]]*DY-yshot)
			+(recpos[3][recpos_loc[4][traceno]]*DZ-zshot)*(recpos[3][recpos_loc[4][traceno]]*DZ-zshot))*scalefac)/0.3048+0.5);		
	} 
	else {
		gcoo[0]=(int)floor(recpos[1][tracl]*DX*scalefac+0.5);
		gcoo[2]=(int)floor(recpos[2][tracl]*DY*scalefac+0.5); /* SOFI3D swaps Y and Z */
		gcoo[1]=(int)floor(recpos[3][tracl]*DZ*scalefac+0.5);
		scoo[0]=(int)floor(xshot*scalefac+0.5);
		scoo[2]=(int)floor(yshot*scalefac+0.5); /* SOFI3D swaps Y and Z */
		scoo[1]=(int)floor(zshot*scalefac+0.5);			
		gwdep=(int)floor(sqrt(z*z+y*y)*scalefac+0.5); /* misused parameter for tunnel geometry */	
		if (nsrc==1) offset=(int)floor((sqrt(
	                 (recpos[1][recpos_loc[4][traceno]]*DX-xshot)*(recpos[1][recpos_loc[4][traceno]]*DX-xshot)
	                +(recpos[2][recpos_loc[4][traceno]]*DY-yshot)*(recpos[2][recpos_loc[4][traceno]]*DY-yshot)
			+(recpos[3][recpos_loc[4][traceno]]*DZ-zshot)*(recpos[3][recpos_loc[4][traceno]]*DZ-zshot))*scalefac)+0.5);		
	}
	gdel=(int)floor(atan2(-y,z)*180.0*scalefac/PI+0.5); /* misused parameter for tunnel geometry */
	/* In case of multiple sources offset does not make much sense, does it? */
	udt=(unsigned short)floor(dt*1.0e6+0.5);
	uns=(unsigned short)ns;
	if (doswap) {
		swap4(&gwdep);
		swap4(&gdel);
		for (m=0;m<3;m++) swap4(&scoo[m]);
		for (m=0;m<3;m++) swap4(&gcoo[m]);
		swap2(&sone);
		swap4(&ione);
		swap2(&scalef);
		swap4(&offset);
		swap2((short *) &udt);
		swap2((short *) &uns);
	}
	/*   1-  4 int tracl;    */ ih=tracl;  if (doswap) swap4(&ih);
	   n+=fwrite(&ih,1,4,outstream); /* trace sequence number within line --- FOR SEISMOGRAMS: RECEIVER NO. FOR THIS SHOT */
	/*   5-  8 int tracr;    */ ih=traceno;    if (doswap) swap4(&ih);
	   n+=fwrite(&ih,1,4,outstream); /* trace sequence number within reel (jth trace in file) */
	/*   9- 12 int fldr;     */ ih=shotno;     if (doswap) swap4(&ih);
	   n+=fwrite(&ih,1,4,outstream); /* field record number --- FOR SEISMOGRAMS: CURENT SHOT NO. */
	/*  13- 16 int tracf;    */ ih=globrecno;  if (doswap) swap4(&ih);
	   n+=fwrite(&ih,1,4,outstream); /* trace number within field record --- RECEIVER NO. */
	/*  17- 20 int ep;       */ ih=comp;       if (doswap) swap4(&ih); 
	   n+=fwrite(&ih,1,4,outstream); /* energy source point no. --- FOR SEISMOGRAMS: component, WARNING: ORIGINAL MEANING IS SOURCE POINT! */
	/*  21- 24 int cdp;      */ ih=tracl;  if (doswap) swap4(&ih);
	   n+=fwrite(&ih,1,4,outstream); /* CDP ensemble number --- FOR SEISMOGRAMS: RECEIVER NO. FOR THIS SHOT */
	/*  25- 28 int cdpt;     */ 
	   n+=fwrite(&ione,1,4,outstream);  /* trace number within CDP ensemble */
	/*  29- 30 short trid;   */ sh=trid;       if (doswap) swap2(&sh); 
	   n+=fwrite(&sh,1,2,outstream); /* trace identification code:
			1 = seismic data (SEISMOGRAM)
			2 = dead (RECEIVER OUTSIDE MODEL)
			3 = dummy (SOURCE and possibly receiver OUTSIDE MODEL)
			4 = time break (SOURCE TIMES)
			6 = sweep (SOURCE SIGNATURE)
		FOR SU DATA: CWP id flags:
			43 = Seismic Data, Vertical Component 
			44 = Seismic Data, Horizontal Component 1 
			45 = Seismic Data, Horizontal Component 2 
			46 = Seismic Data, Radial Component
			47 = Seismic Data, Transverse Component  
		FOR SEG-Y DATA:
			11 = Seismic pressure sensor  
			12 = Multicomponent seismic sensor - Vertical component
			13 = Multicomponent seismic sensor - Cross-line component
			14 = Multicomponent seismic sensor - In-line component
			15 = Rotated multicomponent seismic sensor - Vertical component
			16 = Rotated multicomponent seismic sensor - Transverse component
			17 = Rotated multicomponent seismic sensor - Radial component */
	/*  31- 32 short nvs;    */ n+=fwrite(&sone,1,2,outstream); /* number of vertically summed traces (see vscode in bhed structure) */
	/*  33- 34 short nhs;    */ n+=fwrite(&sone,1,2,outstream); /* number of horizontally summed traces (see vscode in bhed structure) */
	/*  35- 36 short duse;   */ n+=fwrite(&snull,1,2,outstream); /* data use: 1 = production, 2 = test */
	/*  37- 40 int offset;   */ n+=fwrite(&offset,1,4,outstream); /* distance from source point to receiver group 
	                                 (CAUTION: should be negative if opposite to direction in which the line was shot) */
	/*  41- 44 int gelev;    */ n+=fwrite(&gcoo[2],1,4,outstream); /* receiver group elevation from sea level (above sea level is positive)
									  --- CAUTION: SEISMOGRAMS: Y rec. coordinate */
	/*  45- 48 int selev;    */ n+=fwrite(&scoo[2],1,4,outstream); /* source elevation from sea level (above sea level is positive)
									  --- CAUTION: SEISMOGRAMS: Y source coordinate  */
	/*  49- 52 int sdepth;   */ n+=fwrite(&scoo[2],1,4,outstream); /* source depth (positive) --- SEISMOGRAMS: Y source coordinate  */
	/*  53- 56 int gdel;     */ n+=fwrite(&gdel,1,4,outstream);    /* datum elevation at receiver group 
									  --- WARNING: SEISMOGRAMS: DIRECTIVITY NON-STANDARD */
	/*  57- 60 int sdel;     */ n+=fwrite(&inull,1,4,outstream);    /* datum elevation at source */
	/*  61- 64 int swdep;    */ n+=fwrite(&inull,1,4,outstream); /* water depth at source */
	/*  65- 68 int gwdep;    */ n+=fwrite(&gwdep,1,4,outstream); /* water depth at receiver group 
									--- WARNING: SEISMOGRAMS: DIRECTIVITY NON-STANDARD */
	/*  69- 70 short scalel; */ n+=fwrite(&scalef,1,2,outstream); /* scale factor for previous 7 entries with value plus or minus 10 to the
			                                               power 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) */
	/*  71- 72 short scalco; */ n+=fwrite(&scalef,1,2,outstream); /* scale factor for next 4 entries with value plus or minus 10 to the
			                                               power 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) */
	/*  73- 76 int  sx;      */ n+=fwrite(&scoo[0],1,4,outstream); /* X source coordinate */
	/*  77- 80 int  sy;      */ n+=fwrite(&scoo[1],1,4,outstream); /* Y source coordinate */
	/*  81- 84 int  gx;      */ n+=fwrite(&gcoo[0],1,4,outstream); /* X group coordinate  */
	/*  85- 88 int  gy;      */ n+=fwrite(&gcoo[1],1,4,outstream); /* Y group coordinate  */
	/*  89- 90 short counit; */ n+=fwrite(&sone,1,2,outstream); /* coordinate units code: for previous four entries
				1 = length (meters or feet),
				2 = seconds of arc (X = longitude; Y = latitude,
				    a positive value designates the number of seconds east of Greenwich	or north of the equator) */
	/*  91- 92 short wevel;  */ n+=fwrite(&snull,1,2,outstream); /* weathering velocity */
	/*  93- 94 short swevel; */ n+=fwrite(&snull,1,2,outstream); /* subweathering velocity */
	/*  95- 96 short sut;    */ n+=fwrite(&snull,1,2,outstream); /* uphole time at source */
	/*  97- 98 short gut;    */ n+=fwrite(&snull,1,2,outstream); /* uphole time at receiver group */
	/*  99-100 short sstat;	 */ n+=fwrite(&snull,1,2,outstream); /* source static correction */
	/* 101-102 short gstat;	 */ n+=fwrite(&snull,1,2,outstream); /* group static correction */
	/* 103-104 short tstat;	 */ n+=fwrite(&snull,1,2,outstream); /* total static applied */
	sh=(int)-floor(dtshift*1e6+0.5); if (doswap) swap2(&sh); /* time break = strat of the modeling = 0 */
	/* 105-106 short laga;   */ n+=fwrite(&sh,1,2,outstream);
			/* lag time A, time in ms between end of 240-
			   byte trace identification header and time
			   break, positive if time break occurs after
			   end of header, time break is defined as
			   the initiation pulse which maybe recorded
			   on an auxiliary trace or as otherwise
			   specified by the recording system */
	sh=(int)floor(srctime*1e6+0.5); if (doswap) swap2(&sh);		   
	/* 107-108 short lagb;   */ n+=fwrite(&sh,1,2,outstream);
			/* lag time B, time in ms between the time break
			   and the initiation time of the energy source,
			   may be positive or negative */
	sh=(int)floor((dtshift-srctime)*1e6+0.5); if (doswap) swap2(&sh);		   
	/* 109-110 short delrt;  */ n+=fwrite(&sh,1,2,outstream);
			/* delay recording time, time in ms between
			   initiation time of energy source and time
			   when recording of data samples begins
			   (for deep water work if recording does not
			   start at zero time) */
	/* 111-112 short muts;   */ n+=fwrite(&snull,1,2,outstream); /* mute time--start */
	/* 113-114 short mute;   */ n+=fwrite(&snull,1,2,outstream); /* mute time--end */
	/* 115-116 unsigned short ns; */ n+=fwrite(&uns,1,2,outstream); /* number of samples in this trace */
	/* 117-118 unsigned short dt; */ n+=fwrite(&udt,1,2,outstream); /* sample interval; in micro-seconds */
	/* 119-120 short gain;	 */ n+=fwrite(&snull,1,2,outstream); /* gain type of field instruments code:
				1 = fixed,2 = binary, 3 = floating point, 4 ---- N = optional use */
	/* 121-122 short igc;	 */ n+=fwrite(&snull,1,2,outstream); /* instrument gain constant */
	/* 123-124 short igi;	 */ n+=fwrite(&snull,1,2,outstream); /* instrument early or initial gain */
	/* 125-126 short corr;	 */ n+=fwrite(&sone,1,2,outstream); /* correlated: 1 = no, 2 = yes */
	/* 127-128 short sfs;	 */ n+=fwrite(&snull,1,2,outstream); /* sweep frequency at start */
	/* 129-130 short sfe;	 */ n+=fwrite(&snull,1,2,outstream); /* sweep frequency at end */
	/* 131-132 short slen;	 */ n+=fwrite(&snull,1,2,outstream); /* sweep length in ms */
	/* 133-134 short styp;	 */ n+=fwrite(&snull,1,2,outstream); /* sweep type code:
				1 = linear, 2 = parapolic,   3 = exponential, 4 = other */
	/* 135-136 short stas;   */ n+=fwrite(&snull,1,2,outstream); /* sweep trace taper length at start in ms */
	/* 137-138 short stae;   */ n+=fwrite(&snull,1,2,outstream); /* sweep trace taper length at end in ms */
	/* 139-140 short tatyp;  */ n+=fwrite(&snull,1,2,outstream); /* sweep taper type code : 1=linear, 2=cos^2, 3=other */
	/* 141-142 short afilf;  */ n+=fwrite(&snull,1,2,outstream); /* alias filter frequency if used */
	/* 143-144 short afils;  */ n+=fwrite(&snull,1,2,outstream); /* alias filter slope */
	/* 145-146 short nofilf; */ n+=fwrite(&snull,1,2,outstream); /* notch filter frequency if used */
	/* 147-148 short nofils; */ n+=fwrite(&snull,1,2,outstream); /* notch filter slope */
	/* 149-150 short lcf;	 */ n+=fwrite(&snull,1,2,outstream); /* low cut frequency if used */
	/* 151-152 short hcf;	 */ n+=fwrite(&snull,1,2,outstream); /* high cut frequncy if used */
	/* 153-154 short lcs;	 */ n+=fwrite(&snull,1,2,outstream); /* low cut slope */
	/* 155-156 short hcs;	 */ n+=fwrite(&snull,1,2,outstream); /* high cut slope */
	/* 157-158 short year;   */ n+=fwrite(&snull,1,2,outstream); /* year data recorded */
	/* 159-160 short day;	 */ n+=fwrite(&snull,1,2,outstream); /* day of year */
	/* 161-162 short hour;   */ n+=fwrite(&snull,1,2,outstream); /* hour of day (24 hour clock) */
	/* 163-164 short minute; */ n+=fwrite(&snull,1,2,outstream); /* minute of hour */
	/* 165-166 short sec;	 */ n+=fwrite(&snull,1,2,outstream); /* second of minute */
	/* 167-168 short timbas; */ n+=fwrite(&snull,1,2,outstream); /* time basis code: 1 = local, 2 = GMT, 3 = other */
	/* 169-170 short trwf;	 */ n+=fwrite(&snull,1,2,outstream); /* trace weighting factor, i.e. 1/2^N  volts for the least sigificant bit */
	/* 171-172 short grnors; */ n+=fwrite(&snull,1,2,outstream); /* geophone group number of roll switch  position one */
	/* 173-174 short grnofr; */ n+=fwrite(&snull,1,2,outstream); /* geophone group number of trace one within original field record */
	/* 175-176 short grnlof; */ n+=fwrite(&snull,1,2,outstream); /* geophone group number of last trace within original field record */
	/* 177-178 short gaps;	 */ n+=fwrite(&snull,1,2,outstream); /* gap size (total number of groups dropped) */
	/* 179-180 short otrav;	 */ n+=fwrite(&snull,1,2,outstream); /* overtravel taper code: 1 = down (or behind), 2 = up (or ahead) */
	if (susegy==1){
	 ih=(int)floor(xcoo*scalefac+0.5); if (doswap) swap4(&ih); /* 1st horizontal coordinate for models (in SOFI3D ususally X) */
	 /* 181-184 int ???;      */ n+=fwrite(&ih,1,4,outstream);   /* X coordinate of ensemble (CDP) position of this trace (scalar in Trace
									 Header bytes 71-72 applies). The coordinate reference system should be
									 identified through an extended header Location Data stanza. */
	 ih=(int)floor(ycoo*scalefac+0.5); if (doswap) swap4(&ih); /* 2nd horizontal coordinate for models (in SOFI3D ususally Z) */
	 /* 185-188 int ???;      */ n+=fwrite(&ih,1,4,outstream);   /* Y coordinate of ensemble (CDP) position of this trace (scalar in Trace
									 Header bytes 71-72 applies). The coordinate reference system should be
									 identified through an extended header Location Data stanza. */
	 /* 189-192 int ???;      */ n+=fwrite(&inull,1,4,outstream);   /* For 3-D poststack data, this field should be used for the in-line
	 								 number. If one in-line per SEG Y file is being recorded, this value 
									 should be the same for all traces in the file and the same value will 
									 be recorded in bytes 3205-3208 of the Binary File Header. */
	 /* 193-196 int ???;      */ n+=fwrite(&inull,1,4,outstream);   /* For 3-D poststack data, this field should be used for the cross-line
									 number.  This will typically be the same value as the ensemble (CDP)
									 number in Trace Header bytes 21-24, but this does not have to be the 
									 case. */
	 ih=shotno; if (doswap) swap4(&ih);
	 /* 197-200 int ???;      */ n+=fwrite(&ih,1,4,outstream);     /* Shotpoint number - This is probably only applicable to 2-D poststack 
									 data. - Note that it is assumed that the shotpoint number refers to the
									 source location nearest to the ensemble (CDP) location for a particular
									 trace.  If this is not the case, there should be a comment in the Textual
									 File Header explaining what the  shotpoint number actually refers to. */
	 /* 201-202 short ???;    */ n+=fwrite(&snull,1,2,outstream);   /* Scalar to be applied to the shotpoint number in Trace Header bytes 
									 197-200 to  give the real value. If positive, scalar is used as 
									 multiplier; if negative as a divisor; if zero the shotpoint number is not
									 scaled (i.e. it is an integer. A typical value will be -10, allowing
									 shotpoint numbers with one decimal digit to the right of the decimal
									 point). */
									 
	 /* 203-204 short ???;    */ n+=fwrite(&snull,1,2,outstream);   /* Trace value measurement unit:  
									 -1 = Other (should be described in Data Sample Measurement Units Stanza) 
								          0 = Unknown    
									  1 = Pascal (Pa)    
									  2 = Volts (v)  
									  3 = Millivolts (mV) 
  									  4 = Amperes (A) 
 									  5 = Meters (m) 
									  6 = Meters per second (m/s) 
									  7 = Meters per second squared (m/s2) 
									  8 = Newton (N) 
									  9 = Watt (W) 	*/
									  
	 /* 205-208 int ???;    */ n+=fwrite(&inull,1,4,outstream);     /* see bytes 209-210! */
	 /* 209-210 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Transduction Constant - The multiplicative constant used to convert the 
									 Data  Trace samples to the Transduction Units (specified in Trace Header
									 bytes 211- 212).  The constant is encoded as a four-byte, two's 
									 complement integer (bytes 205-208) which is the mantissa and a two-byte,
									 two's complement integer (bytes 209-210) which is the power of ten
									 exponent (i.e. Bytes 205-208 * 10**Bytes  209-210).  */
									 
	 /* 211-212 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Transduction Units - The unit of measurement of the Data Trace samples
									 after  they have been multiplied by the Transduction Constant specified 
									 in Trace  Header bytes 205-210.
									 0 = Unknown    
									 1 = Pascal (Pa)
									 2 = Volts (v)
									 3 = Millivolts (mV) 
									 4 = Amperes (A) 
									 5 = Meters (m) 
									 6 = Meters per second (m/s) 
									 7 = Meters per second squared (m/s2)
									 8 = Newton (N) 
									 9 = Watt (W)  */
	 /* 213-214 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Device/Trace Identifier ? The unit number or id number of the device
									 associated  with the Data Trace (i.e. 4368 for vibrator serial number 
									 4368 or 20316 for gun 16 on string 3 on vessel 2).  This field allows
									 traces to be associated across trace ensembles independently of the 
									 trace number (Trace Header bytes 25-28). */
	 /* 215-216 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Scalar to be applied to times specified in Trace Header bytes 95-114 to
									 give the  true time value in milliseconds.  Scalar = 1, +10, +100, +1000,
									 or +10,000.  If positive, scalar is used as a multiplier; if negative,
									 scalar is used as divisor.  A  value of zero is assumed to be a scalar
									 value of 1. */
									 
	 /* 217-218 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Source Type/Orientation ? Defines the type and the orientation of the
									 energy  source.  The terms vertical, cross-line and in-line refer to the
									 three axes of an  orthogonal coordinate system.  The absolute azimuthal
									 orientation of the  coordinate system axes can be defined in the Bin Grid
									 Definition Stanza.  
									 -1 to -n = Other (should be described in Source Type/Orientation stanza) 
									 0 = Unknown 
									 1 = Vibratory - Vertical orientation 	
									 2 = Vibratory - Cross-line orientation
									 3 = Vibratory - In-line orientation
									 4 = Impulsive - Vertical orientation 
									 5 = Impulsive - Cross-line orientation 
									 6 = Impulsive - In-line orientation 
									 7 = Distributed Impulsive - Vertical orientation 
									 8 = Distributed Impulsive - Cross-line orientation 
									 9 = Distributed Impulsive - In-line orientation */
 	 /* 219-224 ?!?!? ???;  */ n+=fwrite(&snull,1,2,outstream);
		    		  n+=fwrite(&inull,1,4,outstream);     /* Source Energy Direction with respect to the source orientation  - The
									 positive  orientation direction is defined in Bytes 217-218 of the Trace
									 Header.  The energy direction is encoded in tenths of degrees 
									 (i.e. 347.8� is encoded as 3478). */ 							 
	 /* 225-228 int ???;    */ n+=fwrite(&inull,1,4,outstream);     /* see bytes 229-230! */
	 /* 229-230 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Measurement - Describes the source effort used to generate the
									 trace.   The measurement can be simple, qualitative measurements such as
									 the total  weight of explosive used or the peak air gun pressure or the
									 number of vibrators  	times the sweep duration.  Although these simple
									 measurements are acceptable,  it is preferable to use true measurement
									 units of energy or work.  The constant is encoded as a four-byte, two's
									 complement integer (bytes 225-228) which is the mantissa and a two-byte,
									 two's complement integer (bytes 209-230) which is the power of ten
									 exponent (i.e. Bytes 225-228 * 10**Bytes 229- 230). */
									 
	 /* 231-232 short ???;  */ n+=fwrite(&snull,1,2,outstream);     /* Source Measurement Unit - The unit used for the Source Measurement, 
	 								   Trace header bytes 225-230.  
									 -1 = Other (should be described in Source Measurement Unit stanza)   
									  0 = Unknown    
									  1 = Joule (J)
									  2 = Kilowatt (kW)  
									  3 = Pascal (Pa) 
									  4 = Bar (Bar) 
									  4 = Bar-meter (Bar-m) 
									  5 = Newton (N) 
									  6 = Kilograms (kg) */
  
	 /* 233-240 int[2] ???;    */ n+=fwrite(&inull,1,4,outstream); 
	 			     n+=fwrite(&inull,1,4,outstream);  /* Unassigned ? For optional information. */
	}
	else { /* local SU assignments */
		/* 181-184 float d1;     */ fh=dt; native2float(&fh, 1, lbendian, ieeeibm, 0);
	   	   			    n+=fwrite(&fh,1,4,outstream); /* sample spacing for non-seismic data */
		/* 185-188 float f1;     */ n+=fwrite(&fnull,1,4,outstream); /* first sample location for non-seismic data */
		/* 189-192 float d2;     */ n+=fwrite(&fnull,1,4,outstream); /* sample spacing between traces */
		/* 193-296 float f2;     */ n+=fwrite(&fnull,1,4,outstream); /* first trace location */
		/* 197-200 float ungpow; */ n+=fwrite(&fnull,1,4,outstream); /* negative of power used for dynamic range compression */
		/* 201-204 float unscale;*/ n+=fwrite(&fnull,1,4,outstream); /* reciprocal of scaling factor to normalize range */
		/* 205-208 int ntr;      */ ih=ntraces;  if (doswap) swap4(&ih);
		   n+=fwrite(&ih,1,4,outstream); /* number of traces --- NUMBER OF TRACE IN FILE */
		/* 209-210 short mark;   */ n+=fwrite(&snull,1,2,outstream); /* mark selected traces */
       		/* 211-212 short shortpad; */ n+=fwrite(&snull,1,2,outstream); /* alignment padding */
		/* 213-240 short unass[14] */ for (m=0;m<14;m++) n+=fwrite(&snull,1,2,outstream); /* unassigned */
	}
	if (n!=240) fprintf(curerr,"Warning [wtraceh]: segy-header consists of %d bytes instead of 240?!\n",n);
	return n;
}

	
/* reading and writing binary data */

int DDN_rbindata(FILE * instream, int inlen, float * outdata, int outlen, int first, int step, int padding, int lbendian, int ieeeibm, int meterfeet){
	/* CAUTION: first sample is zero, not one !!! */
	float * indata;
	int l, n, m=0;
	if (curerr==NULL) curerr=FP;
	if (instream==NULL) return 0;
	if (outdata==NULL) if ((outdata=(float *)malloc(4*outlen))==NULL) return 0;
	if ((indata=(float *)malloc(4*inlen))==NULL) {free(outdata); outdata=NULL; return 0;}
	if (fread(indata, 4, inlen, instream)!=inlen){free(outdata); outdata=NULL; free(indata); return 0;}
	float2native(indata, inlen, lbendian, ieeeibm, meterfeet);
	if (!step) {
		fprintf(curerr,"Warning [rdata]: reading every 0th sample 0 times?! -> reading all samples once.");
		step=1;
	}
	if (step>0){ /* read every (step)th sample */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	else { /* if (step<0): read every sample -step times */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}		
	free(indata);
	return m; /* samples read to data vector */
}

int DDN_wbindata(FILE * outstream, int outlen, float * indata, int inlen, int first, int step, int padding, int lbendian, int ieeeibm, int meterfeet){
	/* CAUTION: first sample is zero, not one !!! */
	float * outdata;
	int l, n, m=0;
	if (curerr==NULL) curerr=FP;
	if (((outdata=(float *)malloc(4*outlen))==NULL)||(outstream==0)||(indata==NULL)) return 0;
	if (!step) {
		fprintf(curerr,"Warning [rdata]: writing every 0th sample 0 times?! -> writing all samples once.");
		step=1;
	}	
	if (step>0){ /* write every (step)th sample */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	else { /* if (step<0): write every sample -step times */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	/* if (m!=outlen) fprintf(curerr,"Warning [wbindata]: Only %d samples available for output insteasd of %d?!\n",m,outlen); */
	native2float(outdata, outlen, lbendian, ieeeibm, meterfeet);
	m=fwrite(outdata, 4, outlen, outstream);
	free(outdata);
	return m; /* samples written to stream */
}

/* reading and writing textual data */

int DDN_rtxtdata(FILE * instream, int inlen, float * outdata, int outlen, int first, int step, int padding, int asciiebcdic, int meterfeet, char * dataformat, char * dataseperator, int * readdata){
	/* CAUTION: first sample is zero, not one !!! */
	float * indata;
	int l, n, m=0;
	char formatstring[STRING_SIZE],dseperator[]="\t", ddataformat[]="%e";
	if (readdata==NULL) return 0;
	if (dataseperator==NULL) dataseperator=dseperator;
	if (dataformat==NULL) dataformat=ddataformat;
	sprintf(formatstring,"%s%s",dseperator,ddataformat);	
	if (instream==NULL) return 0;
	if (outdata==NULL) if ((outdata=(float *)malloc(4*outlen))==NULL) return 0;
	if ((indata=(float *)malloc(4*inlen))==NULL) {free(outdata); outdata=NULL; return 0;}
	if (ASCIIEBCDIC==asciiebcdic) {
		if ((l=fscanf(instream,dataformat,indata[0]))!=EOF) {
			m=l;
			for (n=1;n<inlen;n++) { 
				l=fscanf(instream,formatstring,indata[n]);
				if (l==EOF) {
					free(indata);
					*readdata=m;
					break;
				}
				else m+=l;
			}
		}
		else {
			free(indata);
			*readdata=0;
			return 0;
		}
	} else  {
		fscanf(curerr,"Error [wtxtdata]: cannot read data in non-native format. \n");
		free(indata);
		*readdata=0;
		return 0;
	}
	inlen=m;	
	if (meterfeet) ft2m(indata, inlen);
	if (!step) {
		fprintf(curerr,"Warning [rtxtdata]: reading every 0th sample 0 times?! -> reading all samples once.");
		step=1;
	}
	if (step>0){ /* read every (step)th sample */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	else { /* if (step<0): read every sample -step times */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	free(indata);
	return m; /* samples read to data vector */
}

int DDN_wtxtdata(FILE * outstream, int outlen, float * indata, int inlen, int first, int step, int padding, int asciiebcdic, int meterfeet, char * dataformat, char * dataseperator, char * dataendmark){
	/* CAUTION: first sample is zero, not one !!! */
	float * outdata;
	char formatstring[STRING_SIZE], dseperator[]="\t", ddataendmark[]="\n", ddataformat[]="%e", tmpoutstring[STRING_SIZE];
	int l, n, m=0;
	if (curerr==NULL) curerr=FP;
	if (dataseperator==NULL) dataseperator=dseperator;
	if (dataendmark==NULL) dataendmark=ddataendmark;
	if (dataformat==NULL) dataformat=ddataformat;
	sprintf(formatstring,"%s%s",dseperator,ddataformat);	
	if (((outdata=(float *)malloc(4*outlen))==NULL)||(outstream==0)||(indata==NULL)) return 0;
	if (!step) {
		fprintf(curerr,"Warning [rdata]: writing every 0th sample 0 times?! -> writing all samples once.");
		step=1;
	}	
	if (step>0){ /* write every (step)th sample */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n+=step) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n+=step) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	else { /* if (step<0): write every sample -step times */
	    if (padding) {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[0];
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=indata[inlen];	
	    }
	    else {
		for (n=first;(m<outlen)&&(n<0);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=0;
		for (n=mymax(n,first);(n<inlen)&&(m<outlen);n++) for (l=0;(m<outlen)&&(l<-step);l++) outdata[m++]=indata[n];
		for (;m<outlen;) outdata[m++]=0;
	    }
	}
	/* if (m!=outlen) fprintf(curerr,"Warning [wtxtdata]: Only %d samples available for output insteasd of %d?!\n",m,outlen); */
	if (meterfeet) m2ft(outdata, outlen);
	if (ASCIIEBCDIC==asciiebcdic) {
		m=fprintf(outstream,dataformat,outdata[0]);
		for (n=1;n<outlen;n++) m+=fprintf(outstream,formatstring,outdata[n]);
		m+=fprintf(outstream,dataendmark);
	} else if ((ASCIIEBCDIC==0)||(asciiebcdic==1)){
		l=sprintf(tmpoutstring,dataformat,outdata[0]);
		m=fwrite(a2ev(tmpoutstring,l),1,l,outstream);
		for (n=1;n<outlen;n++) {
			l=sprintf(tmpoutstring,dataformat,outdata[n]);
			m+=fwrite(a2ev(tmpoutstring,l),1,l,outstream);
		}
		l=sprintf(tmpoutstring,dataendmark);
		m+=fwrite(a2ev(tmpoutstring,l),1,l,outstream);
	} else if ((ASCIIEBCDIC==1)||(asciiebcdic==0)){ /* It would be very surprising if this works! */
		l=sprintf(tmpoutstring,dataformat,outdata[0]);
		m=fwrite(e2av(tmpoutstring,l),1,l,outstream);
		for (n=1;n<outlen;n++) {
			l=sprintf(tmpoutstring,dataformat,outdata[n]);
			m+=fwrite(e2av(tmpoutstring,l),1,l,outstream);
		}
		l=sprintf(tmpoutstring,dataendmark);
		m+=fwrite(e2av(tmpoutstring,l),1,l,outstream);
	} else  {
		fprintf(curerr,"Error [wtxtdata]: unknown character-set?!\n");
		m=0;
	}
	free(outdata);
	return m; /* charcaters written to stream */
}
	
/* reading and writing SEG-Y headers */

int DDN_rsegytxth(FILE * instream, int * asciiebcdic, FILE * outstream){
	/* this subroutine tests some formal features of the textual header, only */
	char segyhbuf[3201], line[81];
	int rin, j, cl=0, el=0;
	if (curerr==NULL) curerr=FP;
	if (instream==NULL) return 0;
	if (ASCIIEBCDIC) fprintf(curerr,"Warning [wsegybinh]: flag for native character format is not ASCII?! (Ignoring this flag ...)\n");
	rin=fread(segyhbuf,1,3200,instream);
	line[80]='\0';	
	if (*asciiebcdic==1) {
		if ((segyhbuf[0]=='c')||(segyhbuf[0]=='C')) {
			*asciiebcdic=0;
			fprintf(curerr,"Caution [rsegytxth]: wrong charset=EBCDIC flag?! Now set to ASCII.\n");
		}
		else for (j=0;j<3201;j++) segyhbuf[j]=e2a(segyhbuf[j]);
	}
	else if ((segyhbuf[0]==a2e('c'))||(segyhbuf[0]==a2e('C'))) {
		*asciiebcdic=1;
		fprintf(curerr,"Caution [rsegytxth]: wrong charset=ASCII flag?! Now set to EBCDIC.\n");
		for (j=0;j<3201;j++) segyhbuf[j]=e2a(segyhbuf[j]);
	}
	for (j=1;j<=40;j++) {
		if ((segyhbuf[(j-1)*80]!='C')&&(segyhbuf[(j-1)*80]!='c')) cl=1;
		if ((segyhbuf[j*80-1]!='\n')) {el=1; segyhbuf[j*80-1]='\n';}
	}
	if (cl) fprintf(curerr,"Caution [rsegytxth]: textual header line(s) do not start with a C?!\n");
	if (el) fprintf(curerr,"Caution [rsegytxth]: textual header line(s) not ended by end-of-line(s)?!\n");
	strncpy(line,&segyhbuf[3120],80);
	if ((!strncasecmp(line,"C40 END TEXTUAL HEADER                                                         \n",79))&&
	    (!strncasecmp(line,"C40 END EBCDIC                                                                 \n",79)))
		fprintf(curerr,"Caution [rsegytxth]: textual header neither ended by END TEXTUAL HEADER nor END EBCDIC tag?!\n");
	if (outstream!=NULL) if (fwrite(segyhbuf,1,3200,outstream)!=3200) 
		fprintf(curerr,"Caution [rsegytxth]: could not write textual header.");
	return rin;
}

int DDN_wsegytxth(FILE * outstream, int asciiebcdic, char * kindofdata, char * infilename, int ns, float dt, int ndt){
	/* this subroutine wirtes an almost empty textual header - some entries starting with "SOFI3D-" are non-standard */
	char segyhbuf[3201];
	int j;
	if (curerr==NULL) curerr=FP;
	if (outstream==NULL) return 0;
	if (ASCIIEBCDIC) fprintf(curerr,"Caution [wsegybinh]: flag for native character format is not ASCII?!\n");	
	/* SEGY-textual header scheme (empty) */
	sprintf(&segyhbuf[0]   ,"%-79.79s\n","C 1 CLIENT                        COMPANY                       CREW NO        ");
	sprintf(&segyhbuf[80]  ,"%-79.79s\n","C 2 LINE            AREA                        MAP ID                         ");
	sprintf(&segyhbuf[160] ,"%-79.79s\n","C 3 REEL NO           DAY-START OF REEL     YEAR      OBSERVER                 ");
	sprintf(&segyhbuf[240] ,"%-79.79s\n","C 4 INSTRUMENT: MFG            MODEL            SERIAL NO                      ");
	sprintf(&segyhbuf[320] ,"%-79.79s\n","C 5 DATA TRACES/RECORD        AUXILIARY TRACES/RECORD         CDP FOLD         ");
	sprintf(&segyhbuf[400] ,"%-79.79s\n","C 6 SAMPLE INTERNAL         SAMPLES/TRACE       BITS/IN      BYTES/SAMPLE      ");
	sprintf(&segyhbuf[480] ,"%-79.79s\n","C 7 RECORDING FORMAT        FORMAT THIS REEL        MEASUREMENT SYSTEM         ");
	sprintf(&segyhbuf[560] ,"%-79.79s\n","C 8 SAMPLE CODE: FLOATING PT     FIXED PT     FIXED PT-GAIN     CORRELATED     ");
	sprintf(&segyhbuf[640] ,"%-79.79s\n","C 9 GAIN  TYPE: FIXED     BINARY     FLOATING POINT     OTHER                  ");
	sprintf(&segyhbuf[720] ,"%-79.79s\n","C10 FILTERS: ALIAS     HZ  NOTCH     HZ  BAND    -     HZ  SLOPE    -    DB/OCT");
	sprintf(&segyhbuf[800] ,"%-79.79s\n","C11 SOURCE: TYPE            NUMBER/POINT        POINT INTERVAL                 ");
	sprintf(&segyhbuf[880] ,"%-79.79s\n","C12     PATTERN:                           LENGTH        WIDTH                 ");
	sprintf(&segyhbuf[960] ,"%-79.79s\n","C13 SWEEP: START     HZ  END     HZ  LENGTH      MS  CHANNEL NO     TYPE       ");
	sprintf(&segyhbuf[1040],"%-79.79s\n","C14 TAPER: START LENGTH       MS  END LENGTH       MS  TYPE                    ");
	sprintf(&segyhbuf[1120],"%-79.79s\n","C15 SPREAD: OFFSET        MAX DISTANCE        GROUP INTERVAL                   ");
	sprintf(&segyhbuf[1200],"%-79.79s\n","C16 GEOPHONES: PER GROUP     SPACING     FREQUENCY     MFG          MODEL      ");
	sprintf(&segyhbuf[1280],"%-79.79s\n","C17     PATTERN:                           LENGTH        WIDTH                 ");
	sprintf(&segyhbuf[1360],"%-79.79s\n","C18 TRACES SORTED BY: RECORD     CDP     OTHER                                 ");
	sprintf(&segyhbuf[1440],"%-79.79s\n","C19 AMPLITUDE RECOVEY: NONE      SPHERICAL DIV       AGC    OTHER              ");
	sprintf(&segyhbuf[1520],"%-79.79s\n","C20 MAP PROJECTION                      ZONE ID       COORDINATE UNITS         ");
	sprintf(&segyhbuf[1600],"%-79.79s\n","C21 PROCESSING: ");
	sprintf(&segyhbuf[1680],"%-79.79s\n","C22 PROCESSING: ");
	/* SOFI3D-header */
	sprintf(&segyhbuf[1760],"C23 SOFI3D-DATA: %s",kindofdata); sprintf(&segyhbuf[1760],"%-79.79s\n",&segyhbuf[1760]);
	sprintf(&segyhbuf[1840],"C24 SOFI3D-INPUT: %s",infilename); sprintf(&segyhbuf[1840],"%-79.79s\n",&segyhbuf[1840]);
	sprintf(&segyhbuf[1920],"C25 SOFI3D-PAR: ns: %d, dt: %d * %e s",ns,ndt,dt); sprintf(&segyhbuf[1920],"%-79.79s\n",&segyhbuf[1920]);
	sprintf(&segyhbuf[2000],"%-79.79s\n","C26 SOFI3D-CONTACT: Thomas Bohlen <thomas.bohlen@kit.edu>      ");
	sprintf(&segyhbuf[2080],"%-79.79s\n","C27 SOFI3D-CONTACT: Thomas Bohlen <thomas.bohlen@kit.edu>        ");
	sprintf(&segyhbuf[2160],"%-79.79s\n","C28 SOFI3D-CONTACT: Thomas Bohlen <thomas.bohlen@kit.edu>       ");
	/* completing SEGY-textual header with almost empty lines, rev.-tag and end-tag */
	for (j=29; j<39;j++) sprintf(&segyhbuf[(j-1)*80],"C%2d%-76.76s\n",j," ");	
	sprintf(&segyhbuf[3040],"%-79.79s\n","C39 SEG Y REV1");
	sprintf(&segyhbuf[3120],"%-79.79s\n","C40 END TEXTUAL HEADER");
	if (asciiebcdic==1) for (j=0;j<3201;j++) segyhbuf[j]=a2e(segyhbuf[j]);
		/* 1: ebcdic textual header, 0 (and default): ascii textual header, */
	return fwrite(segyhbuf,1,3200,outstream);	
}

int DDN_rsegybinh(FILE * instream, int * lbendian, int * ieeeibm, int * meterfeet, int * ns, float * dt){ /* reads only some of the header words */
	int n=0;
	int i;
	int intnull;
	short sintnull, smf, ssh;
	unsigned short usintnull, usdt, usns;
	int doswap=0; if (*lbendian!=LITTLEBIG) doswap=1;
	if (curerr==NULL) curerr=FP;
	   n+=fread(&intnull,1,4,instream); /* job ident no. */
	   n+=fread(&intnull,1,4,instream); /* line no.*/
	   n+=fread(&intnull,1,4,instream); /* reel no. */
	   n+=fread(&sintnull,1,2,instream); /* data traces per ensamble*/
	   n+=fread(&sintnull,1,2,instream); /* dauxillary traces per ensamble. */
	n+=fread(&usdt,1,2,instream); /* sampling interval [mu s] */
	  n+=fread(&usintnull,1,2,instream); /* original sampling interval [mu s] */	
	n+=fread(&usns,1,2,instream); /* no. of samples */
	  n+=fread(&usintnull,1,2,instream); /* no. of samples in the original recording */
	n+=fread(&sintnull,1,2,instream); /* data sample format code (1: 4-byte IBM floats, 5: 4-byte IEEE floats) */  
	if (doswap) swap2(&sintnull);
	if (sintnull==1) *ieeeibm=1;  
	else if (sintnull==5) *ieeeibm=0;
	else {
		if (sintnull>=256) {
			fprintf(curerr,"Caution [rsegybinh]: data format flag > 255: Correct endian?!\n");
			if (LITTLEBIG==1) fprintf(curerr,"\t flag for machine is BIG ENDIAN.\n");
			else if (LITTLEBIG==0) fprintf(curerr,"\t flag for machine is LITTLE ENDIAN.\n");
			else  fprintf(curerr,"\t flag for machine %d is unknown?!\n",LITTLEBIG);
			if (*lbendian==1) fprintf(curerr,"\t flag for data is BIG ENDIAN.\n");
			else if (*lbendian==0) fprintf(curerr,"\t flag for data is LITTLE ENDIAN.\n");
			else  fprintf(curerr,"\t flag for data %d is unknown?!\n",*lbendian);
		}
		if (sintnull==256) {
			fprintf(curerr,"Caution [rsegybinh]: data format flag %d unknown, assuming IBM data with wrong endian.\n",sintnull);
			*lbendian=!*lbendian;
			doswap=!doswap;
			*ieeeibm=1; 
		}
		else if (sintnull==1280) {
			fprintf(curerr,"Caution [rsegybinh]: data format flag %d unknown, assuming IEEE data with wrong endian.\n",sintnull);	
			*lbendian=!*lbendian;			
			doswap=!doswap;
			*ieeeibm=0;				
		}		
		else {
			fprintf(curerr,"Caution [rsegybinh]: data format flag %d unknown/not implemented, assuming IEEE data.\n",sintnull);
			*ieeeibm=0;
		}
	}
	   n+=fread(&sintnull,1,2,instream); /* ensemble fold */
	   n+=fread(&sintnull,1,2,instream); /* trace sorting code 0=unknown */
	   n+=fread(&sintnull,1,2,instream); /* vertical sum code*/
	   for (i=0;i<7;i++) n+=fread(&sintnull,1,2,instream); /* 7 sweep parameters */
	   n+=fread(&sintnull,1,2,instream); /* taper type */
	   n+=fread(&sintnull,1,2,instream); /* corellated data traces: 1=no, 2=yes */
	   n+=fread(&sintnull,1,2,instream); /* binary gain recovered:  1=yes, 2=no */
	   n+=fread(&sintnull,1,2,instream); /* amplitude recovery method:  1=none */
	n+=fread(&smf,1,2,instream); /*2: feet, 1: meter (SOFI3D assumes that SI units are used) */
	   n+=fread(&sintnull,1,2,instream); /* impulse signal polarity pressure */
	   /* not defined: increase in pressure or downward movement gives positive number on tape */
	   /* 1 = increase in pressure or upward movement gives negative number on tape */
	   /* 2 = increase in pressure or upward movement gives positive number on tape */		
	   n+=fread(&sintnull,1,2,instream); /* vibrator polarity code */
	   for (i=0;i<60;i++) n+=fread(&intnull,1,4,instream); /* unassingned */
	  n+=fread(&usintnull,1,2,instream); /* SEG-Y rev. coded as (1byte).(1byte) */
	  n+=fread(&sintnull,1,2,instream); /* all traces consist of the same no. samples at the same sampling interval */
	 n+=fread(&ssh,1,2,instream); /* number of 3200-byte extended textual headers following the binary header */
	   n+=fread(&sintnull,1,2,instream); /* unassingned */
	   for (i=0;i<23;i++) n+=fread(&intnull,1,4,instream); /* unassingned */	   
	if (doswap) {
		*dt= ((float)*swap2((short*)&usdt))*1.0e-6;
		*ns= *swap2((short*)&usns);
		*meterfeet= *swap2(&smf);
		swap2((short*)&ssh);
	}
	else	{
		*dt=((float)usdt)*1.0e-6;
		*ns=usns;
		*meterfeet=smf;
	}	
	if (ssh) {
		if (fseek(instream, 3200*ssh, SEEK_CUR)==(3200*ssh)) 
			fprintf(curerr,"Caution [rsegybinh]: skipping %d * 3200 bytes additional textual header.\n",ssh);
		else fprintf(curerr,"Error [rsegybinh]: could not skip %d * 3200 bytes additional textual header.\n",ssh);
	}
	switch (*meterfeet){
		case 0 :
		   fprintf(curerr,"Caution [rsegybinh]: key for units not set, assuming meter.\n"); 
		   *meterfeet=0; 
		   break;
		case 1 :
		   *meterfeet=0;
		   break;
		case 2 : 
		   *meterfeet=1;
		   break;
		default: 
		   fprintf(curerr,"Caution [rsegybinh]: key for units %d unknown, assuming meter.\n",*meterfeet);
		   *meterfeet=0;
		   break;
	}
	return n+ssh*3200;
}

int DDN_wsegybinh(FILE * outstream, int lbendian, int ieeeibm, int meterfeet, int ns, float dt, int ndt, int ntrpr, int nart){ /* almost empty binary header */
	int n=0;
	int i;
	short h;
	unsigned short uh;
	const int intnull=0;
	const short sintnull=0;
	int doswap=0; if (lbendian!=LITTLEBIG) doswap=1;
	if (curerr==NULL) curerr=FP;
	   n+=fwrite(&intnull,1,4,outstream); /* job ident no. */
	   n+=fwrite(&intnull,1,4,outstream); /* line no. */
	   n+=fwrite(&intnull,1,4,outstream); /* reel no. */
	if (ntrpr==0) fprintf(curerr,"Caution [wsegybinh]: Number of traces (ntrpr) in binary SEG-Y header is zero.");
	h=ntrpr; if (doswap) swap2(&h);
	 n+=fwrite(&h,1,2,outstream); /* data traces per ensamble */
	h=nart;	if (doswap) swap2(&h);
	 n+=fwrite(&h,1,2,outstream); /* auxillary traces per ensamble. */
	uh=(unsigned short)floor(0.5+((float)ndt*dt)*1.0e6); if (doswap) swap2((short*)&uh);
	 n+=fwrite(&uh,1,2,outstream); /* sampling interval [mu s] */
	uh=(unsigned short)floor(0.5+dt*1.0e6); if (doswap) swap2((short*)&uh);
	 n+=fwrite(&uh,1,2,outstream); /* original sampling interval [mu s] */		
	uh=(unsigned short)ns; if (doswap) swap2((short*)&uh);
	 n+=fwrite(&uh,1,2,outstream); /* no. of samples */
	uh=(unsigned short)floor(0.5+(float)dt/(float)ndt); if (doswap) swap2((short*)&uh);
	 n+=fwrite(&uh,1,2,outstream); /* no. of samples in the original recording */
	if (ieeeibm==1) h=1; else h=5; if (doswap) swap2(&h);	
	 n+=fwrite(&h,1,2,outstream); /* data sample format code (1: 4-byte IBM floats, 5: 4-byte IEEE floats) */
	   n+=fwrite(&sintnull,1,2,outstream); /* ensemble fold */
	   n+=fwrite(&sintnull,1,2,outstream); /* trace sorting code 0=unknown */
	   n+=fwrite(&sintnull,1,2,outstream); /* vertical sum code*/
	   for (i=0;i<7;i++) n+=fwrite(&sintnull,1,2,outstream); /* 7 sweep parameters */
	   n+=fwrite(&sintnull,1,2,outstream); /* taper type */
	   n+=fwrite(&sintnull,1,2,outstream); /* corellated data traces: 1=no, 2=yes */
	   n+=fwrite(&sintnull,1,2,outstream); /* binary gain recovered:  1=yes, 2=no */
	   n+=fwrite(&sintnull,1,2,outstream); /* amplitude recovery method:  1=none */
	if (meterfeet==1) h=2; /* feet */ 
	else h=1; if (doswap) swap2(&h); /* meter (SOFI3D assumes that SI units are used) */
	 n+=fwrite(&h,1,2,outstream);
	   n+=fwrite(&sintnull,1,2,outstream); /* impulse signal polarity pressure */
	   /* not defined: increase in pressure or downward movement gives positive number on tape */
	   /* 1 = increase in pressure or upward movement gives negative number on tape */
	   /* 2 = increase in pressure or upward movement gives positive number on tape */		
	   n+=fwrite(&sintnull,1,2,outstream); /* vibrator polarity code */
	   for (i=0;i<60;i++) n+=fwrite(&intnull,1,4,outstream); /* unassingned */
	if (LITTLEBIG==1) uh=256; /* SEG-Y rev. 1.0 : 256*byte1+byte2*/ else uh=1; /* SEG-Y rev. 1.0 256*byte2+byte1 */
	n+=fwrite(&uh,1,2,outstream); /* SEG-Y rev. coded as (1byte).(1byte) */
	h=1; if (doswap) swap2(&h);
	 n+=fwrite(&h,1,2,outstream); /* all traces consist of the same no. samples at the same sampling interval */
	   n+=fwrite(&sintnull,1,2,outstream); /* number of 3200-byte extended textual headers following the binary header */
	   n+=fwrite(&sintnull,1,2,outstream); /* unassingned */
	   for (i=0;i<23;i++) n+=fwrite(&intnull,1,4,outstream); /* unassingned */
	return n;
}



/* converting feet into meter and vice versa */

float * m2ft(float * iovec, int len){
	int n;
	/* 1ft = 0.3048 m (exactly!) */
	for (n=0;n<len;n++)  iovec[n]=iovec[n]/0.3048;
	return iovec;
}

float * ft2m(float * iovec, int len){
	int n;
	/* 1ft = 0.3048 m (exactly!) */
	for (n=0;n<len;n++)  iovec[n]=iovec[n]*0.3048;
	return iovec;
}


/* swapping bytes */

short * swap2(short * ip) { /* swap bytes of any 2-bytes pointer */
	*ip=((*ip>>8)&0xff)|((*ip&0xff)<<8);
	return ip; 
}

int * swap4(int * ip) { /* swap bytes of any 4-bytes pointer */
	*ip=((*ip>>24)&0xff)|((*ip&0xff)<<24)|((*ip>>8)&0xff00)|((*ip&0xff00)<<8);
	return ip; 
}

float * swap4fv(float * fp, int len){/* swap bytes of a vector of 4-bytes floats */
	int n;
	for (n=0;n<len;n++) swap4((int *)&fp[n]);
	return (float *)fp;
}


/* converting any 4-byte floats to native 4-byte floats and vice versa (meter (SI system) = native) */

float * float2native(float * indata, int len, int lbendian, int ieeeibm, int meterfeet){
	int n;
	if (lbendian!=LITTLEBIG) for (n=0;n<len;n++) swap4((int *) &indata[n]);
	if ((!ieeeibm)&&IEEEIBM) ieee2ibmv((int*) indata, len);
	else if ((ieeeibm)&&(!IEEEIBM)) ibm2ieeev((int*) indata, len);
	if (meterfeet) m2ft(indata, len);
	return indata;
}

float * native2float(float * outdata, int len, int lbendian, int ieeeibm, int meterfeet){
	int n;
	if (meterfeet) ft2m(outdata, len);
	if ((!ieeeibm)&&IEEEIBM) ibm2ieeev((int*) outdata, len);
	else if ((ieeeibm)&&(!IEEEIBM)) ieee2ibmv((int*) outdata, len);
	if (lbendian!=LITTLEBIG) for (n=0;n<len;n++) swap4((int *) &outdata[n]);
	return outdata;
}


/* converting ieee floats into ibm floats and vice versa */

int * ieee2ibmv(int * ieee, int n) /* converting a vector of ieee floats into ibm floats */
{
	int h, m, i, d;
	if (curerr==NULL) curerr=FP;
	if (sizeof(int)!=4) {
		fprintf(curerr,"Error [ieee2ibm]: sizeof(int)!=4. \n");	return NULL;}
	if (sizeof(float)!=4) {
		fprintf(curerr,"Error [ieee2ibm]: sizeof(float)!=4. \n"); return NULL;}			
	for (i=0;i<n;++i) {
		h=ieee[i];
		if (h) {
	    		m=(0x007fffff&h)|0x00800000;
	    		d=(int)((0x7f800000&h)>>23)-126;
	    		while (d&0x3) {m>>=1; d++;}
	    		h=(0x80000000&h)|(((d>>2)+64)<< 24)|m;
		}
		ieee[i]=h;
	}
	return ieee;
}


int * ibm2ieeev(int * ibm, int n) /* converting a vector of ibm floats into ieee floats */
{
	int h, m, i, d;
	if (curerr==NULL) curerr=FP;
	if (sizeof(int)!=4) {
		fprintf(curerr,"Error [ibm2ieee]: sizeof(int)!=4. \n");	return NULL;}
	if (sizeof(float)!=4) {
		fprintf(curerr,"Error [ibm2ieee]: sizeof(float)!=4. \n"); return NULL;}		
	for (i=0;i<n;++i) {
		h=ibm[i];		
	if (h) {
		m=0x00ffffff&h;
		if (m) {	
			d=(int)((0x7f000000&h)>>22)-130;
			while (!(m&0x00800000)) {
				m<<=1;			
				d--;
			}
			if (d>254) h=(0x80000000&h)|0x7f7fffff;
			else if (d<=0) h=0;
			else h=(0x80000000&h)|(d<<23)|(0x007fffff&m);
		}
		else {
			h=0;
			fprintf(curerr,"Caution [ibm2ieee]: data point %d set to zero. \n",i);
		}
	}
	ibm[i]=h;
	}
	return ibm;
}



/* converting first 128 ascii characters to ebcdic and ebcdic to ascii (Caution: There seem to exist differnt tables for ebcdic?!) */

const int a2etab[128]={
/* 0-32 according to http://publib.boulder.ibm.com/infocenter/lnxpcomp/v8v101/index.jsp?topic=/com.ibm.xlf101l.doc/xlflr/asciit.htm */
  0,  1,  2,  3, 55, 45, 46, 47, 22,  5, /*0-9*/
 37, 11, 12, 13, 14, 15, 16, 17, 18, 19, /*10-19*/
 60, 61, 50, 38, 24, 25, 63, 39, /*20-27*/
 28, 29, 30, 31, /*28-31: maybe only similar */
 64, 90,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96, 75, 97, /*32-47*/
240,241,242,243,244,245,246,247,248,249, /*48-57*/
122, 94, 76,126,110,111,124, /*58-64*/
193,194,195,196,197,198,199,200,201, /*65-73*/
209,210,211,212,213,214,215,216,217, /*74-82*/
226,227,228,229,230,231,232,233, /*83-90*/
186,224,187,176,109,121, /*91-96*/
129,130,131,132,133,134,135,136,137, /*97-105*/
145,146,147,148,149,150,151,152,153, /*106-114*/
162,163,164,165,166,167,168,169, /*115-122*/
192, 79,208,161,  7}; /*123-127*/


unsigned char a2e(unsigned char ain){ /* ascii 2 ebcdic */
	int n, c;
	c=(int)(unsigned char)ain;
	for (n=0;n<128;n++) if (n==c) return (unsigned char) a2etab[n];
	return 32; /* default: blanc character */
}

unsigned char e2a(unsigned char ein){ /* ebcdic 2 ascii */
	int n,c;
	c=(int)(unsigned char)ein;
	for (n=0;n<128;n++) if (a2etab[n]==c) return (unsigned char) n;
	return 64; /* default: blanc character (equal to ascii @)*/
}

char * a2ev(char * inout, int len){ /* ascii 2 ebcdic for vectors */
	int m;
	for (m=0;m<len;m++){
		inout[m]=(char)a2e((unsigned char) inout[m]);
	}
	return inout;
}

char * e2av(char * inout, int len){ /* ebcdic 2 ascii for vectors */
	int m;
	for (m=0;m<len;m++){
		inout[m]=(char)e2a((unsigned char) inout[m]);
	}
	return inout;
}
