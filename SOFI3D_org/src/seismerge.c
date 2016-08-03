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
 *  merging SEG-Y files with non-native endian quick and dirty (The SEG-Y header of
 *	first file is preserved, except for number of traces.)
 *  seismerge <segyfile-name_without_cpu-no> <maximum_cpu-no> <output-filename>
 *  ----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

int main(int argc, char **argv){

	char basename[256], *outname, *inname;
	FILE * infiles, * outfile;
	int n, maxcpu=0, m=0, l, b=0, c;
	unsigned short nsamp, msamp, ntraces=0, mtraces; /* mixed short and unsigned short */
	n=argc;

	if (n==1) {
		fprintf(stdout," usage of seismerge (merging SEG-Y files):\n");
		fprintf(stdout," seismerge segyfile-name_without_.cpu-no. [maximum_cpu-no. [output-filename]]\n\n");
		fprintf(stdout," Path length (input and output file) must not exceed 255 characters.\n");		
		fprintf(stdout," Default number for maximum-cpu-no. : 1000.\n");
		return 0;
	}
	strncpy(basename,argv[1],256);
	if (strlen(basename)>255) {
		fprintf(stdout," seismerge cannot handle path with more than 255 characters.\n");
		return 0;
	}
	if (n>2) {
		sscanf(argv[2], "%d", &maxcpu);
		if (maxcpu==0) fprintf(stdout,"Error: seismerge could not read max. number of cpus.\n");
	}
	if (maxcpu==0) fprintf(stdout," Warning: Setting maximum CPU-no. to default value: %d.\n",maxcpu=1000);	
	if (n>3) outname=argv[3]; else outname=basename;	
	if ((outfile=fopen(outname,"w+"))!=NULL)
		fprintf(stderr, "\t File %s opend for writing. ... \n",outname);
	else {
		fprintf(stderr, "\t Outout file %s could not be opened for writing.\n",outname);
		return -1;
	}
	inname=(char *)malloc(strlen(basename)+12);
	n=0;
	sprintf(inname, "%s%s%d", basename,".",n);
	if((infiles=fopen(inname,"r+"))!=NULL){
		fseek(infiles,3220,SEEK_SET);
		fread(&nsamp,2,1,infiles);
		fseek(infiles,3714,SEEK_SET);
		fread(&msamp,2,1,infiles);
		if (nsamp!=msamp)
			fprintf(stdout,"Warning: Number of samples in SEG-Y header and 1st trace header differ?!\n");
		if (nsamp!=0) { /* -> ns from SEG-Y header */
			fseek(infiles,-3600,SEEK_END);
			msamp=nsamp;
			nsamp=((nsamp>>8)&0xff)|((nsamp&0xff)<<8);
			ntraces=ftell(infiles)/(nsamp*4+240);
			fseek(infiles,3714,SEEK_SET);
			fwrite(&msamp,2,1,infiles);
			for (l=1;l<ntraces;l++){
				fseek(infiles,nsamp*4+238,SEEK_CUR);
				fwrite(&msamp,2,1,infiles);
			}					
		} 
		else { /* -> ns from trace headers */
			fseek(infiles,3220,SEEK_SET);
			fwrite(&msamp,2,1,infiles);
			msamp=((msamp>>8)&0xff)|((msamp&0xff)<<8);
			fseek(infiles,0,SEEK_END);
			l=ftell(infiles);
			fseek(infiles,3716,SEEK_SET);
			while (ftell(infiles)<l-(msamp*4+238)){
				fseek(infiles,msamp*4+238,SEEK_CUR);
				fread(&msamp,2,1,infiles);
				msamp=((msamp>>8)&0xff)|((msamp&0xff)<<8);
				ntraces++;
			}
			ntraces++;	
		}
		fseek(infiles,0,SEEK_END);
		rewind(infiles);
		while((c=fgetc(infiles))!=EOF) {
			fputc(c,outfile);
			b++;
		}
		fclose(infiles);
		fprintf(stderr, "\t File %s with %d traces (%d bytes) has been copied to file %s. ... \n",inname,ntraces,b,outname);
		m=1;
	}
	for (n++;n<maxcpu + 1;n++){
		mtraces=0;
		b=0;
		sprintf(inname, "%s%s%d", basename,".",n);
		if((infiles=fopen(inname,"r+"))!=NULL){
			fseek(infiles,3220,SEEK_SET);
			fread(&nsamp,2,1,infiles);
			fseek(infiles,3714,SEEK_SET);
			fread(&msamp,2,1,infiles);
			if (nsamp!=msamp)
				fprintf(stdout,"Warning: Number of samples in SEG-Y header and 1st trace header differ?!\n");
			if (nsamp!=0) { /* -> ns from SEG-Y header */
				fseek(infiles,-3600,SEEK_END);
				msamp=nsamp;
				nsamp=((nsamp>>8)&0xff)|((nsamp&0xff)<<8);
				mtraces=ftell(infiles)/(nsamp*4+240);
				fseek(infiles,3714,SEEK_SET);
				fwrite(&msamp,2,1,infiles);
				for (l=1;l<mtraces;l++){
					fseek(infiles,nsamp*4+238,SEEK_CUR);				
					fwrite(&msamp,2,1,infiles);
				}
			} 
			else { /* -> ns from trace headers */
				msamp=((msamp>>8)&0xff)|((msamp&0xff)<<8);
				fseek(infiles,0,SEEK_END);
				l=ftell(infiles);
				fseek(infiles,3716,SEEK_SET);
				while (ftell(infiles)<l-(msamp*4+238)){
					fseek(infiles,msamp*4+238,SEEK_CUR);
					fread(&msamp,2,1,infiles);
					msamp=((msamp>>8)&0xff)|((msamp&0xff)<<8);
					mtraces++;
				}
				mtraces++;
			}
			ntraces+=mtraces;
			fseek(infiles,3600,SEEK_SET);
			while((c=fgetc(infiles))!=EOF) {
				fputc(c,outfile);
				b++;
			}
			fclose(infiles);
			fprintf(stderr, "\t File %s has been appended to file %s (%d bytes appended, now %d traces). ... \n",inname,outname,b,ntraces);
			m++;
		}
	}
	fseek(outfile,3212,SEEK_SET);
	mtraces=((ntraces>>8)&0xff)|((ntraces&0xff)<<8);
	fwrite(&mtraces,2,1,outfile);
	fseek(outfile,3604,SEEK_SET);
	for (l=1;l<=ntraces;l++){
		n=((l>>24)&0xff)|((l&0xff)<<24)|((l>>8)&0xff00)|((l&0xff00)<<8);
		fwrite(&n,4,1,outfile);
		fseek(outfile,106,SEEK_CUR);	
		fread(&msamp,2,1,outfile);
		msamp=((msamp>>8)&0xff)|((msamp&0xff)<<8);
		fseek(outfile,msamp*4+128,SEEK_CUR);
	}
	fclose(outfile);
	free(inname);
	fprintf(stderr, " %d traces in %d files (with basename %s) \n merged into one file (%s). \n",ntraces,m,basename,outname);
	return 0;
}
