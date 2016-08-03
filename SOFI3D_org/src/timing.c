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
 *   output timing information, i.e. real times for wavefield updates and exchange
 *    and some statitics (average, standard deviation)
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void timing(double * time_v_update,  double * time_s_update, double * time_s_exchange, double * time_v_exchange,
            double * time_timestep, int ishot){


extern int NT, RUN_MULTIPLE_SHOTS;
extern FILE *FP; 
extern char LOG_FILE[STRING_SIZE];

char T_FILE[STRING_SIZE];
FILE *fp;
int nt;
double 	time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0, time_av_s_exchange=0.0, time_av_timestep=0.0, 
	stdev_v_update=0.0, stdev_s_update=0.0, stdev_s_exchange=0.0, stdev_v_exchange=0.0, stdev_timestep=0.0;



	for (nt=1;nt<=NT;nt++){
		time_av_v_update+=time_v_update[nt];
		time_av_s_update+=time_s_update[nt];
		time_av_s_exchange+=time_s_exchange[nt];
		time_av_v_exchange+=time_v_exchange[nt];
		time_av_timestep+=time_timestep[nt];
	}

	time_av_v_update=time_av_v_update/(double)NT;
	time_av_s_update=time_av_s_update/(double)NT;
	time_av_v_exchange=time_av_v_exchange/(double)NT;
	time_av_s_exchange=time_av_s_exchange/(double)NT;
	time_av_timestep=time_av_timestep/(double)NT;

	/* calculate standard deviation */
	for (nt=1;nt<=NT;nt++){
		stdev_v_update+=pow((time_v_update[nt]-time_av_v_update),2.0);
		stdev_s_update+=pow((time_s_update[nt]-time_av_s_update),2.0);
		stdev_s_exchange+=pow((time_s_exchange[nt]-time_av_s_exchange),2.0);
		stdev_v_exchange+=pow((time_v_exchange[nt]-time_av_v_exchange),2.0);
		stdev_timestep+=pow((time_timestep[nt]-time_av_timestep),2.0);

	}	

	stdev_v_update=sqrt(stdev_v_update/(NT-1))/(time_av_v_update*0.01);
	stdev_s_update=sqrt(stdev_s_update/(NT-1))/(time_av_s_update*0.01);
	stdev_s_exchange=sqrt(stdev_s_exchange/(NT-1))/(time_av_s_exchange*0.01);
	stdev_v_exchange=sqrt(stdev_v_exchange/(NT-1))/(time_av_v_exchange*0.01);
	stdev_timestep=sqrt(stdev_timestep/(NT-1))/(time_av_timestep*0.01);

	fprintf(FP,"**Info from function timing (written by PE 0) \n");
	fprintf(FP," Average times and standard deviations (in percent) for \n");
	fprintf(FP,"   velocity update:  \t %6.3f seconds +/- %6.3f per cent \n",time_av_v_update, stdev_v_update);
	fprintf(FP,"   stress update:  \t %6.3f seconds +/- %6.3f per cent \n",time_av_s_update, stdev_s_update);
	fprintf(FP,"   velocity exchange:  \t %6.3f seconds +/- %6.3f per cent   \n",time_av_v_exchange, stdev_v_exchange);
	fprintf(FP,"   stress exchange:  \t %6.3f seconds +/- %6.3f per cent  \n",time_av_s_exchange, stdev_s_exchange);
	fprintf(FP,"   timestep:  \t\t %6.3f seconds +/- %6.3f per cent  \n",time_av_timestep, stdev_timestep);


	/* output of timings to ASCII file (for later performance analysis) */
	if (RUN_MULTIPLE_SHOTS) 
		sprintf(T_FILE,"%s.timings.shot%d",LOG_FILE,ishot);  
	else 
		sprintf(T_FILE,"%s.timings",LOG_FILE);  
	fprintf(FP," PE 0 is writing timing information in ASCII format to %s \n", T_FILE);
	fp=fopen(T_FILE,"w");
	fprintf(fp," %% timestep \t  velocity update [s] \t stress update [s] \t velocity exchange [s] \t stress exchange [s]\t total [s]\n");
	for (nt=1;nt<=NT;nt++)
		fprintf(fp," %d \t\t %e \t\t %e \t\t %e \t\t %e \t\t %e \n", 
			nt, time_v_update[nt],time_s_update[nt],time_s_exchange[nt],time_v_exchange[nt], time_timestep[nt]);
	fclose(fp);
}
