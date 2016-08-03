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
 * For the averaging of material properties each process requires values
 * at the indices 0 and NX+1 etc. These lie on the neighbouring processes.
 * Thus, they have to be copied which is done by this function.
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void matcopy_acoustic(float *** rho, float *** pi){


	extern int MYID, NX, NY, NZ, INDEX[7];
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;
	extern FILE *FP;

	MPI_Status status;	
	double time1=0.0, time2=0.0;	
	int i, j, k;
	float *** bufferlef_to_rig, *** bufferrig_to_lef;
	float *** buffertop_to_bot, *** bufferbot_to_top;
	float *** bufferfro_to_bac, *** bufferbac_to_fro;


	bufferlef_to_rig = f3tensor(0,NY+1,0,NZ+1,1,2);
	bufferrig_to_lef = f3tensor(0,NY+1,0,NZ+1,1,2);
	buffertop_to_bot = f3tensor(0,NX+1,0,NZ+1,1,2);
	bufferbot_to_top = f3tensor(0,NX+1,0,NZ+1,1,2);
	bufferfro_to_bac = f3tensor(0,NY+1,0,NX+1,1,2);
	bufferbac_to_fro = f3tensor(0,NY+1,0,NX+1,1,2);

	
	if (MYID==0){
		fprintf(FP,"\n\n **Message from matcopy (written by PE %d):",MYID);
		fprintf(FP,"\n Copy material properties at inner boundaries ... \n");
		time1=MPI_Wtime();
	}	


	

	/* top-bottom -----------------------------------------------------------*/	
	
	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){

			/* storage of top of local volume into buffer */
			buffertop_to_bot[i][k][1]  =  rho[1][i][k];
			buffertop_to_bot[i][k][2]  =  pi[1][i][k];
			
		}
	}


	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){

			
			/* storage of bottom of local volume into buffer */
			bufferbot_to_top[i][k][1]  =  rho[NY][i][k];
			bufferbot_to_top[i][k][2]  =  pi[NY][i][k];
			
		}
	}
	
	
	MPI_Bsend(&buffertop_to_bot[0][0][1],(NX+2)*(NZ+2)*2,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&buffertop_to_bot[0][0][1], (NX+2)*(NZ+2)*2,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferbot_to_top[0][0][1],(NX+2)*(NZ+2)*2,MPI_FLOAT,INDEX[4],TAG6,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferbot_to_top[0][0][1], (NX+2)*(NZ+2)*2,MPI_FLOAT,INDEX[3],TAG6,MPI_COMM_WORLD,&status);   
	

	
	
	
	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){

			rho[NY+1][i][k] = buffertop_to_bot[i][k][1];
			pi[NY+1][i][k] = buffertop_to_bot[i][k][2];

			
			
		}
	}

	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){

			rho[0][i][k] = bufferbot_to_top[i][k][1];
			pi[0][i][k] = bufferbot_to_top[i][k][2];

		}
	}

	
	/* left-right -----------------------------------------------------------*/	
	
	
	
	
	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){


			/* storage of left edge of local volume into buffer */
			bufferlef_to_rig[j][k][1] =  rho[j][1][k];
			bufferlef_to_rig[j][k][2] =  pi[j][1][k];

		}
	}


	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){
			/* storage of right edge of local volume into buffer */
			bufferrig_to_lef[j][k][1] =  rho[j][NX][k];
			bufferrig_to_lef[j][k][2] =  pi[j][NX][k];

		}
	}


	
 	MPI_Bsend(&bufferlef_to_rig[0][0][1],(NY+2)*(NZ+2)*2,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferlef_to_rig[0][0][1], (NY+2)*(NZ+2)*2,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferrig_to_lef[0][0][1],(NY+2)*(NZ+2)*2,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferrig_to_lef[0][0][1], (NY+2)*(NZ+2)*2,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	
	
	
	
	

	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){

			rho[j][NX+1][k] = bufferlef_to_rig[j][k][1];
			pi[j][NX+1][k] = bufferlef_to_rig[j][k][2];

		}
	}

	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){
			rho[j][0][k] = bufferrig_to_lef[j][k][1];
			pi[j][0][k] = bufferrig_to_lef[j][k][2];

		
		}
	}	
	
		

	
	
		/* front-back -----------------------------------------------------------*/	
	
	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){


			/* storage of front side of local volume into buffer */
			bufferfro_to_bac[j][i][1]  =  rho[j][i][1];
			bufferfro_to_bac[j][i][2]  =  pi[j][i][1];

			
		}
	}


	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){
			
			/* storage of back side of local volume into buffer */
			bufferbac_to_fro[j][i][1]  =  rho[j][i][NZ];
			bufferbac_to_fro[j][i][2]  =  pi[j][i][NZ];

		}
	}


	MPI_Bsend(&bufferfro_to_bac[0][0][1],(NX+2)*(NY+2)*2,MPI_FLOAT,INDEX[5],TAG3,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferfro_to_bac[0][0][1], (NX+2)*(NY+2)*2,MPI_FLOAT,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferbac_to_fro[0][0][1],(NX+2)*(NY+2)*2,MPI_FLOAT,INDEX[6],TAG4,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferbac_to_fro[0][0][1], (NX+2)*(NY+2)*2,MPI_FLOAT,INDEX[5],TAG4,MPI_COMM_WORLD,&status);




	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){

			rho[j][i][NZ+1] = bufferfro_to_bac[j][i][1];
			pi[j][i][NZ+1] = bufferfro_to_bac[j][i][2];


		}
	}


	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){
			
			rho[j][i][0] = bufferbac_to_fro[j][i][1];
			pi[j][i][0] = bufferbac_to_fro[j][i][2];
			

		}
	}


	if (MYID==0){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}

	free_f3tensor(bufferlef_to_rig,0,NY+1,0,NZ+1,1,2);
	free_f3tensor(bufferrig_to_lef,0,NY+1,0,NZ+1,1,2);
	free_f3tensor(buffertop_to_bot,0,NX+1,0,NZ+1,1,2);
	free_f3tensor(bufferbot_to_top,0,NX+1,0,NZ+1,1,2);
	free_f3tensor(bufferfro_to_bac,0,NY+1,0,NX+1,1,2);
	free_f3tensor(bufferbac_to_fro,0,NY+1,0,NX+1,1,2);



}




















