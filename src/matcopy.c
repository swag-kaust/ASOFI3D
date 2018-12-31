/*------------------------------------------------------------------------
 * For the averaging of material properties each process requires values
 * at the indices 0 and NX+1 etc. These lie on the neighbouring processes.
 * Thus, they have to be copied which is done by this function.
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"


#define NUMPARAM 14

// 5 isotropic + 9 Cij - needs to be corrected VK

void matcopy(float *** rho, float *** pi, float *** u,
        float *** C11, float *** C12, float *** C13, float *** C22, float *** C23, float *** C33,
        float *** C44, float *** C55, float *** C66,
		float *** taus, float *** taup){


	extern int MYID, NX, NY, NZ, L, INDEX[7];
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;
	extern FILE *FP;

	MPI_Status status;
	double time1=0.0, time2=0.0;
	int i, j, k;
	float *** bufferlef_to_rig, *** bufferrig_to_lef;
	float *** buffertop_to_bot, *** bufferbot_to_top;
	float *** bufferfro_to_bac, *** bufferbac_to_fro;


	// 12345 - rho,pi,u,taus,taup 678... C11,C12,C13,C22,C23,C33,C44,C55,C66

	bufferlef_to_rig = f3tensor(0,NY+1,0,NZ+1,1,NUMPARAM);
	bufferrig_to_lef = f3tensor(0,NY+1,0,NZ+1,1,NUMPARAM);
	buffertop_to_bot = f3tensor(0,NX+1,0,NZ+1,1,NUMPARAM);
	bufferbot_to_top = f3tensor(0,NX+1,0,NZ+1,1,NUMPARAM);
	bufferfro_to_bac = f3tensor(0,NY+1,0,NX+1,1,NUMPARAM);
	bufferbac_to_fro = f3tensor(0,NY+1,0,NX+1,1,NUMPARAM);


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
			buffertop_to_bot[i][k][3]  =  u[1][i][k];

			// anisotropic (orthorhombic) parameters in alphabetic order

			buffertop_to_bot[i][k][6]  =  C11[1][i][k];
			buffertop_to_bot[i][k][7]  =  C12[1][i][k];
			buffertop_to_bot[i][k][8]  =  C13[1][i][k];
			buffertop_to_bot[i][k][9]  =  C22[1][i][k];
			buffertop_to_bot[i][k][10]  =  C23[1][i][k];
			buffertop_to_bot[i][k][11]  =  C33[1][i][k];
			buffertop_to_bot[i][k][12]  =  C44[1][i][k];
			buffertop_to_bot[i][k][13]  =  C55[1][i][k];
			buffertop_to_bot[i][k][14]  =  C66[1][i][k];


			if (L){
				buffertop_to_bot[i][k][4]  = taus[1][i][k];
				buffertop_to_bot[i][k][5]  = taup[1][i][k];
			}

		}
	}


	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){


			/* storage of bottom of local volume into buffer */
			bufferbot_to_top[i][k][1]  =  rho[NY][i][k];
			bufferbot_to_top[i][k][2]  =  pi[NY][i][k];
			bufferbot_to_top[i][k][3]  =  u[NY][i][k];

			// anisotropic (orthorhombic) parameters in alphabetic order

			bufferbot_to_top[i][k][6]  =  C11[NY][i][k];
			bufferbot_to_top[i][k][7]  =  C12[NY][i][k];
			bufferbot_to_top[i][k][8]  =  C13[NY][i][k];
			bufferbot_to_top[i][k][9]  =  C22[NY][i][k];
			bufferbot_to_top[i][k][10]  =  C23[NY][i][k];
			bufferbot_to_top[i][k][11]  =  C33[NY][i][k];
			bufferbot_to_top[i][k][12]  =  C44[NY][i][k];
			bufferbot_to_top[i][k][13]  =  C55[NY][i][k];
			bufferbot_to_top[i][k][14]  =  C66[NY][i][k];

			if (L){
				bufferbot_to_top[i][k][4]  = taus[NY][i][k];
				bufferbot_to_top[i][k][5]  = taup[NY][i][k];
			}
		}
	}


	MPI_Bsend(&buffertop_to_bot[0][0][1],(NX+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&buffertop_to_bot[0][0][1], (NX+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferbot_to_top[0][0][1],(NX+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[4],TAG6,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferbot_to_top[0][0][1], (NX+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[3],TAG6,MPI_COMM_WORLD,&status);





	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){

			rho[NY+1][i][k] = buffertop_to_bot[i][k][1];
			pi[NY+1][i][k] = buffertop_to_bot[i][k][2];
			u[NY+1][i][k] = buffertop_to_bot[i][k][3];

			// anisotropic (orthorhombic) parameters in alphabetic order
			// note that y direction is vertical and positive downwards
			// top to bottom  - top was loaded at the previous block

            C11[NY+1][i][k] = buffertop_to_bot[i][k][6];
            C12[NY+1][i][k] = buffertop_to_bot[i][k][7];
			C13[NY+1][i][k] = buffertop_to_bot[i][k][8];
			C22[NY+1][i][k] = buffertop_to_bot[i][k][9];
			C23[NY+1][i][k] = buffertop_to_bot[i][k][10];
			C33[NY+1][i][k] = buffertop_to_bot[i][k][11];
			C44[NY+1][i][k] = buffertop_to_bot[i][k][12];
			C55[NY+1][i][k] = buffertop_to_bot[i][k][13];
			C66[NY+1][i][k] = buffertop_to_bot[i][k][14];

			if (L){
				taus[NY+1][i][k] = buffertop_to_bot[i][k][4];
				taup[NY+1][i][k] = buffertop_to_bot[i][k][5];
			}


		}
	}

	for (i=0;i<=NX+1;i++){
		for (k=0;k<=NZ+1;k++){

			rho[0][i][k] = bufferbot_to_top[i][k][1];
			pi[0][i][k] = bufferbot_to_top[i][k][2];
			u[0][i][k] = bufferbot_to_top[i][k][3];

			// anisotropic (orthorhombic) parameters in alphabetic order
			// note that y direction is vertical and positive downwards
			// bottom to top  - bottom was loaded at the previous block

            C11[0][i][k] = bufferbot_to_top[i][k][6];
            C12[0][i][k] = bufferbot_to_top[i][k][7];
			C13[0][i][k] = bufferbot_to_top[i][k][8];
			C22[0][i][k] = bufferbot_to_top[i][k][9];
			C23[0][i][k] = bufferbot_to_top[i][k][10];
			C33[0][i][k] = bufferbot_to_top[i][k][11];
			C44[0][i][k] = bufferbot_to_top[i][k][12];
			C55[0][i][k] = bufferbot_to_top[i][k][13];
			C66[0][i][k] = bufferbot_to_top[i][k][14];


			if (L){
				taus[0][i][k] = bufferbot_to_top[i][k][4];
				taup[0][i][k] = bufferbot_to_top[i][k][5];
			}
		}
	}


	/* left-right -----------------------------------------------------------*/




	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){


			/* storage of left edge of local volume into buffer */
			bufferlef_to_rig[j][k][1] =  rho[j][1][k];
			bufferlef_to_rig[j][k][2] =  pi[j][1][k];
			bufferlef_to_rig[j][k][3] =  u[j][1][k];


			// anisotropic (orthorhombic) parameters in alphabetic order

			bufferlef_to_rig[j][k][6]  =  C11[j][1][k];
			bufferlef_to_rig[j][k][7]  =  C12[j][1][k];
			bufferlef_to_rig[j][k][8]  =  C13[j][1][k];
			bufferlef_to_rig[j][k][9]  =  C22[j][1][k];
			bufferlef_to_rig[j][k][10]  =  C23[j][1][k];
			bufferlef_to_rig[j][k][11]  =  C33[j][1][k];
			bufferlef_to_rig[j][k][12]  =  C44[j][1][k];
			bufferlef_to_rig[j][k][13]  =  C55[j][1][k];
			bufferlef_to_rig[j][k][14]  =  C66[j][1][k];



			if (L){
				bufferlef_to_rig[j][k][4] =  taus[j][1][k];
				bufferlef_to_rig[j][k][5] =  taup[j][1][k];
			}
		}
	}


	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){
			/* storage of right edge of local volume into buffer */
			bufferrig_to_lef[j][k][1] =  rho[j][NX][k];
			bufferrig_to_lef[j][k][2] =  pi[j][NX][k];
			bufferrig_to_lef[j][k][3] =  u[j][NX][k];

			// anisotropic (orthorhombic) parameters in alphabetic order

			bufferrig_to_lef[j][k][6]  =  C11[j][NX][k];
			bufferrig_to_lef[j][k][7]  =  C12[j][NX][k];
			bufferrig_to_lef[j][k][8]  =  C13[j][NX][k];
			bufferrig_to_lef[j][k][9]  =  C22[j][NX][k];
			bufferrig_to_lef[j][k][10]  =  C23[j][NX][k];
			bufferrig_to_lef[j][k][11]  =  C33[j][NX][k];
			bufferrig_to_lef[j][k][12]  =  C44[j][NX][k];
			bufferrig_to_lef[j][k][13]  =  C55[j][NX][k];
			bufferrig_to_lef[j][k][14]  =  C66[j][NX][k];


			if (L){
				bufferrig_to_lef[j][k][4] =  taus[j][NX][k];
				bufferrig_to_lef[j][k][5] =  taup[j][NX][k];
			}
		}
	}



	MPI_Bsend(&bufferlef_to_rig[0][0][1],(NY+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferlef_to_rig[0][0][1], (NY+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferrig_to_lef[0][0][1],(NY+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferrig_to_lef[0][0][1], (NY+2)*(NZ+2)*NUMPARAM,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&status);






	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){

			rho[j][NX+1][k] = bufferlef_to_rig[j][k][1];
			pi[j][NX+1][k] = bufferlef_to_rig[j][k][2];
			u[j][NX+1][k] = bufferlef_to_rig[j][k][3];

			// anisotropic (orthorhombic) parameters in alphabetic order
			// note that y direction is vertical and positive downwards
			// bottom to top  - bottom was loaded at the previous block

            C11[j][NX+1][k] = bufferlef_to_rig[j][k][6];
            C12[j][NX+1][k] = bufferlef_to_rig[j][k][7];
			C13[j][NX+1][k] = bufferlef_to_rig[j][k][8];
			C22[j][NX+1][k] = bufferlef_to_rig[j][k][9];
			C23[j][NX+1][k] = bufferlef_to_rig[j][k][10];
			C33[j][NX+1][k] = bufferlef_to_rig[j][k][11];
			C44[j][NX+1][k] = bufferlef_to_rig[j][k][12];
			C55[j][NX+1][k] = bufferlef_to_rig[j][k][13];
			C66[j][NX+1][k] = bufferlef_to_rig[j][k][14];

			if (L){
				taus[j][NX+1][k] = bufferlef_to_rig[j][k][4];
				taup[j][NX+1][k] = bufferlef_to_rig[j][k][5];
			}
		}
	}

	for (j=0;j<=NY+1;j++){
		for (k=0;k<=NZ+1;k++){
			rho[j][0][k] = bufferrig_to_lef[j][k][1];
			pi[j][0][k] = bufferrig_to_lef[j][k][2];
			u[j][0][k] = bufferrig_to_lef[j][k][3];

			// anisotropic (orthorhombic) parameters in alphabetic order
			// note that y direction is vertical and positive downwards
			// bottom to top  - bottom was loaded at the previous block

            C11[j][0][k] = bufferrig_to_lef[j][k][6];
            C12[j][0][k] = bufferrig_to_lef[j][k][7];
			C13[j][0][k] = bufferrig_to_lef[j][k][8];
			C22[j][0][k] = bufferrig_to_lef[j][k][9];
			C23[j][0][k] = bufferrig_to_lef[j][k][10];
			C33[j][0][k] = bufferrig_to_lef[j][k][11];
			C44[j][0][k] = bufferrig_to_lef[j][k][12];
			C55[j][0][k] = bufferrig_to_lef[j][k][13];
			C66[j][0][k] = bufferrig_to_lef[j][k][14];

			if (L){
				taus[j][0][k] = bufferrig_to_lef[j][k][4];
				taup[j][0][k] = bufferrig_to_lef[j][k][5];
			}

		}
	}





	/* front-back -----------------------------------------------------------*/

	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){


			/* storage of front side of local volume into buffer */
			bufferfro_to_bac[j][i][1]  =  rho[j][i][1];
			bufferfro_to_bac[j][i][2]  =  pi[j][i][1];
			bufferfro_to_bac[j][i][3]  =  u[j][i][1];

			// anisotropic (orthorhombic) parameters in alphabetic order

			bufferfro_to_bac[j][i][6]  =  C11[j][i][1];
			bufferfro_to_bac[j][i][7]  =  C12[j][i][1];
			bufferfro_to_bac[j][i][8]  =  C13[j][i][1];
			bufferfro_to_bac[j][i][9]  =  C22[j][i][1];
			bufferfro_to_bac[j][i][10]  =  C23[j][i][1];
			bufferfro_to_bac[j][i][11]  =  C33[j][i][1];
			bufferfro_to_bac[j][i][12]  =  C44[j][i][1];
			bufferfro_to_bac[j][i][13]  =  C55[j][i][1];
			bufferfro_to_bac[j][i][14]  =  C66[j][i][1];

			if (L){
				bufferfro_to_bac[j][i][4]  = taus[j][i][1];
				bufferfro_to_bac[j][i][5]  = taup[j][i][1];
			}

		}
	}


	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){

			/* storage of back side of local volume into buffer */
			bufferbac_to_fro[j][i][1]  =  rho[j][i][NZ];
			bufferbac_to_fro[j][i][2]  =  pi[j][i][NZ];
			bufferbac_to_fro[j][i][3]  =  u[j][i][NZ];

			// anisotropic (orthorhombic) parameters in alphabetic order

			bufferbac_to_fro[j][i][6]  =  C11[j][i][NZ];
			bufferbac_to_fro[j][i][7]  =  C12[j][i][NZ];
			bufferbac_to_fro[j][i][8]  =  C13[j][i][NZ];
			bufferbac_to_fro[j][i][9]  =  C22[j][i][NZ];
			bufferbac_to_fro[j][i][10]  =  C23[j][i][NZ];
			bufferbac_to_fro[j][i][11]  =  C33[j][i][NZ];
			bufferbac_to_fro[j][i][12]  =  C44[j][i][NZ];
			bufferbac_to_fro[j][i][13]  =  C55[j][i][NZ];
			bufferbac_to_fro[j][i][14]  =  C66[j][i][NZ];

			if (L){
				bufferbac_to_fro[j][i][4]  = taus[j][i][NZ];
				bufferbac_to_fro[j][i][5]  = taup[j][i][NZ];
			}
		}
	}


	MPI_Bsend(&bufferfro_to_bac[0][0][1],(NX+2)*(NY+2)*NUMPARAM,MPI_FLOAT,INDEX[5],TAG3,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferfro_to_bac[0][0][1], (NX+2)*(NY+2)*NUMPARAM,MPI_FLOAT,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferbac_to_fro[0][0][1],(NX+2)*(NY+2)*NUMPARAM,MPI_FLOAT,INDEX[6],TAG4,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferbac_to_fro[0][0][1], (NX+2)*(NY+2)*NUMPARAM,MPI_FLOAT,INDEX[5],TAG4,MPI_COMM_WORLD,&status);




	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){

			rho[j][i][NZ+1] = bufferfro_to_bac[j][i][1];
			pi[j][i][NZ+1] = bufferfro_to_bac[j][i][2];
			u[j][i][NZ+1] = bufferfro_to_bac[j][i][3];

			C11[j][i][NZ+1] = bufferfro_to_bac[j][i][6];
            C12[j][i][NZ+1] = bufferfro_to_bac[j][i][7];
			C13[j][i][NZ+1] = bufferfro_to_bac[j][i][8];
			C22[j][i][NZ+1] = bufferfro_to_bac[j][i][9];
			C23[j][i][NZ+1] = bufferfro_to_bac[j][i][10];
			C33[j][i][NZ+1] = bufferfro_to_bac[j][i][11];
			C44[j][i][NZ+1] = bufferfro_to_bac[j][i][12];
			C55[j][i][NZ+1] = bufferfro_to_bac[j][i][13];
			C66[j][i][NZ+1] = bufferfro_to_bac[j][i][14];

			if (L){
				taus[j][i][NZ+1] = bufferfro_to_bac[j][i][4];
				taup[j][i][NZ+1] = bufferfro_to_bac[j][i][5];
			}

		}
	}


	for (i=0;i<=NX+1;i++){
		for (j=0;j<=NY+1;j++){

			rho[j][i][0] = bufferbac_to_fro[j][i][1];
			pi[j][i][0] = bufferbac_to_fro[j][i][2];
			u[j][i][0] = bufferbac_to_fro[j][i][3];

			C11[j][i][0] = bufferbac_to_fro[j][i][6];
            C12[j][i][0] = bufferbac_to_fro[j][i][7];
			C13[j][i][0] = bufferbac_to_fro[j][i][8];
			C22[j][i][0] = bufferbac_to_fro[j][i][9];
			C23[j][i][0] = bufferbac_to_fro[j][i][10];
			C33[j][i][0] = bufferbac_to_fro[j][i][11];
			C44[j][i][0] = bufferbac_to_fro[j][i][12];
			C55[j][i][0] = bufferbac_to_fro[j][i][13];
			C66[j][i][0] = bufferbac_to_fro[j][i][14];

			if (L){
				taus[j][i][0] = bufferbac_to_fro[j][i][4];
				taup[j][i][0] = bufferbac_to_fro[j][i][5];
			}

		}
	}


	if (MYID==0){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}

	free_f3tensor(bufferlef_to_rig,0,NY+1,0,NZ+1,1,NUMPARAM);
	free_f3tensor(bufferrig_to_lef,0,NY+1,0,NZ+1,1,NUMPARAM);
	free_f3tensor(buffertop_to_bot,0,NX+1,0,NZ+1,1,NUMPARAM);
	free_f3tensor(bufferbot_to_top,0,NX+1,0,NZ+1,1,NUMPARAM);
	free_f3tensor(bufferfro_to_bac,0,NY+1,0,NX+1,1,NUMPARAM);
	free_f3tensor(bufferbac_to_fro,0,NY+1,0,NX+1,1,NUMPARAM);



}




















