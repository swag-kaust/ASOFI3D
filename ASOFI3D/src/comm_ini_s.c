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
 *   initialisation of repeated comunications. This may reduce the
 *   network overhead. Communication is started later
 *   with MPI_START(request)
  *  ----------------------------------------------------------------------*/

#include "fd.h"

void comm_ini_s(float *** bufferlef_to_rig, float *** bufferrig_to_lef, 
float *** buffertop_to_bot, float *** bufferbot_to_top, 
float *** bufferfro_to_bac, float *** bufferbac_to_fro, 
MPI_Request *req_send, MPI_Request *req_rec){


	extern int NX, NY, NZ, INDEX[7], FDORDER;
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	int nf1, nf2;

	/* comunication initialisation for persistent communication */

	/* buffer arrays are copied into local buffers using buffered send (bsend),
	  Actually send is activated by MPI_START(request) within time loop.
	  MPI_BSEND and MPI_RECV (see below) are then non-blocking.
	*/
	
	
	/* number of wavefield parameters that need to be exchanged - see exchange_v.c */
	nf1=3*FDORDER/2 -1;
	nf2=nf1-1;
	
	
	MPI_Bsend_init(&bufferlef_to_rig[1][1][1],NY*NZ*nf2,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD,&req_send[0]);
	MPI_Bsend_init(&bufferrig_to_lef[1][1][1],NY*NZ*nf1,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD,&req_send[1]);
	MPI_Bsend_init(&bufferfro_to_bac[1][1][1],NX*NY*nf2,MPI_FLOAT,INDEX[5],TAG3,MPI_COMM_WORLD,&req_send[2]);
	MPI_Bsend_init(&bufferbac_to_fro[1][1][1],NX*NY*nf1,MPI_FLOAT,INDEX[6],TAG4,MPI_COMM_WORLD,&req_send[3]);
	MPI_Bsend_init(&buffertop_to_bot[1][1][1],NX*NZ*nf2,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD,&req_send[4]);
	MPI_Bsend_init(&bufferbot_to_top[1][1][1],NX*NZ*nf1,MPI_FLOAT,INDEX[4],TAG6,MPI_COMM_WORLD,&req_send[5]);

	/* initialising of receive of buffer arrays. Same arrays for send and receive
	   are used. Thus, before starting receive, it is necessary to check if
	   Bsend has copied data into local buffers, i.e. has completed.
	*/
	MPI_Recv_init(&bufferlef_to_rig[1][1][1],NY*NZ*nf2,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&req_rec[0]);
	MPI_Recv_init(&bufferrig_to_lef[1][1][1],NY*NZ*nf1,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&req_rec[1]);
	MPI_Recv_init(&bufferfro_to_bac[1][1][1],NX*NY*nf2,MPI_FLOAT,INDEX[6],TAG3,MPI_COMM_WORLD,&req_rec[2]);
	MPI_Recv_init(&bufferbac_to_fro[1][1][1],NX*NY*nf1,MPI_FLOAT,INDEX[5],TAG4,MPI_COMM_WORLD,&req_rec[3]);
	MPI_Recv_init(&buffertop_to_bot[1][1][1],NX*NZ*nf2,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&req_rec[4]);
	MPI_Recv_init(&bufferbot_to_top[1][1][1],NX*NZ*nf1,MPI_FLOAT,INDEX[3],TAG6,MPI_COMM_WORLD,&req_rec[5]);

}
