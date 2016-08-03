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
 *  fd.h - include file for sofi3D
 *
 *  ---------------------------------------------------------------------*/


/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>


#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)

#define PI (3.141592653589793)
#define NPAR 50
//#define STRING_SIZE 74 //previous value, sometimes not enough to handle longer file names
#define STRING_SIZE 256
#define REQUEST_COUNT 6
#define NPROCX_MAX 100
#define NPROCY_MAX 100
#define NPROCZ_MAX 100


/* declaration of functions */


void absorb(float *** absorb_coeff);

void absorb_PML(float *** absorb_coeffx, float *** absorb_coeffy, float *** absorb_coeffz);

void CPML_coeff(float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
		float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
		float * K_y, float * alpha_prime_y, float * a_y, float * b_y,
		float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half,
		float * K_z, float * alpha_prime_z, float * a_z, float * b_z,
		float * K_z_half, float * alpha_prime_z_half, float * a_z_half, float * b_z_half);

void CPML_ini_elastic(int * xb, int * yb, int * zb);

void av_mat(float *** rho, float *** pi, float *** u,
		float *** taus, float *** taup,
		float  *** uipjp, float *** ukpkp, float *** uipkp, float *** tausipjp,
		float  *** tausjpkp, float  *** tausipkp, float  *** rjp, float  *** rkp, float  *** rip );

void av_mat_acoustic(float *** rho, float  *** rjp, float  *** rkp, float  *** rip );

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns);

void checkfd(FILE *fp, float *** prho, float *** ppi, float *** pu,
		float *** ptaus, float *** ptaup, float *peta, float **srcpos, int nsrc, int **recpos, int ntr);

void checkfd_acoustic(FILE *fp, float *** prho, float *** ppi, float **srcpos, int nsrc, int **recpos, int ntr);

void checkfd_rsg(FILE *fp, float *** prho, float *** ppi, float *** pu,
		float *** ptaus, float *** ptaup, float *peta);

void comm_ini(float *** bufferlef_to_rig,
		float *** bufferrig_to_lef, float *** buffertop_to_bot,
		float *** bufferbot_to_top, float *** bufferfro_to_bac,
		float *** bufferbac_to_fro, MPI_Request *req_send, MPI_Request *req_rec);

void comm_ini_s(float *** bufferlef_to_rig,
		float *** bufferrig_to_lef, float *** buffertop_to_bot,
		float *** bufferbot_to_top, float *** bufferfro_to_bac,
		float *** bufferbac_to_fro, MPI_Request *req_send, MPI_Request *req_rec);

void comm_ini_acoustic(float *** bufferlef_to_rig,
		float *** bufferrig_to_lef, float *** buffertop_to_bot,
		float *** bufferbot_to_top, float *** bufferfro_to_bac,
		float *** bufferbac_to_fro, MPI_Request *req_send, MPI_Request *req_rec);


void eqsource(int nt, float *** sxx, float *** syy, float *** szz,
		float *** sxy, float *** syz, float *** sxz,
		float **  srcpos_loc, float ** signals, int nsrc, int * stype,
		float amom, float str, float dip, float rake);

void exchange_par(void);

void exchange_s_rsg(float *** sxx, float *** syy, float *** szz,
		float *** sxy, float *** syz, float *** sxz,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro);

double exchange_s(int nt, float *** sxx, float *** syy, float *** szz,
		float *** sxy, float *** syz, float *** sxz,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro, MPI_Request * req_send, MPI_Request * req_rec);

double exchange_s_acoustic(int nt, float *** sxx,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro, MPI_Request * req_send, MPI_Request * req_rec);

void exchange_v_rsg(int nt, float *** vx, float *** vy, float *** vz,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro);

double exchange_v(int nt, float *** vx, float *** vy, float *** vz,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro, MPI_Request * req_send, MPI_Request * req_rec);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,float *** syz, float *** sxz,
		float *** rxx, float *** ryy,float *** rzz, float *** rxy, float *** ryz, float *** rxz,
		float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y,
		float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z,
		float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy,
		float *** psi_vxz, float *** psi_vyz, float *** psi_vzz);

void info(FILE *fp);

void initproc(void);

int initsour(int nxs,int nys, int  nzs, int *nxsl,int *nysl, int  *nzsl );

void matcopy(float *** rho, float *** pi, float *** u,
		float *** taus, float *** taup);

void matcopy_acoustic(float *** rho, float *** pi);

void merge(int nsnap, int type);

void mergemod(char modfile[STRING_SIZE], int format);

void model_visco(float  ***  rho, float ***  pi, float ***  u,
		float ***  taus, float ***  taup, float *  eta);

void model_elastic(float  ***  rho, float ***  pi, float ***  u,
		float ***  taus, float ***  taup, float *  eta);

void model_acoustic(float  ***  rho, float ***  pi);

void note(FILE *fp);

void outseis(FILE *fp, FILE *fpdata, float **section,
		int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
		int nsrc, int ns, int seis_form[6], int ishot, int comp);

void outseis_glob(FILE *fp, FILE *fpdata, float **section,
		int **recpos, int **recpos_loc, int ntr, float ** srcpos,
		int nsrc, int ns, int seis_form[6], int ishot, int comp);

void output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

int plane_wave(float *** force_points);

void PML_ini_acoustic(int * xa, int * xb, int * ya, int * yb, int * za, int * zb);

void psource(int nt, float *** sxx, float *** syy, float *** szz, float **  srcpos_loc, float ** signals, int nsrc, int * stype);

void psource_acoustic(int nt, float *** sxx, float **  srcpos_loc, float ** signals, int nsrc, int * stype);

void psource_rsg(int nt, float *** sxx, float *** syy, float *** szz,
		float **  srcpos_loc, float ** signals, int nsrc);

float ** pwsources(int *nsrc, int *stype);

float *rd_sour(int *nts,FILE* fp_source);

void readbufs(float *** sxx, float *** syy, float *** szz,
		float *** sxy, float *** syz, float *** sxz,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro);

void readbufv(float *** vx, float *** vy, float *** vz, 
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro);

float readdsk(FILE *fp_in, int format);


void read_par_json(FILE *fp, char *fileinp);

void readmod_acoustic(float  ***  rho, float ***  pi, int ishot);

void readmod(float  ***  rho, float ***  pi, float ***  u, float ***  taus, float ***  taup, float *  eta);

int **receiver(FILE *fp, int *ntr);


void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionvz,
		float **sectionp, float **sectioncurl, float **sectiondiv,
		int  **recpos, int  **recpos_loc, int ntr, float ** srcpos, int nsrc,int ns);

void saveseis_glob(FILE *fp, float **sectiondata,
		int  **recpos, int  **recpos_loc, int ntr, float ** srcpos, int nsrc,int ns, int sectiondatatype);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,
		float *** rxx, float *** ryy,float *** rzz, float *** rxy, float *** ryz, float *** rxz,
		float *** syz, float *** sxz, float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, 
		float *** psi_sxy_y, float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z, 
		float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy, 
		float *** psi_vxz, float *** psi_vyz, float *** psi_vzz);

void seismo_acoustic(int lsamp, int ntr, int **recpos, float **sectionvx, float **sectionvy, 
		float **sectionvz, float **sectiondiv, float **sectioncurl, float **sectionp,
		float ***vx, float ***vy, float ***vz, float ***sxx, float ***pi);

void seismo(int lsamp, int ntr, int **recpos, float **sectionvx, float **sectionvy, 
		float **sectionvz, float **sectiondiv, float **sectioncurl, float **sectionp,
		float ***vx, float ***vy, float ***vz,
		float ***sxx, float ***syy, float ***szz, float ***pi, float ***u);

void seismo_rsg(int lsamp, int ntr, int **recpos, float **sectionvx, float **sectionvy, 
		float **sectionvz, float **sectiondiv, float **sectioncurl, float **sectionp,
		float ***vx, float ***vy, float ***vz,
		float ***sxx, float ***syy, float ***szz, float ***pi, float ***u);

void snap_acoustic(FILE *fp, int nt, int nsnap, int format, int type, 
		float ***vx, float ***vy, float ***vz, float ***sxx, float ***pi,
		int idx, int idy, int idz, int nx1, int ny1, int nz1, int nx2,
		int ny2, int nz2);

void snap(FILE *fp, int nt, int nsnap, int format, int type, 
		float ***vx, float ***vy, float ***vz, float ***sxx, float ***syy, float ***szz, float ***u, float ***pi,
		int idx, int idy, int idz, int nx1, int ny1, int nz1, int nx2,
		int ny2, int nz2);


void snap_rsg(FILE *fp, int nt, int nsnap, int format, int type, 
		float ***vx, float ***vy, float ***vz, float ***sxx, float ***syy, float ***szz,
		float ***u, float ***pi,
		int idx, int idy, int idz, int nx1, int ny1, int nz1, int nx2,
		int ny2, int nz2);

void snapmerge(int nsnap);

float **sources(FILE * fpsrc, int *nsrc, int * stype);

int **splitrec(int **recpos,int *ntr_loc, int ntr,int *recswitch);

float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc, int * stype_loc, int *stype);

void surface(int ndepth, float *** u, float *** pi, float ***taus, float *** taup,
		float * eta, float *** sxx, float ***syy, float ***szz, float *** sxy,float *** syz,
		float *** rxx, float *** ryy, float ***rzz, float *** vx, float *** vy,
		float *** vz, float * K_x, float * a_x, float * b_x, float * K_z, float * a_z, float * b_z, 
		float *** psi_vxx, float *** psi_vzz );

void surface_elastic(int ndepth, float *** u, float *** pi,
		float *** sxx, float ***syy, float ***szz, float *** sxy,float *** syz,
		float *** vx, float *** vy, float *** vz,
		float * K_x, float * a_x, float * b_x, float * K_z, float * a_z, float * b_z, 
		float *** psi_vxx, float *** psi_vzz );

void surface_acoustic(int ndepth,  float *** pi, float *** sxx, float *** vx, float *** vy, float *** vz);

void timing(double * time_v_update,  double * time_s_update, double * time_s_exchange, double * time_v_exchange,
		double * time_timestep, int ishot);


double update_s(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
                float *** vx, float *** vy, float *** vz,
                float *** sxx, float *** syy, float *** szz, float *** sxy,
                float *** syz, float *** sxz, float *** rxx, float *** ryy,
                float *** rzz, float *** rxy, float *** ryz, float *** rxz,
                float ***  pi, float ***  u, float ***  uipjp, float ***  ujpkp, float ***  uipkp,
                float  ***  taus, float  ***  tausipjp, float  ***  tausjpkp, float  ***  tausipkp,
                float  ***  taup, float *  eta,
                float *** vxyyx,float *** vyzzy,float *** vxzzx,float *** vxxyyzz,float *** vyyzz,float *** vxxzz,float *** vxxyy,
                float *** vxyyx_2,float *** vyzzy_2,float *** vxzzx_2,float *** vxxyyzz_2,float *** vyyzz_2,float *** vxxzz_2,float *** vxxyy_2,
                float *** vxyyx_3,float *** vyzzy_3,float *** vxzzx_3,float *** vxxyyzz_3,float *** vyyzz_3,float *** vxxzz_3,float *** vxxyy_3,
                float *** vxyyx_4,float *** vyzzy_4,float *** vxzzx_4,float *** vxxyyzz_4,float *** vyyzz_4,float *** vxxzz_4,float *** vxxyy_4,
                float *** rxx_2, float *** ryy_2,
                float *** rzz_2, float *** rxy_2, float *** ryz_2, float *** rxz_2,float *** rxx_3, float *** ryy_3,
                float *** rzz_3, float *** rxy_3, float *** ryz_3, float *** rxz_3,float *** rxx_4, float *** ryy_4,
                float *** rzz_4, float *** rxy_4, float *** ryz_4, float *** rxz_4);


double update_s_elastic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
                        float *** vx, float *** vy, float *** vz,
                        float *** sxx, float *** syy, float *** szz, float *** sxy,
                        float *** syz, float *** sxz, float *** rxx, float *** ryy,
                        float *** rzz, float *** rxy, float *** ryz, float *** rxz,
                        float ***  pi, float ***  u, float ***  uipjp, float ***  ujpkp, float ***  uipkp,
                        float  ***  taus, float  ***  tausipjp, float  ***  tausjpkp, float  ***  tausipkp,
                        float  ***  taup, float *  eta,
                        float *** vxyyx,float *** vyzzy,float *** vxzzx,float *** vxxyyzz,float *** vyyzz,float *** vxxzz,float *** vxxyy,
                        float *** vxyyx_2,float *** vyzzy_2,float *** vxzzx_2,float *** vxxyyzz_2,float *** vyyzz_2,float *** vxxzz_2,float *** vxxyy_2,
                        float *** vxyyx_3,float *** vyzzy_3,float *** vxzzx_3,float *** vxxyyzz_3,float *** vyyzz_3,float *** vxxzz_3,float *** vxxyy_3,
                        float *** vxyyx_4,float *** vyzzy_4,float *** vxzzx_4,float *** vxxyyzz_4,float *** vyyzz_4,float *** vxxzz_4,float *** vxxyy_4);

/*double update_s_PML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
 float *** vx, float *** vy, float *** vz, float *** sxx, float *** syy, float *** szz, float *** sxy,
 float *** syz, float *** sxz, float *** vx1, float *** vy1, float *** vz1, float *** sxx1, float *** syy1, float *** szz1, float *** sxy1,
 float *** syz1, float *** sxz1, float *** vx2, float *** vy2, float *** vz2, float *** sxx2, float *** syy2, float *** szz2, float *** sxy2,
 float *** syz2, float *** sxz2, float *** vx3, float *** vy3, float *** vz3, float *** sxx3, float *** syy3, float *** szz3, float *** sxy3,
 float *** syz3, float *** sxz3, float *** rxx, float *** ryy, float *** rzz, float *** rxy, float *** ryz, float *** rxz,float ***  pi,
 float ***  u, float ***  uipjp, float ***  ujpkp, float ***  uipkp, float  ***  taus, float  ***  tausipjp, float  ***  tausjpkp,
 float  ***  tausipkp, float  ***  taup, float *  eta, float *** absorb_coeffx, float *** absorb_coeffy, float *** absorb_coeffz);*/

double update_s_CPML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz, float *** rxx, float *** ryy,
		float *** rzz, float *** rxy, float *** ryz, float *** rxz, float ***  pi, float ***  u, float ***  uipjp, float ***  ujpkp, float ***  uipkp,
		float  ***  taus, float  ***  tausipjp, float  ***  tausjpkp, float  ***  tausipkp,  float  ***  taup, float *  eta,
		float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		float * K_z, float * a_z, float * b_z, float * K_z_half, float * a_z_half, float * b_z_half,
		float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy, float *** psi_vxz, float *** psi_vyz, float *** psi_vzz);

double update_s_CPML_elastic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz, float ***  pi, float ***  u, float ***  uipjp,
		float ***  ujpkp, float ***  uipkp, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		float * K_z, float * a_z, float * b_z, float * K_z_half, float * a_z_half, float * b_z_half,
		float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy, float *** psi_vxz, float *** psi_vyz, float *** psi_vzz);

double update_s_acoustic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
		float *** vx, float *** vy, float *** vz, float *** sxx, float ***  pi);

double update_s_acoustic_PML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
		float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** sxx1, float *** sxx2, float *** sxx3,
		float ***  pi, float ***absorb_coeffx, float ***absorb_coeffy, float ***absorb_coeffz);

void update_s_rsg(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt,
		float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,
		float *** syz, float *** sxz, float *** rxx, float *** ryy,
		float *** rzz, float *** rxy, float *** ryz, float *** rxz,
		float ***  pi, float ***  u,
		float  ***  taus, float  ***  taup, float *  eta);

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,
		float *** syz, float *** sxz, float  ***  rho,  float  *** rjp, float  *** rkp, float  *** rip,
        float **  srcpos_loc, float ** signals, int nsrc, float ***absorb_coeff, int * stype, float *** svx, float *** svy, float *** svz,
        float *** svx_2, float *** svy_2, float *** svz_2, float *** svx_3, float *** svy_3, float *** svz_3,
        float *** svx_4, float *** svy_4, float *** svz_4);

/*double update_v_PML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *** vx, float *** vy, float *** vz, 
float *** sxx, float *** syy, float *** szz, float *** sxy,float *** syz, float *** sxz, float *** vx1, float *** vy1, 
float *** vz1, float *** sxx1, float *** syy1, float *** szz1, float *** sxy1, float *** syz1, float *** sxz1, float *** vx2, 
float *** vy2, float *** vz2, float *** sxx2, float *** syy2, float *** szz2, float *** sxy2, float *** syz2, float *** sxz2, 
float *** vx3, float *** vy3, float *** vz3, float *** sxx3, float *** syy3, float *** szz3, float *** sxy3, float *** syz3, 
float *** sxz3,float  ***  rho, float **  srcpos_loc, float ** signals, int nsrc, float *** absorb_coeffx, float *** absorb_coeffy, float *** absorb_coeffz, int * stype);*/

double update_v_CPML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,
		float *** syz, float *** sxz, float  ***  rho,  float  *** rjp, float  *** rkp, float  *** rip,
		float **  srcpos_loc, float ** signals, int nsrc, float *** absorb_coeff, int * stype,
		float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		float * K_z, float * a_z, float * b_z, float * K_z_half, float * a_z_half, float * b_z_half,
		float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y, float *** psi_syy_y,
		float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z);

double update_v_acoustic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float  ***  rho, float **  srcpos_loc, float ** signals, int nsrc, float ***absorb_coeff, int * stype);

double update_v_acoustic_PML(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz, float *** sxx,  float *** vx1,
		float *** vy2, float *** vz3, float  ***  rho, float **  srcpos_loc,
		float ** signals, int nsrc, float ***absorb_coeffx, float ***absorb_coeffy, float ***absorb_coeffz, int * stype);

void update_v_rsg(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
		int nt, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy,
		float *** syz, float *** sxz, float  ***  rho,
		float **  srcpos_loc, float ** signals, int nsrc, float ***absorb_coeff);

float ** wavelet(float **srcpos_loc, int nsrc);

void writebufs(float *** sxx, float *** syy, float *** szz, 
		float *** sxy, float *** syz, float *** sxz,
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro);

void writebufv(float *** vx, float *** vy, float *** vz, 
		float *** bufferlef_to_rig, float *** bufferrig_to_lef,
		float *** buffertop_to_bot, float *** bufferbot_to_top,
		float *** bufferfro_to_bac, float *** bufferbac_to_fro);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float *** rho, int format);

void writepar(FILE *fp, int ns);

void zero_acoustic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *** vx, float *** vy, float *** vz,
		float *** sxx);

void zero_elastic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz,
                  float *** vxyyx,float *** vyzzy,float *** vxzzx,float *** vxxyyzz,float *** vyyzz,float *** vxxzz,float *** vxxyy,
                  float *** vxyyx_2,float *** vyzzy_2,float *** vxzzx_2,float *** vxxyyzz_2,float *** vyyzz_2,float *** vxxzz_2,float *** vxxyy_2,
                  float *** vxyyx_3,float *** vyzzy_3,float *** vxzzx_3,float *** vxxyyzz_3,float *** vyyzz_3,float *** vxxzz_3,float *** vxxyy_3,
                  float *** vxyyx_4,float *** vyzzy_4,float *** vxzzx_4,float *** vxxyyzz_4,float *** vyyzz_4,float *** vxxzz_4,float *** vxxyy_4,
                  float *** svx, float *** svy, float *** svz,
                  float *** svx_2, float *** svy_2, float *** svz_2, float *** svx_3, float *** svy_3, float *** svz_3,
                  float *** svx_4, float *** svy_4, float *** svz_4);

void zero(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *** vx, float *** vy, float *** vz,
                  float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz,
                  float *** vxyyx,float *** vyzzy,float *** vxzzx,float *** vxxyyzz,float *** vyyzz,float *** vxxzz,float *** vxxyy,
                  float *** vxyyx_2,float *** vyzzy_2,float *** vxzzx_2,float *** vxxyyzz_2,float *** vyyzz_2,float *** vxxzz_2,float *** vxxyy_2,
                  float *** vxyyx_3,float *** vyzzy_3,float *** vxzzx_3,float *** vxxyyzz_3,float *** vyyzz_3,float *** vxxzz_3,float *** vxxyy_3,
                  float *** vxyyx_4,float *** vyzzy_4,float *** vxzzx_4,float *** vxxyyzz_4,float *** vyyzz_4,float *** vxxzz_4,float *** vxxyy_4,
                  float *** svx, float *** svy, float *** svz,
                  float *** svx_2, float *** svy_2, float *** svz_2, float *** svx_3, float *** svy_3, float *** svz_3,
          float *** svx_4, float *** svy_4, float *** svz_4, float *** rxx, float *** ryy, float *** rzz,
          float *** rxy, float *** ryz, float *** rxz,float *** rxx_2, float *** ryy_2,
          float *** rzz_2, float *** rxy_2, float *** ryz_2, float *** rxz_2,float *** rxx_3, float *** ryy_3,
          float *** rzz_3, float *** rxy_3, float *** ryz_3, float *** rxz_3,float *** rxx_4, float *** ryy_4,
          float *** rzz_4, float *** rxy_4, float *** ryz_4, float *** rxz_4);

void zero_elastic_CPML( int NX, int NY, int NZ, float *** vx, float *** vy, float *** vz,float *** sxx, float *** syy,
		float *** szz, float *** sxy, float *** syz, float *** sxz, float *** rxx, float *** ryy, float *** rzz,
		float *** rxy, float *** ryz, float *** rxz, float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x,
		float *** psi_sxy_y, float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z,
		float *** psi_szz_z, float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy,
                       float *** psi_vzy, float *** psi_vxz, float *** psi_vyz, float *** psi_vzz,float *** rxx_2, float *** ryy_2,
                       float *** rzz_2, float *** rxy_2, float *** ryz_2, float *** rxz_2,float *** rxx_3, float *** ryy_3,
                       float *** rzz_3, float *** rxy_3, float *** ryz_3, float *** rxz_3,float *** rxx_4, float *** ryy_4,
                       float *** rzz_4, float *** rxy_4, float *** ryz_4, float *** rxz_4);


/* declaration of functions for json parser in json_parser.c*/
int read_objects_from_intputfile(FILE *fp, char input_file[STRING_SIZE],char ** varname_list,char ** value_list);

void print_objectlist_screen(FILE *fp, int number_readobject,char ** varname_list,char ** value_list);

int count_occure_charinstring(char stringline[STRING_SIZE], char teststring[]);

void copy_str2str_uptochar(char string_in[STRING_SIZE], char string_out[STRING_SIZE], char teststring[]);

int get_int_from_objectlist(char string_in[STRING_SIZE], int number_readobject, int * int_buffer,
		char ** varname_list,char ** value_list);

int get_float_from_objectlist(char string_in[STRING_SIZE], int number_readobject, float * double_buffer,
		char ** varname_list,char ** value_list);

int get_string_from_objectlist(char string_in[STRING_SIZE], int number_readobject, char string_buffer[STRING_SIZE],
		char ** varname_list,char ** value_list);

int is_string_blankspace(char string_in[STRING_SIZE]);

void remove_blankspaces_around_string(char string_in[STRING_SIZE] );

void add_object_tolist(char string_name[STRING_SIZE],char string_value[STRING_SIZE], int * number_read_object,
		char ** varname_list,char ** value_list );


/* utility functions (defined in file util.c)*/
void err(char error_text[]);
void err2(char errformat[],char errfilename[]);
void warning(char warn_text[]);

double maximum(float **a, int nx, int ny);

float *vector(int nl, int nh);
float **fmatrix(int nrl, int nrh, int ncl, int nch);

int *ivector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);

float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);

void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl,
		int ndh);

double *dvector(int nl, int nh);
void free_dvector(double *v, int nl, int nh);

