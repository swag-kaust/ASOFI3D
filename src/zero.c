/*------------------------------------------------------------------------
 *  Initialization of the wave field with zero values (zero wavefield)
 *
 *  ----------------------------------------------------------------------*/

#include "data_structures.h"
#include "fd.h"
#include "globvar.h"


void zero(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, Velocity *v,
          Tensor3d *s,
          VelocityDerivativesTensor *dv,
          VelocityDerivativesTensor *dv_2,
          VelocityDerivativesTensor *dv_3,
          VelocityDerivativesTensor *dv_4,
          StressDerivativesWrtVelocity *ds_dv,
          StressDerivativesWrtVelocity *ds_dv_2,
          StressDerivativesWrtVelocity *ds_dv_3,
          StressDerivativesWrtVelocity *ds_dv_4,
          Tensor3d *r,
          Tensor3d *r_2,
          Tensor3d *r_3,
          Tensor3d *r_4) {
    
    
    extern int FDORDER_TIME;
    extern int FDORDER;
    register int i, j, k;

    float ***vx = v->x;
    float ***vy = v->y;
    float ***vz = v->z;

    float ***sxx = s->xx;
    float ***syy = s->yy;
    float ***szz = s->zz;
    float ***sxy = s->xy;
    float ***syz = s->yz;
    float ***sxz = s->xz;

    float ***vxyyx   = dv->xyyx;
    float ***vyzzy   = dv->yzzy;
    float ***vxzzx   = dv->xzzx;
    float ***vxxyyzz = dv->xxyyzz;
    float ***vyyzz   = dv->yyzz;
    float ***vxxzz   = dv->xxzz;
    float ***vxxyy   = dv->xxyy;

    float ***vxyyx_2   = dv_2->xyyx;
    float ***vyzzy_2   = dv_2->yzzy;
    float ***vxzzx_2   = dv_2->xzzx;
    float ***vxxyyzz_2 = dv_2->xxyyzz;
    float ***vyyzz_2   = dv_2->yyzz;
    float ***vxxzz_2   = dv_2->xxzz;
    float ***vxxyy_2   = dv_2->xxyy;

    float ***vxyyx_3   = dv_3->xyyx;
    float ***vyzzy_3   = dv_3->yzzy;
    float ***vxzzx_3   = dv_3->xzzx;
    float ***vxxyyzz_3 = dv_3->xxyyzz;
    float ***vyyzz_3   = dv_3->yyzz;
    float ***vxxzz_3   = dv_3->xxzz;
    float ***vxxyy_3   = dv_3->xxyy;

    float ***vxyyx_4   = dv_4->xyyx;
    float ***vyzzy_4   = dv_4->yzzy;
    float ***vxzzx_4   = dv_4->xzzx;
    float ***vxxyyzz_4 = dv_4->xxyyzz;
    float ***vyyzz_4   = dv_4->yyzz;
    float ***vxxzz_4   = dv_4->xxzz;
    float ***vxxyy_4   = dv_4->xxyy;

    float ***svx = ds_dv->x;
    float ***svy = ds_dv->y;
    float ***svz = ds_dv->z;

    float ***svx_2 = ds_dv_2->x;
    float ***svy_2 = ds_dv_2->y;
    float ***svz_2 = ds_dv_2->z;

    float ***svx_3 = ds_dv_3->x;
    float ***svy_3 = ds_dv_3->y;
    float ***svz_3 = ds_dv_3->z;

    float ***svx_4 = ds_dv_4->x;
    float ***svy_4 = ds_dv_4->y;
    float ***svz_4 = ds_dv_4->z;

    float ***rxx = r->xx;
    float ***ryy = r->yy;
    float ***rzz = r->zz;
    float ***rxy = r->xy;
    float ***ryz = r->yz;
    float ***rxz = r->xz;

    float ***rxx_2 = r_2->xx;
    float ***ryy_2 = r_2->yy;
    float ***rzz_2 = r_2->zz;
    float ***rxy_2 = r_2->xy;
    float ***ryz_2 = r_2->yz;
    float ***rxz_2 = r_2->xz;

    float ***rxx_3 = r_3->xx;
    float ***ryy_3 = r_3->yy;
    float ***rzz_3 = r_3->zz;
    float ***rxy_3 = r_3->xy;
    float ***ryz_3 = r_3->yz;
    float ***rxz_3 = r_3->xz;

    float ***rxx_4 = r_4->xx;
    float ***ryy_4 = r_4->yy;
    float ***rzz_4 = r_4->zz;
    float ***rxy_4 = r_4->xy;
    float ***ryz_4 = r_4->yz;
    float ***rxz_4 = r_4->xz;
    

    /* Set memory variables to zero */
    /* Need to be done seperately due to r is allocated different */
    for (j=1;j<=ny2-FDORDER/2;j++){
        for (i=1;i<=nx2-FDORDER/2;i++){
            for (k=1;k<=nz2-FDORDER/2;k++){
                rxx[j][i][k]=0.0;
                ryy[j][i][k]=0.0;
                rzz[j][i][k]=0.0;
                rxy[j][i][k]=0.0;
                ryz[j][i][k]=0.0;
                rxz[j][i][k]=0.0;
                if (FDORDER_TIME>2) {
                    rxx_2[j][i][k]=0.0;
                    ryy_2[j][i][k]=0.0;
                    rzz_2[j][i][k]=0.0;
                    rxy_2[j][i][k]=0.0;
                    ryz_2[j][i][k]=0.0;
                    rxz_2[j][i][k]=0.0;
                    rxx_3[j][i][k]=0.0;
                    ryy_3[j][i][k]=0.0;
                    rzz_3[j][i][k]=0.0;
                    rxy_3[j][i][k]=0.0;
                    ryz_3[j][i][k]=0.0;
                    rxz_3[j][i][k]=0.0;
                    if (FDORDER_TIME==4) {
                        rxx_4[j][i][k]=0.0;
                        ryy_4[j][i][k]=0.0;
                        rzz_4[j][i][k]=0.0;
                        rxy_4[j][i][k]=0.0;
                        ryz_4[j][i][k]=0.0;
                        rxz_4[j][i][k]=0.0;
                    }
                }

            }
        }
    }
    
    for (j=ny1;j<=ny2;j++){
        for (i=nx1;i<=nx2;i++){
            for (k=nz1;k<=nz2;k++){
                vx[j][i][k]=0.0;
                vy[j][i][k]=0.0;
                vz[j][i][k]=0.0;
                sxx[j][i][k]=0.0;
                syy[j][i][k]=0.0;
                szz[j][i][k]=0.0;
                sxy[j][i][k]=0.0;
                syz[j][i][k]=0.0;
                sxz[j][i][k]=0.0;
                
                if(FDORDER_TIME > 2){
                    vxyyx[j][i][k]=0.0;
                    vyzzy[j][i][k]=0.0;
                    vxzzx[j][i][k]=0.0;
                    vxxyyzz[j][i][k]=0.0;
                    vyyzz[j][i][k]=0.0;
                    vxxzz[j][i][k]=0.0;
                    vxxyy[j][i][k]=0.0;
                    
                    vxyyx_2[j][i][k]=0.0;
                    vyzzy_2[j][i][k]=0.0;
                    vxzzx_2[j][i][k]=0.0;
                    vxxyyzz_2[j][i][k]=0.0;
                    vyyzz_2[j][i][k]=0.0;
                    vxxzz_2[j][i][k]=0.0;
                    vxxyy_2[j][i][k]=0.0;
                    
                    vxyyx_3[j][i][k]=0.0;
                    vyzzy_3[j][i][k]=0.0;
                    vxzzx_3[j][i][k]=0.0;
                    vxxyyzz_3[j][i][k]=0.0;
                    vyyzz_3[j][i][k]=0.0;
                    vxxzz_3[j][i][k]=0.0;
                    vxxyy_3[j][i][k]=0.0;
                    
                    
                    svx[j][i][k]=0.0;
                    svy[j][i][k]=0.0;
                    svz[j][i][k]=0.0;
                    
                    svx_2[j][i][k]=0.0;
                    svy_2[j][i][k]=0.0;
                    svz_2[j][i][k]=0.0;
                    
                    svx_3[j][i][k]=0.0;
                    svy_3[j][i][k]=0.0;
                    svz_3[j][i][k]=0.0;
                    
                    if (FDORDER_TIME==4){
                        vxyyx_4[j][i][k]=0.0;
                        vyzzy_4[j][i][k]=0.0;
                        vxzzx_4[j][i][k]=0.0;
                        vxxyyzz_4[j][i][k]=0.0;
                        vyyzz_4[j][i][k]=0.0;
                        vxxzz_4[j][i][k]=0.0;
                        vxxyy_4[j][i][k]=0.0;
                        
                        svx_4[j][i][k]=0.0;
                        svy_4[j][i][k]=0.0;
                        svz_4[j][i][k]=0.0;
                    }
                }
            }
        }
    }
    
}
