#include "data_structures.h"
#include "fd.h"

void init_velocity(
        Velocity *v,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    v->x = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    v->y = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    v->z = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
}

void free_velocity(
        Velocity *v,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    free_f3tensor(v->x, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(v->y, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(v->z, nrl, nrh, ncl, nch, ndl, ndh);
}

void init_tensor3d(
        Tensor3d *t,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    t->xy = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    t->yz = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    t->xz = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    t->xx = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    t->yy = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    t->zz = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
}

void free_tensor3d(
        Tensor3d *t,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    free_f3tensor(t->xy, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(t->yz, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(t->xz, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(t->xx, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(t->yy, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(t->zz, nrl, nrh, ncl, nch, ndl, ndh);
}

void init_velocity_derivatives_tensor(
        VelocityDerivativesTensor *dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    dv->xyyx = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    dv->yzzy = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    dv->xzzx = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    dv->yyzz = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    dv->xxzz = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    dv->xxyy = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    dv->xxyyzz = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
}

void free_velocity_derivatives_tensor(
        VelocityDerivativesTensor *dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    free_f3tensor(dv->xyyx, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(dv->yzzy, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(dv->xzzx, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(dv->yyzz, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(dv->xxzz, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(dv->xxyy, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(dv->xxyyzz, nrl, nrh, ncl, nch, ndl, ndh);
}

void init_stress_derivatives_wrt_velocity(
        StressDerivativesWrtVelocity *ds_dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {


    ds_dv->x = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    ds_dv->y = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    ds_dv->z = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
}

void free_stress_derivatives_wrt_velocity(
        StressDerivativesWrtVelocity *ds_dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {

    free_f3tensor(ds_dv->x, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(ds_dv->y, nrl, nrh, ncl, nch, ndl, ndh);
    free_f3tensor(ds_dv->z, nrl, nrh, ncl, nch, ndl, ndh);
}
