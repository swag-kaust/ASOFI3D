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
