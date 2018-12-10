#ifndef __DATA_STRUCTURES__
#define __DATA_STRUCTURES__

// Structure that contains velocity components.
typedef struct {
    float ***x;
    float ***y;
    float ***z;
} Velocity;

typedef struct {
    float ***xy;
    float ***yz;
    float ***xz;
    float ***xx;
    float ***yy;
    float ***zz;
} Tensor3d;

/* Tensor containing derivatives of the velocity.
 * Naming logic is the following: all letters are in pairs,
 * with the first letter in pair denoting velocity component,
 * while the second letter in a pair denoting variable wrt which
 * derivative of the velocity component is taken.
 * All pairs are summed.
 * For example, `xyyx` means `d(v_x)/dy + d(v_y)/dx`.
 */
typedef struct {
    float ***xyyx;
    float ***yzzy;
    float ***xzzx;
    float ***xxyyzz;
    float ***yyzz;
    float ***xxzz;
    float ***xxyy;
} VelocityDerivativesTensor;

/* Velocity derivatives.
 */
typedef struct {
    float xx;
    float xy;
    float xz;
    float yy;
    float yz;
    float zz;
} Strain_ijk;

/* Derivatives of the stress with respect to the velocity components.
 * For example, x-component of this structure is the derivative
 * with respect to the x-component of the velocity.
 */
typedef struct {
    float ***x;
    float ***y;
    float ***z;
} StressDerivativesWrtVelocity;

/*
 * Anisotropic material parameters.
 * Elastic constants and density.
 * Cij are the elements of 6x6 symmetric matrix. They are related to the 
 * components of the fourth-order stiffness tensor Cijkl in the following way:
 * Cijkl -> Cmn,
 * where
 * m = i if i = j,
 * m = 9 - (i + j) if i /= j
 * and the same relation is held between n and kl indices.
 * See, for example, the Appendix of
 *     Cheadle S.P. et al.
 *     Orthorhombic anistropy: A physical seismic modeling study.
 *     Geophysics,  vol. 56, 1991.
 *     https://library.seg.org/doi/pdf/10.1190/1.1442971
 * for further details.
 */
typedef struct {
    float ***C11;
    float ***C22;
    float ***C33;
    float ***C44;
    float ***C55;
    float ***C66;
    float ***C12;
    float ***C13;
    float ***C23;

    // Material density.
    float ***rho;
    // Next three fields are parameters on the half-integer grid
    // (p stands for "+1/2").
    float ***C66ipjp;
    float ***C44jpkp;
    float ***C55ipkp;
} OrthoPar;

/* ****************************************************************************
   Allocation and deallocation operations.
*/
void init_velocity(
        Velocity *v,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void free_velocity(
        Velocity *v,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void init_tensor3d(
        Tensor3d *t,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void free_tensor3d(
        Tensor3d *t,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void init_velocity_derivatives_tensor(
        VelocityDerivativesTensor *dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void free_velocity_derivatives_tensor(
        VelocityDerivativesTensor *dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void init_stress_derivatives_wrt_velocity(
        StressDerivativesWrtVelocity *ds_dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void free_stress_derivatives_wrt_velocity(
        StressDerivativesWrtVelocity *ds_dv,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

#endif
