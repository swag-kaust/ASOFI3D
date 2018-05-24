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

/* Derivatives of the stress with respect to the velocity components.
 * For example, x-component of this structure is the derivative
 * with respect to the x-component of the velocity.
 */
typedef struct {
    float ***x;
    float ***y;
    float ***z;
} StressDerivativesWrtVelocity;

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
