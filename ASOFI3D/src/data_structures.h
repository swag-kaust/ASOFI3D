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

/* ****************************************************************************
   Allocation and deallocation operations.
*/
void init_velocity(
        Velocity *v,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void free_velocity(
        Velocity *v,
        int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

#endif
