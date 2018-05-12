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

#endif
