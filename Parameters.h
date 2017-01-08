//
//  Parameters.h
//  StaticComplete
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cblas.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

//#include <immintrin.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define EP 0.9624
#define THETA0 31.7

#ifndef M_PI
#define M_PI 3.1415926535
#endif


time_t StartTime;
uint16_t iteration;

int HImode;
int AVGmode;
double AVGfactor;

/////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    double x,y,z;
} Vector3D_t;

typedef struct {
    uint64_t x,y,z;
} bin_t;

typedef struct {
    double rx,ry,rz;
    double fx,fy,fz;
} Chain_t;

/*#if defined(LINEARIZE)
typedef struct {
    double xd,yd,zd;//factor align to left
    bin_t C000;
} TriLinear_t;
#endif*/

/* Mobility Matrix */
typedef struct{
    float *DiffMatrix; //matrix calculated in current iteration
    float *selfMatrix; //self matrix calculated in current iteration
    uint32_t *gridCounting; //# of pairs in each grid point is calculated
    
    float *DiffMatrixAvg; //averaged matrix read from the binary file
    float *selfMatrixAvg; //averaged self matrix read from the binary file
/*#if defined(LINEARIZE)
    float *DiffMatrixFine; //Fined matrix (self matrix not included)
#endif*/
    uint32_t *counter; //gridCounting for the pre-averaged matrix
} Matrix_t;

/* Geyer-Winter(TEA) Parameters */
typedef struct{
    double beta_ij;
    float *C_i;
}gwParm_t;
/////////////////////////////////////////////////////////////////////////////
// Gloable
/////////////////////////////////////////////////////////////////////////////
long *idum; //SEED

uint32_t N,NP;
double L; //box size

uint64_t binmax,bmax;

double c_normal; //normalized concentration

double ****m2; //Ewald K space

char* c_normalchar;

Chain_t *Chain;
Matrix_t *matrix;
gwParm_t *gwParm;

double ext;
uint32_t step;

Vector3D_t *R_CoM;
Vector3D_t **distCoM;

double MSD_CoM;

/////////////////////////////////////////////////////////////////////////////
// Others
/////////////////////////////////////////////////////////////////////////////

int numThreads;
int trajStep;

char* outputName;
char* trajName;
char* directory;

bool ReadInitFromFile;
bool CheckOverlap;
bool EwaldSum;
bool MatrixCalculation;
bool isLinearize;

#if defined(CLUSTER_DEBUG)
FILE *debug;
#endif

#endif /* Parameters_h */
