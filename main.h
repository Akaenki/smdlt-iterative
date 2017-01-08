//
//  main.h
//  StaticComplete
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef main_h
#define main_h

#include "Parameters.h"

/* Initialized forces */
void initForce();

/* Calculate bonding forces */
void forceBond();

/* Calculate exclusive volume interactions */
void forceLJ();

void updateChains();

void applyPBC();

/* Pass the 3N by 3N matrix in bin to *D_ij */
void locateD_ij(bin_t bin, float *D_ij);
//void locateD_ijFine(bin_t bin, float *D_ij);

/* Will cahnge the format of storing the matrix 
 * In the large matrix (matrix->?) the matrix is stored in the order of xx xy xz yx ... zy zz components of each 3*3 matrix from top left to bottom right
 * In D_ij the whole 3N*3N matrix is stored from top left to bottom right */
void locateSelf(float *D_ij);

/* calcualte CoM by chain index */
Vector3D_t CenterOfMass(int i);

/* store the distances of CoM between any chains in a 2D array */
void distCenterOfMass();

void printTrajectory(int t);

double gasdev(long *idum);
long long timer();

/* Debugging Functions */
#if defined(CLUSTER_DEBUG)
void BondLengthCheck(int t);
void printTotalForces();
void initDispl();
void forceLJMod();
#endif

#endif /* main_h */
