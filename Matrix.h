//
//  Matrix.h
//  StaticComplete
//
//  Created by Linling Miao on 6/13/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h

#include "Parameters.h"

void updateMatrix();

/* Calculate mobility matrices */
void regularRPY(int i,int j,bin_t bin);
void EwaldSumRPY(int i,int j,bin_t bin);

/* will create a .txt file to record values of beta_ij*/
void GeyerWinterCalculation();
void printGeyerWinterParameters();
/* beta_ij (double)*1 8bytes
 * C_i (float)*3N 4bytes each
 */

void printMatrix();
/* Header:
    * Ntotal: total number of beads (uint32)*1
    * N: number of beads per chain (uint32)*1
    * binmax: max number of bins in one dimension (uint64)*1
    * dr: binsize (double)*1
 
 * Body:
    * DiffMatrix: full averaged matrix (float)*9*N*N*binmax**3
    * selfMatrix: self interaction matrix (float)*9*N*N
    * gridCounting: counter (uint32)*binmax**3
 */

/* use to print part of the matrix in .txt file */
void printMatrixPart();

/* get the nearest image distances from real distances */
Vector3D_t getNID(double dx,double dy,double dz);

/* calculate the bin number from the distance of CoM of two chains */
bin_t binCoM(Vector3D_t dCoM_ij);
//bin_t binCoMFine(int i, int j);

#endif /* Matrix_h */
