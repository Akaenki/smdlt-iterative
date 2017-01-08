//
//  Properties.h
//  StaticComplete
//
//  Created by Linling Miao on 6/13/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Properties_h
#define Properties_h

#include "Parameters.h"

/* store CoM of all chains */
void storeCoM();

/* Calculate self-diffusivity */
void D_short(int t);
//void printR_CoM(int t);

/* Calculate end to end vectors */
void Rendtoend(int t);

/* Calculate NID by bead indices */
Vector3D_t NIDbyIndex(int i,int j);

//void velocityField();

#endif /* Properties_h */
