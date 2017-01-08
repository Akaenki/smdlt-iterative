//
//  Properties.c
//  StaticComplete
//
//  Created by Linling Miao on 6/13/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* This file contains functions calculating different properties 
 * The print functions are removed. Please create your own print functions */

#include "Properties.h"
#include "main.h"

void storeCoM(){
    if(R_CoM==NULL){
        R_CoM = calloc(NP,sizeof(Vector3D_t));
    }
    
    for(int i = 0; i<NP; ++i){
        R_CoM[i] = CenterOfMass(i);
    }
}

/*void printR_CoM(int t){
    
    if(t%1000==0){
        FILE *Rg = fopen(RgName,"a");
        int i;
        for(i = 0; i<NP; ++i){
            fprintf(Rg,"%d: %lf %lf %lf\n",i,R_CoM[i].x,R_CoM[i].y,R_CoM[i].z);
        }
        fclose(Rg);
    }
}*/

void D_short(int t){

    for(int i = 0; i<NP; ++i){
        
        Vector3D_t dCoM;
        Vector3D_t new_CoM = CenterOfMass(i);
        dCoM.x = new_CoM.x - R_CoM[i].x;
        dCoM.y = new_CoM.y - R_CoM[i].y;
        dCoM.z = new_CoM.z - R_CoM[i].z;
        
        dCoM.x -= round(dCoM.x/L)*L;
        dCoM.y -= round(dCoM.y/L)*L;
        dCoM.z -= round(dCoM.z/L)*L;
        
        MSD_CoM += dCoM.x*dCoM.x + dCoM.y*dCoM.y + dCoM.z*dCoM.z;
        
    }
    
    /*if(t>0 && t%10000==0){
        
        FILE *diffusivity = fopen(DiffName,"a");
        fprintf(diffusivity, "%lf \n",MSD_CoM/(6.0*t*DT)/NP);
        fclose(diffusivity);
        
    }*/
}

void Rendtoend(int t){
    if(t%1000==0){
        //FILE *ree = fopen(ReeName,"a");
        
        for(int i = 0; i<NP; ++i){
            double dx = 0;
            double dy = 0;
            double dz = 0;
            
            for(int j = 1; j<N; ++j){
                Vector3D_t NID = NIDbyIndex(i*N+j,i*N+j-1);
                dx += NID.x;
                dy += NID.y;
                dz += NID.z;
            }
            
            //fprintf(ree,"%lf %lf %lf\n",dx,dy,dz);
        }
    }
}

Vector3D_t NIDbyIndex(int i,int j){
    Vector3D_t NID;
    NID.x = Chain[i].rx - Chain[j].rx;
    NID.y = Chain[i].ry - Chain[j].ry;
    NID.z = Chain[i].rz - Chain[j].rz;
    
    NID.x -= round(NID.x/L)*L;
    NID.y -= round(NID.y/L)*L;
    NID.z -= round(NID.z/L)*L;
    
    return NID;
}

