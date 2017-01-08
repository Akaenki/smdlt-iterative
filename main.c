//
//  main.c
//  StaticComplete
//
//  Created by Linling Miao on 6/28/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "main.h"
#include "Initialization.h"
#include "Matrix.h"
#include "Properties.h"

int main(int argc, char * argv[]) {
    Initialization(argc,argv);
    
    int t;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  BEGAIN OF TIME LOOP
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(t = 0; t<TMAX; ++t){
        /*if(t==0){
            triger = 1;
            fprintf(debug,">t = 0 out reference\n");
        }
        if(t>=180447){
            if(t==180447) triger_ct = 0;
            triger_ct++;
            triger=1;
            fprintf(debug,">t = %d track start",t);
        }*/
        //long time1 = timer();
        
        //initDispl();//debug
        
        initForce();
        forceBond();
        forceLJ();
        //forceLJMod();
        
        //if(triger) printTotalForces();
        
        //long time2 = timer() - time1;
        storeCoM();
        distCenterOfMass();
        
        if(t%100==0){
            updateMatrix();
        }
        
        updateChains();

        applyPBC();
        
        //printR_CoM(t);
        //Rendtoend(t);
        //D_short(t);
        
        if(t%trajStep==0) printTrajectory(t);
        
        /*if(triger){
            printTrajectory(t);
            triger = 0;
        }
        
        if(triger_ct==4){
            fclose(debug);
            exit(1);
        }*/

        //BondLengthCheck(t);
        /*if(t%1==0){
            long end = timer() - StartTime;
            long second = end / CLOCKS_PER_SEC;
            long nsecond = end % CLOCKS_PER_SEC;
            printf(">Total: %d, %ldms%ld\nUpdateChains: %ldms%ld\n",t,second,nsecond/1000,time2/CLOCKS_PER_SEC,(time2%CLOCKS_PER_SEC)/1000);
            StartTime = timer();
            exit(1);
        }*/
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  END OF TIME LOOP
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(MatrixCalculation){
        printMatrix();
    }
    printIteration();
   
    if(iteration < MAXITER){
        char* command = malloc(200*sizeof(char));
        sprintf(command,"bash run.sh %s %lf 1",c_normalchar,L);
        system(command);
        
        free(command);
    }
    
    return 0;
}


void initForce(){
    for(int i = 0; i<N*NP; ++i){
        Chain[i].fx = 0.0;
        Chain[i].fy = 0.0;
        Chain[i].fz = 0.0;
    }
}

void forceBond(){
//#pragma omp parallel for private(i) schedule(dynamic)
    for(int i = 0; i<NP; ++i){
        for(int j = 1; j<N; ++j){
            double dx = Chain[i*N+j].rx - Chain[i*N+j-1].rx;
            double dy = Chain[i*N+j].ry - Chain[i*N+j-1].ry;
            double dz = Chain[i*N+j].rz - Chain[i*N+j-1].rz;
            
            Vector3D_t NID = getNID(dx,dy,dz);
            
            double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
            double r = sqrt(rr);
            double Fs = -KAPPA*(r-2.0);
            
            Chain[i*N+j].fx += Fs*NID.x/r;
            Chain[i*N+j].fy += Fs*NID.y/r;
            Chain[i*N+j].fz += Fs*NID.z/r;
            Chain[i*N+j-1].fx -= Fs*NID.x/r;
            Chain[i*N+j-1].fy -= Fs*NID.y/r;
            Chain[i*N+j-1].fz -= Fs*NID.z/r;
            
        }
    }
}

void forceLJ(){
//#pragma omp parallel for private(i) schedule(dynamic)
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<N; ++j){
            int temp;
            if(j<N) temp = 2; else temp = 1;
            for(int k = i*N+j+temp; k<N*NP; ++k){
                double dx = Chain[i*N+j].rx - Chain[k].rx;
                double dy = Chain[i*N+j].ry - Chain[k].ry;
                double dz = Chain[i*N+j].rz - Chain[k].rz;
                
                Vector3D_t NID = getNID(dx,dy,dz);
                
                double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                
                if(rr<25){
                    double ratio = 4.00/rr;
                    double r6 = ratio*ratio*ratio;
                    if(r6>3) r6 = 3;
                    double coeff = (12*EPSILON/rr)*(r6*r6-r6);
                    Chain[i*N+j].fx += coeff*NID.x;
                    Chain[i*N+j].fy += coeff*NID.y;
                    Chain[i*N+j].fz += coeff*NID.z;
                    Chain[k].fx -= coeff*NID.x;
                    Chain[k].fy -= coeff*NID.y;
                    Chain[k].fz -= coeff*NID.z;
                }
            }
        }
    }
}

#if defined(CLUSTER_DEBUG)
void forceLJMod(){
//#pragma omp parallel for private(i) schedule(dynamic)
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<N; ++j){
            int temp;
            if(j<N) temp = 2; else temp = 1;
            for(int k = i*N+j+temp; k<N*NP; ++k){
                double dx = Chain[i*N+j].rx - Chain[k].rx;
                double dy = Chain[i*N+j].ry - Chain[k].ry;
                double dz = Chain[i*N+j].rz - Chain[k].rz;
                
                Vector3D_t NID = getNID(dx,dy,dz);
                
                double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                
                if(rr<4) {
                    double m = 5.0;
                    double coeff = m/sqrt(rr);
                    
                    Chain[i*N+j].fx += coeff*NID.x;
                    Chain[i*N+j].fy += coeff*NID.y;
                    Chain[i*N+j].fz += coeff*NID.z;
                    Chain[k].fx -= coeff*NID.x;
                    Chain[k].fy -= coeff*NID.y;
                    Chain[k].fz -= coeff*NID.z;
                }
                else if(rr<25){
                    double ratio = 4.00/rr;
                    double r6 = ratio*ratio*ratio;
                    if(r6>3) r6 = 3;
                    double coeff = (12*EPSILON/rr)*(r6*r6-r6);
                    Chain[i*N+j].fx += coeff*NID.x;
                    Chain[i*N+j].fy += coeff*NID.y;
                    Chain[i*N+j].fz += coeff*NID.z;
                    Chain[k].fx -= coeff*NID.x;
                    Chain[k].fy -= coeff*NID.y;
                    Chain[k].fz -= coeff*NID.z;
                }
            }
        }
    }
}
#endif

void updateChains(){
    int nn = 3*N;
    double p = sqrt(2*DT);
    
    float *RR = calloc(3*N*NP,sizeof(float));
    for(int i = 0; i<3*N*NP; ++i){
        RR[i] = gasdev(idum);
    }
    
    if(HImode){ //Preaveraged HI
#pragma omp parallel for schedule(static) //change schedule to dynamic may increase parallel overhead
        for(int i = 0; i<NP; ++i){
            float *D_ij = calloc(nn*nn,sizeof(float));
            float *D_force = calloc(nn,sizeof(float));
            float *D_noise = calloc(nn,sizeof(float));
            
            float *force = calloc(nn,sizeof(float));
            float *random = calloc(nn,sizeof(float));
            
            for(int j = 0; j<NP; ++j){
                
                if(i == j){
                    locateSelf(D_ij);
                } else{
                    //bin_t bin = binCoMFine(i,j);
                    //locateD_ijFine(bin, D_ij);
                    
                    //printf("i: %d, j: %d, rc: %lf %lf %lf\n",i,j,distCoM[i][j].dx,distCoM[i][j].dy,distCoM[i][j].dz);
                    
                    if(0) {
                        //LinearizeMatrix(i, j, D_ij);
                    } else{
                        bin_t bin = binCoM(distCoM[i][j]);
                        locateD_ij(bin, D_ij);
                    }
                    
                    //exit(EXIT_SUCCESS);
                    
                }
                
                for(int k = 0; k<N; ++k){
                    force[3*k] = Chain[j*N+k].fx;
                    force[3*k+1] = Chain[j*N+k].fy;
                    force[3*k+2] = Chain[j*N+k].fz;
                }
                
                for(int k = 0; k<nn; ++k){
                    random[k] = RR[j*N*3+k];
                }
                
                /*  Matrix Multiplication Function grabed from CBLAS library, optimized and multithreaded by openBLAS
                 *  cblas_sgemm for single precision (float), cblas_dgemm for double precision (double)
                 *  Matrix Multiplication Follows the form: C = alpha * ( A dot B ) + beta * C
                 *  ( *1 , *2 , *3 , M , N , K , alpha , A[M*K] , Ida , B[K*N] , Idb , beta , C[M*N] , Idc)  
                 *  Detail see documents for LINPACK and CBLAS */
                cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,D_ij,nn,force,1,1,D_force,1);
                
                if(i == j){  //same as setting beta_ii = 1 for diagonal components
                    for(int ii = 0; ii<nn; ++ii){
                        D_ij[ii*nn+ii] = 0; //D_ij[ii*nn+ii] /= gwParm->beta_ij;
                    }
                }
                
                cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,gwParm->beta_ij,D_ij,nn,random,1,1,D_noise,1);
                
                
                /*for(int ii = 0; ii<nn; ++ii){
                    for(int jj = 0; jj<nn; ++jj){
                        D_force[ii] += force[jj]*D_ij[ii*nn+jj];
                        //if(i==j && ii==jj){
                            //nothing
                        //} else{
                        //    D_noise[ii] += gwParm->beta_ij*RR[j*3*N+jj]*D_ij[ii*nn+jj];
                        //}
                    }
                }*/
            }
            
            free(D_ij);
            free(force);
            free(random);
            
            /* update chains */
            for(int k = 0; k<N; ++k){
                Chain[i*N+k].rx += DT*D_force[k*3] + p*gwParm->C_i[k*3]*(D_noise[k*3]+RR[3*(i*N+k)]);
                Chain[i*N+k].ry += DT*D_force[k*3+1] + p*gwParm->C_i[k*3+1]*(D_noise[k*3+1]+RR[3*(i*N+k)+1]);
                Chain[i*N+k].rz += DT*D_force[k*3+2] + p*gwParm->C_i[k*3+2]*(D_noise[k*3+2]+RR[3*(i*N+k)+2]);
            }
            
            
            free(D_force);
            free(D_noise);
        }
        
    } else{ //Free Draining
#pragma omp parallel for //private(i) schedule(dynamic) //schedule(static,2)
        for(int i = 0; i<NP*N; ++i){
            Chain[i].rx += DT*Chain[i].fx+p*RR[3*i];
            Chain[i].ry += DT*Chain[i].fy+p*RR[3*i+1];
            Chain[i].rz += DT*Chain[i].fz+p*RR[3*i+2];
        }
    }
    
    free(RR);
    
/*#if defined(CLUSTER_DEBUG)
     for(int i = 0; i<N*NP; i++){
        Chain[i].rx += D_total[3*i];
        Chain[i].rx += D_total[3*i+1];
        Chain[i].rx += D_total[3*i+2];
    }
#endif*/
}

void locateD_ij(bin_t bin, float *D_ij){
    int nn = 3*N;
    uint64_t start = bin.x*binmax*binmax*nn*nn+bin.y*binmax*nn*nn+bin.z*nn*nn;
    
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            int start2 = i*9*N+j*9;
            
            for(int k = 0; k<9; ++k){
                D_ij[(i*3+k/3)*nn+j*3+k%3] = matrix->DiffMatrixAvg[start+start2+k];
            }
        }
    }
}

/*#if defined(LINEARIZE)
void locateD_ijFine(bin_t bin, float *D_ij){
    int nn = 3*N;
    uint64_t start = bin.x*bmax*bmax*nn*nn+bin.y*bmax*nn*nn+bin.z*nn*nn;
    
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            int start2 = i*9*N+j*9;
            
            for(int k = 0; k<9; ++k){
                D_ij[(i*3+k/3)*nn+j*3+k%3] = matrix->DiffMatrixFine[start+start2+k];
            }
        }
    }
}
#endif*/

void locateSelf(float *D_ij){
    int nn = 3*N;
    
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            int start = i*9*N+j*9;
            
            for(int k = 0; k<9; ++k){
                D_ij[(i*3+k/3)*nn+j*3+k%3] = matrix->selfMatrixAvg[start+k];
            }
        }
    }
}

void applyPBC(){
    for(int i = 0; i<N*NP; ++i){
        Chain[i].rx -= L*round(Chain[i].rx/L);
        Chain[i].ry -= L*round(Chain[i].ry/L);
        Chain[i].rz -= L*round(Chain[i].rz/L);
    }
}

Vector3D_t CenterOfMass(int i){
    Vector3D_t CoM;
    CoM.x = Chain[i*N].rx;
    CoM.y = Chain[i*N].ry;
    CoM.z = Chain[i*N].rz;
    for(int j = 1; j<N; ++j){
        double dx = Chain[i*N+j].rx - Chain[i*N+j-1].rx;
        double dy = Chain[i*N+j].ry - Chain[i*N+j-1].ry;
        double dz = Chain[i*N+j].rz - Chain[i*N+j-1].rz;
        
        Vector3D_t NID = getNID(dx,dy,dz);
        
        CoM.x += (N-j)*NID.x/N;
        CoM.y += (N-j)*NID.y/N;
        CoM.z += (N-j)*NID.z/N;
    }
    
    CoM.x -= round(CoM.x/L)*L;
    CoM.y -= round(CoM.y/L)*L;
    CoM.z -= round(CoM.z/L)*L;
    
    return CoM;
}

void distCenterOfMass(){
    if(distCoM==NULL){
        distCoM = calloc(NP,sizeof(distCoM[0]));
        int i;
        for(i = 0; i<NP; ++i){
            distCoM[i] = calloc(NP,sizeof(distCoM[0][0]));
        }
    }
    
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<NP; ++j){
            Vector3D_t CoM_i = CenterOfMass(i);
            Vector3D_t CoM_j = CenterOfMass(j);
            
            double dx = CoM_j.x - CoM_i.x;
            double dy = CoM_j.y - CoM_i.y;
            double dz = CoM_j.z - CoM_i.z;

            dx -= round(dx/L)*L;
            dy -= round(dy/L)*L;
            dz -= round(dz/L)*L;
            
            //printf("%d-%d: %lf %lf %lf\n",i,j,dx,dy,dz);
            
            /*bin_t bin;
            if(dx>=0) bin.x = (int)(dx/dr)+binmax/2; else bin.x = floor(dx/dr)+binmax/2;
            if(dy>=0) bin.y = (int)(dy/dr)+binmax/2; else bin.y = floor(dy/dr)+binmax/2;
            if(dz>=0) bin.z = (int)(dz/dr)+binmax/2; else bin.z = floor(dz/dr)+binmax/2;*/
            
            distCoM[i][j].x = dx;
            distCoM[i][j].y = dy;
            distCoM[i][j].z = dz;
            
            /*FILE *check = fopen("./output/check.txt","a");
            //fprintf(check,"%d:%lf %lf %lf %d:%lf %lf %lf\n",i,CoM_i.x,CoM_i.y,CoM_i.z,j,CoM_j.x,CoM_j.y,CoM_j.z);
            fprintf(check,"%d-%d: %llu %llu %llu\n",i,j,bin.x,bin.y,bin.z);
            fclose(check);*/
        }
    }
}

void printTrajectory(int t){
    FILE *Trajectory;
    Trajectory = fopen(trajName, "a");
    fprintf(Trajectory, "%d\n%d\n", N*NP, t);
    for(int i = 0; i<N*NP; ++i){
        fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
    }
    fclose(Trajectory);
}

double gasdev(long *idum){
    double ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0){
        do{
            v1 = 2.0*ran1(idum)-1.0;
            v2 = 2.0*ran1(idum)-1.0;
            rsq = v1*v1+v2*v2;
        }
        while(rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else{
        iset = 0;
        return gset;
    }
} 

long long timer(){
    struct timespec ts;
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return 1000000000*ts.tv_sec + ts.tv_nsec;
}

/*#if defined(CLUSTER_DEBUG)
void printTotalForces(){
    fprintf(debug,">Total Forces: \n");
    
    for(int i = 0; i<N*NP; i++){
        fprintf(debug,"Index: %d, fx: %.4lf, fy: %.4lf, fz: %.4lf\n",i,Chain[i].fx,Chain[i].fy,Chain[i].fz);
    }
}

void BondLengthCheck(int t){
    
    if(test && triger==0 && triger_ct<100){
        printf("t: %d, # offs: %d\n", t, test);
        printf("SEED: %ld\n",*idum);
        //triger = 1;
        triger_ct++;
    }
    else if(test){
        printf("t: %d, # offs: %d\n", t, test);
        triger_ct++;
    }
}
#endif*/