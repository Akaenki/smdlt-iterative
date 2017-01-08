//
//  Matrix.c
//  StaticComplete
//
//  Created by Linling Miao on 6/13/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "Matrix.h"

void updateMatrix(){
    if(MatrixCalculation || HImode==0){
        for(int i = 0; i<NP; ++i){
            for(int j = 0; j<NP; ++j){
                bin_t bin = binCoM(distCoM[i][j]);
                //printf("%llu,%llu,%llu\n",bin.x,bin.y,bin.z);
                if(EwaldSum) EwaldSumRPY(i,j,bin);
                else regularRPY(i,j,bin);
            }
        }
        
        if(step%100==0){
            printf("-----------step: #%d-------------\n",step);
            //printMM(1, 1);
            //if(step>0) printMatrix();
            //printMatrixPart();
        }
    } else{
        if(step%100==0){
            printf("-----------step: #%d-------------\n",step);
        }
    }
    
    step++;
}

void regularRPY(int i,int j,bin_t bin){
    int nn = 3*N;
    if(i==j){
/* Parallel here is not recomanded */
//#pragma omp parallel for private(ii) schedule(dynamic)
        for(int ii = 0; ii<N; ++ii){
            for(int jj = 0; jj<N; ++jj){
                int start0 = ii*9*N+jj*9;
                if(ii==jj){
                    matrix->selfMatrix[start0] += 1.0;
                    matrix->selfMatrix[start0+4] += 1.0;
                    matrix->selfMatrix[start0+8] += 1.0;
                } else{
                    double dx = Chain[i*N+ii].rx - Chain[j*N+jj].rx;
                    double dy = Chain[i*N+ii].ry - Chain[j*N+jj].ry;
                    double dz = Chain[i*N+ii].rz - Chain[j*N+jj].rz;
                    
                    Vector3D_t NID = getNID(dx,dy,dz);
                    //printf("self %lf,%lf,%lf\n",NID.x,NID.y,NID.z);
                    double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                    double r = sqrt(rr);
                    
                    double mm1,mm2;
                    if(r>=2.0){
                        mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                        mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                    } else{
                        mm1 = 1.0-9.0*r/32.0;
                        mm2 = 3.0/(32.0*r);
                    }
                    
                    matrix->selfMatrix[start0] += mm1+mm2*NID.x*NID.x;
                    matrix->selfMatrix[start0+1] += mm2*NID.x*NID.y;
                    matrix->selfMatrix[start0+2] += mm2*NID.x*NID.z;
                    matrix->selfMatrix[start0+3] += mm2*NID.y*NID.x;
                    matrix->selfMatrix[start0+4] += mm1+mm2*NID.y*NID.y;
                    matrix->selfMatrix[start0+5] += mm2*NID.y*NID.z;
                    matrix->selfMatrix[start0+6] += mm2*NID.z*NID.x;
                    matrix->selfMatrix[start0+7] += mm2*NID.z*NID.y;
                    matrix->selfMatrix[start0+8] += mm1+mm2*NID.z*NID.z;
                }
            }
        }
    } else{
        for(int ii = 0; ii<N; ++ii){
            for(int jj = 0; jj<N; ++jj){
                uint64_t start = bin.x*binmax*binmax*nn*nn+bin.y*binmax*nn*nn+bin.z*nn*nn+ii*9*N+jj*9;
                
                double dx = Chain[i*N+ii].rx - Chain[j*N+jj].rx;
                double dy = Chain[i*N+ii].ry - Chain[j*N+jj].ry;
                double dz = Chain[i*N+ii].rz - Chain[j*N+jj].rz;
                
                Vector3D_t NID = getNID(dx,dy,dz);
                //printf("dij %lf,%lf,%lf\n",NID.x,NID.y,NID.z);
                double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                double r = sqrt(rr);
                
                double mm1,mm2;
                if(r>=2.0){
                    mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                    mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                } else{
                    mm1 = 1.0-9.0*r/32.0;
                    mm2 = 3.0/(32.0*r);
                }
                
                matrix->DiffMatrix[start] += mm1+mm2*NID.x*NID.x;
                matrix->DiffMatrix[start+1] += mm2*NID.x*NID.y;
                matrix->DiffMatrix[start+2] += mm2*NID.x*NID.z;
                matrix->DiffMatrix[start+3] += mm2*NID.y*NID.x;
                matrix->DiffMatrix[start+4] += mm1+mm2*NID.y*NID.y;
                matrix->DiffMatrix[start+5] += mm2*NID.y*NID.z;
                matrix->DiffMatrix[start+6] += mm2*NID.z*NID.x;
                matrix->DiffMatrix[start+7] += mm2*NID.z*NID.y;
                matrix->DiffMatrix[start+8] += mm1+mm2*NID.z*NID.z;
            }
        }
        matrix->gridCounting[bin.x*binmax*binmax+bin.y*binmax+bin.z] += 1;
    }
}

void EwaldSumRPY(int i,int j,bin_t bin){
    int nn = 3*N;
    double vol = L*L*L;
    double alpha = 6.4/L;
    //double rcutoff = L/2.0;
    
    //printf("alpha: %lf\n",alpha);
    
    if(i==j){
/* Parallel nor recommanded here */
//#pragma omp parallel for private(ii) schedule(dynamic)
        for(int ii = 0; ii<N; ++ii){
            for(int jj = 0; jj<N; ++jj){
                int start0 = ii*9*N+jj*9;
                
                double dx = Chain[i*N+ii].rx - Chain[j*N+jj].rx;
                double dy = Chain[i*N+ii].ry - Chain[j*N+jj].ry;
                double dz = Chain[i*N+ii].rz - Chain[j*N+jj].rz;
                
                //printf("dx: %lf, dy: %lf, dz: %lf\n",dx,dy,dz);
                
                //Reciprocal Part
                for(int kx = -KMAX; kx<KMAX+1; ++kx){
                    double rkx = 2.0*M_PI*kx/L;
                    for(int ky = -KMAX; ky<KMAX+1; ++ky){
                        double rky = 2.0*M_PI*ky/L;
                        for(int kz = -KMAX; kz<KMAX+1; ++kz){
                            double rkz = 2.0*M_PI*kz/L;
                            double kk = kx*kx+ky*ky+kz*kz;
                            if(kk != 0){
                                double mm4 = cos(rkx*dx+rky*dy+rkz*dz);
                                
                                //printf("cos(kr) = %e\n",mm4);
                                
                                matrix->selfMatrix[start0] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][0]/vol;
                                matrix->selfMatrix[start0+1] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][1]/vol;
                                matrix->selfMatrix[start0+2] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][2]/vol;
                                matrix->selfMatrix[start0+3] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][3]/vol;
                                matrix->selfMatrix[start0+4] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][4]/vol;
                                matrix->selfMatrix[start0+5] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][5]/vol;
                                matrix->selfMatrix[start0+6] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][6]/vol;
                                matrix->selfMatrix[start0+7] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][7]/vol;
                                matrix->selfMatrix[start0+8] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][8]/vol;
                                
                                //printf("xx: %f, yy: %f, zz:%f\n",matrix->selfMatrix[start0],matrix->selfMatrix[start0+4],matrix->selfMatrix[start0+8]);
                            }
                        }
                    }
                }
                
                //Real Part (Cutoff 0.5L)
                dx -= round(dx/L)*L; dy -= round(dy/L)*L; dz -= round(dz/L)*L;
                double rr = dx*dx+dy*dy+dz*dz;
                double r = sqrt(rr);
                
                if(jj == ii){
                    double diag = 1.0-6.0/sqrt(M_PI)*alpha+40.0/3.0/sqrt(M_PI)*alpha*alpha*alpha;
                    
                    matrix->selfMatrix[start0] += diag;
                    matrix->selfMatrix[start0+4] += diag;
                    matrix->selfMatrix[start0+8] += diag;
                    
                    //printf("diagonal: %lf\n\n",diag);
                }
                else{
                    double c1 = 0.75/r+0.5/pow(r,3.0);
                    double c2 = 4.0*pow(alpha,7.0)*pow(r,4.0)+3.0*pow(alpha,3.0)*rr-20.0*pow(alpha,5.0)*rr-4.5*alpha+14.0*pow(alpha,3.0)+alpha/rr;
                    double c3 = 0.75/r-1.5/pow(r,3.0);
                    double c4 = -4.0*pow(alpha,7.0)*pow(r,4.0)-3.0*pow(alpha,3.0)*rr+16.0*pow(alpha,5.0)*rr+1.5*alpha-2.0*pow(alpha,3.0)-3.0*alpha/rr;
                    
                    double mm1,mm2;
                    if(r>=2.0){
                        mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI);
                        mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI);
                    }
                    else{
                        mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI)+(1.0-9.0*r/32.0-0.75/r-0.5/pow(r,3.0));
                        mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI)+(3.0*r/32.0-0.75/r+1.5/pow(r,3.0));
                    }
                    
                    matrix->selfMatrix[start0] += mm1+mm2*dx*dx/rr;
                    matrix->selfMatrix[start0+1] += mm2*dx*dy/rr;
                    matrix->selfMatrix[start0+2] += mm2*dx*dz/rr;
                    matrix->selfMatrix[start0+3] += mm2*dy*dx/rr;
                    matrix->selfMatrix[start0+4] += mm1+mm2*dy*dy/rr;
                    matrix->selfMatrix[start0+5] += mm2*dy*dz/rr;
                    matrix->selfMatrix[start0+6] += mm2*dz*dx/rr;
                    matrix->selfMatrix[start0+7] += mm2*dz*dy/rr;
                    matrix->selfMatrix[start0+8] += mm1+mm2*dz*dz/rr;
                 
                    //printf("xx: %f, yy: %f, zz:%f\n\n",matrix->selfMatrix[start0],matrix->selfMatrix[start0+4],matrix->selfMatrix[start0+8]);
                }
            }
        }
    } else{
//#pragma omp parallel for private(ii) schedule(dynamic)
        for(int ii = 0; ii<N; ++ii){
            for(int jj = 0; jj<N; ++jj){
                uint64_t start = bin.x*binmax*binmax*nn*nn+bin.y*binmax*nn*nn+bin.z*nn*nn+ii*9*N+jj*9;
                
                double dx = Chain[i*N+ii].rx - Chain[j*N+jj].rx;
                double dy = Chain[i*N+ii].ry - Chain[j*N+jj].ry;
                double dz = Chain[i*N+ii].rz - Chain[j*N+jj].rz;
                
                //printf("dx: %lf, dy: %lf, dz: %lf\n",dx,dy,dz);
                //printf("dx:
                
                //Reciprocal Part
                for(int kx = -KMAX; kx<KMAX+1; ++kx){
                    double rkx = 2.0*M_PI*kx/L;
                    for(int ky = -KMAX; ky<KMAX+1; ++ky){
                        double rky = 2.0*M_PI*ky/L;
                        for(int kz = -KMAX; kz<KMAX+1; ++kz){
                            double rkz = 2.0*M_PI*kz/L;
                            double kk = kx*kx+ky*ky+kz*kz;
                            if(kk != 0 && kk <= KMAX*KMAX){
                                double mm4 = cos(rkx*dx+rky*dy+rkz*dz);
                                //printf("%d %d %d:\n",kx,ky,kz);
                                //printf("cos(kr) = %e\n",mm4);
                                
                                matrix->DiffMatrix[start] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][0]/vol;
                                matrix->DiffMatrix[start+1] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][1]/vol;
                                matrix->DiffMatrix[start+2] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][2]/vol;
                                matrix->DiffMatrix[start+3] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][3]/vol;
                                matrix->DiffMatrix[start+4] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][4]/vol;
                                matrix->DiffMatrix[start+5] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][5]/vol;
                                matrix->DiffMatrix[start+6] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][6]/vol;
                                matrix->DiffMatrix[start+7] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][7]/vol;
                                matrix->DiffMatrix[start+8] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][8]/vol;
                                
                                //printf("xx: %f, yy: %f, zz:%f\n",matrix->DiffMatrix[start],matrix->DiffMatrix[start+4],matrix->DiffMatrix[start+8]);
                            }
                        }
                    }
                }
                
                //Real Part (Cutoff 0.5L)
                /*for(int nx = -nmax; nx<nmax+1; ++nx){
                    for(int ny = -nmax; ny<nmax+1; ++ny){
                        for(int nz = -nmax; nz<nmax+1; ++nz){
                            
                            double ddx = dx+nx*L;
                            double ddy = dy+ny*L;
                            double ddz = dz+nz*L;
                            
                            double rr = ddx*ddx+ddy*ddy+ddz*ddz;
                            double r = sqrt(rr);
                            
                            if(r<rcutoff){
                                double c1 = 0.75/r+0.5/pow(r,3.0);
                                double c2 = 4.0*pow(alpha,7.0)*pow(r,4.0)+3.0*pow(alpha,3.0)*rr-20.0*pow(alpha,5.0)*rr-4.5*alpha+14.0*pow(alpha,3.0)+alpha/rr;
                                double c3 = 0.75/r-1.5/pow(r,3.0);
                                double c4 = -4.0*pow(alpha,7.0)*pow(r,4.0)-3.0*pow(alpha,3.0)*rr+16.0*pow(alpha,5.0)*rr+1.5*alpha-2.0*pow(alpha,3.0)-3.0*alpha/rr;
                                
                                printf("\nc1: %lf, c2: %lf, c3: %lf, c4: %lf\n\n",c1,c2,c3,c4);
                                
                                double mm1,mm2;
                                if(r>=2.0){
                                    mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI);
                                    mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI);
                                }
                                else{
                                    mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI)+(1.0-9.0*r/32.0-3.0/4.0/r-0.5*pow(r,-3.0));
                                    mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI)+(3.0*r/32.0-3.0/4.0/r+1.5*pow(r,-3.0));
                                }
                                
                                //printf("mm1: %lf, mm2: %lf\n",mm1,mm2);
                                
                                matrix->DiffMatrix[start] += mm1+mm2*dx*dx/rr;
                                matrix->DiffMatrix[start+1] += mm2*dx*dy/rr;
                                matrix->DiffMatrix[start+2] += mm2*dx*dz/rr;
                                matrix->DiffMatrix[start+3] += mm2*dy*dx/rr;
                                matrix->DiffMatrix[start+4] += mm1+mm2*dy*dy/rr;
                                matrix->DiffMatrix[start+5] += mm2*dy*dz/rr;
                                matrix->DiffMatrix[start+6] += mm2*dz*dx/rr;
                                matrix->DiffMatrix[start+7] += mm2*dz*dy/rr;
                                matrix->DiffMatrix[start+8] += mm1+mm2*dz*dz/rr;
                                
                                printf("%lf, %lf, %lf\n",mm1+mm2*dx*dx/rr,mm1+mm2*dy*dy/rr,mm1+mm2*dz*dz/rr);
                                
                                printf("xx: %f, yy: %f, zz:%f\n\n",matrix->DiffMatrix[start],matrix->DiffMatrix[start+4],matrix->DiffMatrix[start+8]);
                                
                            }
                        }
                    }
                }*/
                
                dx -= round(dx/L)*L; dy -= round(dy/L)*L; dz -= round(dz/L)*L;
                double rr = dx*dx+dy*dy+dz*dz;
                double r = sqrt(rr);
                
                //printf("r: %lf \n",r);
                
                double c1 = 0.75/r+0.5/pow(r,3.0);
                double c2 = 4.0*pow(alpha,7.0)*pow(r,4.0)+3.0*pow(alpha,3.0)*rr-20.0*pow(alpha,5.0)*rr-4.5*alpha+14.0*pow(alpha,3.0)+alpha/rr;
                double c3 = 0.75/r-1.5/pow(r,3.0);
                double c4 = -4.0*pow(alpha,7.0)*pow(r,4.0)-3.0*pow(alpha,3.0)*rr+16.0*pow(alpha,5.0)*rr+1.5*alpha-2.0*pow(alpha,3.0)-3.0*alpha/rr;
                
                //printf("c1: %lf, c2: %lf, c3: %lf, c4: %lf\n",c1,c2,c3,c4);
                
                double mm1,mm2;
                if(r>=2.0){
                    mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI);
                    mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI);
                    
                    //printf("OUT mm1: %lf, mm2: %lf\n",mm1,mm2);
                }
                else{
                    mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI)+(1.0-9.0*r/32.0-3.0/4.0/r-0.5*pow(r,-3.0));
                    mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI)+(3.0*r/32.0-3.0/4.0/r+1.5*pow(r,-3.0));
                    
                    //printf("IN mm1: %lf, mm2: %lf\n",mm1,mm2);
                }
                
                //printf("mm1: %lf, mm2: %lf\n",mm1,mm2);
                
                matrix->DiffMatrix[start] += mm1+mm2*dx*dx/rr;
                matrix->DiffMatrix[start+1] += mm2*dx*dy/rr;
                matrix->DiffMatrix[start+2] += mm2*dx*dz/rr;
                matrix->DiffMatrix[start+3] += mm2*dy*dx/rr;
                matrix->DiffMatrix[start+4] += mm1+mm2*dy*dy/rr;
                matrix->DiffMatrix[start+5] += mm2*dy*dz/rr;
                matrix->DiffMatrix[start+6] += mm2*dz*dx/rr;
                matrix->DiffMatrix[start+7] += mm2*dz*dy/rr;
                matrix->DiffMatrix[start+8] += mm1+mm2*dz*dz/rr;
                
                //printf("%lf, %lf, %lf\n",mm1+mm2*dx*dx/rr,mm1+mm2*dy*dy/rr,mm1+mm2*dz*dz/rr);
            
                //printf("xx: %f, yy: %f, zz:%f\n\n",matrix->DiffMatrix[start],matrix->DiffMatrix[start+4],matrix->DiffMatrix[start+8]);
            }
        }
        matrix->gridCounting[bin.x*binmax*binmax+bin.y*binmax+bin.z] += 1;
    }
}


void GeyerWinterCalculation(){
    uint64_t bintotal = binmax*binmax*binmax;
    int nn = 3*N;
    double eps = 0.0;
    
    for(int i = 0; i<bintotal; ++i){
        for(int j = 0; j<nn*nn; ++j){
            /*if(counter[i]){
                eps += matrix->DiffMatrix[i*nn*nn+j]/gridCounting[i];
                
            }*/
            eps += matrix->DiffMatrixAvg[i*nn*nn+j];
        }
    }
    
    for(int i = 0; i<nn; ++i){
        for(int j = 0; j<nn; ++j){
            //eps += matrix->selfMatrix[i]/step/NP;
            if(i!=j) eps += matrix->selfMatrixAvg[i*nn+j];
        }
    }
    
    //eps -= nn; /* Subtract diagonal terms in matrix->selfMatrix */
    
    uint64_t nt = bintotal*nn; /* Number of total rows of the entire matrix */
    
    printf("-----Geyer-Winter Parameters-----\n");
    
    eps /= nt*nn - nn;
    printf(">epsilon=%e\n", eps);
    
    gwParm->beta_ij = (1-sqrt(1-((nt-1)*eps*eps-(nt-2)*eps)))/((nt-1)*eps*eps-(nt-2)*eps);
    //printf(">>>>row=%llu\n",nt);
    printf(">beta=%f\n", gwParm->beta_ij);
    
    printf("---------------------------------\n");
    
    gwParm->C_i = calloc(nn,sizeof(float));
    
    for(int i = 0; i<bintotal; ++i){
        for(int j = 0; j<N; ++j){
            for(int k = 0; k<N; ++k){
                for(int l= 0; l<9; ++l){
                    if(matrix->gridCounting[i]){
                        gwParm->C_i[j*3+l/3] += (matrix->DiffMatrix[i*nn*nn+j*N*9+k*9+l]/matrix->gridCounting[i])*(matrix->DiffMatrix[i*nn*nn+j*N*9+k*9+l]/matrix->gridCounting[i]);
                    }
                    //gwParm->C_i[j*3+l/3] += matrix->DiffMatrixAvg[i*nn*nn+j*N*9+k*9+l]*matrix->DiffMatrixAvg[i*nn*nn+j*N*9+k*9+l];
                }
            }
        }
    }
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            for(int k = 0; k<9; ++k){
                gwParm->C_i[i*3+k/3] += matrix->selfMatrixAvg[i*9*N+j*9+k]*matrix->selfMatrixAvg[i*9*N+j*9+k];
            }
        }
    }
    for(int i = 0; i<nn; ++i){
        gwParm->C_i[i] -= 1.0;
        gwParm->C_i[i] = gwParm->C_i[i]*gwParm->beta_ij*gwParm->beta_ij+1.0;
        gwParm->C_i[i] = 1.0/sqrt(gwParm->C_i[i]);
        //printf("C_%d=%f\n",i,gwParm->C_i[i]);
    }
    
    char* str = malloc(30*sizeof(char));
    sprintf(str,"%s/log_%s.txt",directory,c_normalchar);
    
    FILE *log = fopen(str,"a");
    fprintf(log,"iteration #%d: beta = %f\n",iteration,gwParm->beta_ij);
    
    fclose(log);
    
    free(str);
}

void printGeyerWinterParameters(){
    /*FILE *log = fopen("./output/log.txt","a"); //MARK: log name
    fprintf(log,"iteration #%d: beta = %f\n",iteration,beta_ij);
    fclose(log);*/
    
    FILE *gwParameters;
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/GeyerWinterParam_%d.bin",directory,iteration);
    
    gwParameters = fopen(str,"wb");
    fwrite(&gwParm->beta_ij,sizeof(double),1,gwParameters);
    fwrite(gwParm->C_i,sizeof(float),3*N,gwParameters);
    fclose(gwParameters);
    free(str);
}

void printMatrix(){
    uint32_t stotal = 9*N*N;
    uint64_t btotal = binmax*binmax*binmax;
    
    FILE *Diffusion;
    char* str = malloc(100*sizeof(char));
    sprintf(str, "%s/Matrix_%s_%u.bin",directory,c_normalchar,iteration);//MARK: name
    
    Diffusion = fopen(str,"wb");
    
    printf(">Generating Binary Averaged Matrix File...\n");
    
    uint32_t Ntotal = N*NP;
    double dr = BIN_SIZE;
    
    fwrite(&Ntotal,sizeof(uint32_t),1,Diffusion);
    fwrite(&N,sizeof(uint32_t),1,Diffusion);
    fwrite(&binmax,sizeof(uint64_t),1,Diffusion);
    fwrite(&dr,sizeof(double),1,Diffusion);
    
    if(AVGmode==0 || HImode==0){
        for(int i = 0; i<btotal; ++i){
            if(matrix->gridCounting[i]){
                for(int j = 0; j<9*N*N; ++j){
                    float value = matrix->DiffMatrix[i*9*N*N+j]/matrix->gridCounting[i];
                    fwrite(&value,sizeof(float),1,Diffusion);
                }
            } else{
                for(int j = 0; j<9*N*N; ++j){
                    float value = matrix->DiffMatrix[i*9*N*N+j];
                    fwrite(&value,sizeof(float),1,Diffusion);
                }
            }
        }
        for(int i = 0; i<stotal; ++i){
            float value = matrix->selfMatrix[i]/step/NP;
            fwrite(&value,sizeof(float),1,Diffusion);
        }
    } else if(AVGmode==1){
        for(int i = 0; i<btotal; ++i){
            if(matrix->gridCounting[i]){
                for(int j = 0; j<9*N*N; ++j){
                    float value = AVGfactor*(matrix->DiffMatrix[i*9*N*N+j]/matrix->gridCounting[i]) + (1.0 - AVGfactor)*matrix->DiffMatrixAvg[i];
                    fwrite(&value,sizeof(float),1,Diffusion);
                }
            } else{
                for(int j = 0; j<9*N*N; ++j){
                    float value = (1.0 - AVGfactor)*matrix->DiffMatrixAvg[i*9*N*N+j];
                    fwrite(&value,sizeof(float),1,Diffusion);
                }
            }
        }
        for(int i = 0; i<stotal; ++i){
            float value = AVGfactor*(matrix->selfMatrix[i]/step/NP) + (1.0 - AVGfactor)*matrix->selfMatrixAvg[i];
            fwrite(&value,sizeof(float),1,Diffusion);
        }
    }
    
    
    fwrite(matrix->gridCounting,sizeof(uint32_t),btotal,Diffusion);
    fwrite(&step,sizeof(uint32_t),1,Diffusion);
    fclose(Diffusion);
    
    free(str);
}

void printMatrixPart(){
    uint64_t xlower = binmax/2; uint64_t xupper = binmax/2+1;
    uint64_t ylower = binmax/2; uint64_t yupper = binmax/2+1;
    uint64_t zlower = binmax/2; uint64_t zupper = binmax/2+1;
    int stotal = 9*N*N;
    uint64_t i,j,k,l;
    FILE *Diffusion;
    char* str = malloc(sizeof(char)*50);
    sprintf(str, "./output/Matrix_%d.txt", NP);
    Diffusion = fopen(str, "w");
    for(i = xlower; i<xupper; ++i){
        for(j = ylower; j<yupper; ++j){
            for(k = zlower; k<zupper; ++k){
                fprintf(Diffusion,"%d ",matrix->gridCounting[i*binmax*binmax+j*binmax+k]);
            }
        }
    }
    for(i = xlower; i<xupper; ++i){
        for(j = ylower; j<yupper; ++j){
            for(k = zlower; k<zupper; ++k){
                uint64_t start = i*binmax*binmax*9*N*N + j*binmax*9*N*N + k*9*N*N;
                fprintf(Diffusion,"\ndx=%lld~%lld dy=%lld~%lld dz=%lld~%lld\n",2*i-binmax,2*i-binmax+2,2*j-binmax,2*j-binmax+2,2*k-binmax,2*k-binmax+2);
                if(matrix->gridCounting[i*binmax*binmax+j*binmax+k]!=0){
                    for(l = 0; l<9*N*N; ++l){
                        double value = matrix->DiffMatrix[start+l]/matrix->gridCounting[i*binmax*binmax+j*binmax+k];
                        fprintf(Diffusion,"%lf ", value);
                    }
                } else{
                    for(l = 0; l<9*N*N; ++l){
                        double value = matrix->DiffMatrix[start+l];
                        fprintf(Diffusion,"%lf ", value);
                    }
                }
            }
        }
    }
    fprintf(Diffusion,"\n\nSelf:\n");
    for(i = 0; i<stotal; ++i){
        fprintf(Diffusion, "%lf ", matrix->selfMatrix[i]/step/NP);
    }
    fclose(Diffusion);
    free(str);
}

Vector3D_t getNID(double dx,double dy,double dz){
    dx -= round(dx/L)*L;
    dy -= round(dy/L)*L;
    dz -= round(dz/L)*L;
    
    Vector3D_t NID;
    NID.x = dx;
    NID.y = dy;
    NID.z = dz;
    
    return NID;
}

bin_t binCoM(Vector3D_t dCoM_ij){
    bin_t bin;
    
    double dx = dCoM_ij.x;
    double dy = dCoM_ij.y;
    double dz = dCoM_ij.z;
    
    bin.x = floor(dx/BIN_SIZE)+binmax/2;
    bin.y = floor(dy/BIN_SIZE)+binmax/2;
    bin.z = floor(dz/BIN_SIZE)+binmax/2;
    
    return bin;
}

/*bin_t binCoMFine(int i, int j){
    bin_t bin;
    
    double dx = distCoM[i][j].dx;
    double dy = distCoM[i][j].dy;
    double dz = distCoM[i][j].dz;
    
    bin.x = floor(dx/binsize)+bmax/2;
    bin.y = floor(dy/binsize)+bmax/2;
    bin.z = floor(dz/binsize)+bmax/2;
    
    return bin;
}*/



