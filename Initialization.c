//
//  Initialization.c
//  StaticComplete
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "Initialization.h"
#include "Parameters.h"
#include "Matrix.h"
#include "util.h"

void Initialization(int argc, char * argv[]){
    ParseInput(argc, argv);
    
    matrix = (Matrix_t*)malloc(sizeof(Matrix_t));
    gwParm = (gwParm_t*)malloc(sizeof(gwParm_t));
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    trajStep = 10000;
        
    //isLinearize = false;
    //if(isLinearize) printf(">Using trilinear interpolation\n");
    
    MatrixCalculation = true;
    
    directory = malloc(50*sizeof(char));
    directory = "./output";
    
    if(HImode){
        readIteration();
        readMatrixAvg();
        GeyerWinterCalculation();
        
        //fineMatrix();
        //printMatrixFine();
        //free(matrix->DiffMatrixAvg);

        iteration++;
        if(EwaldSum) printf(">Runing iteration #%d HI with Ewald Summated RPY...\n",iteration);
        else printf(">Runing iteration #%d HI with Regular RPY...\n",iteration);
        
        printIteration();
    } else{
        iteration = 1;
        if(EwaldSum) printf(">Runing iteration #1 FD calcualting HI with Ewald Summated RPY...\n");
        else printf(">Runing iteration #1 FD calcualting HI with Regular RPY...\n");
        
        printIteration();
    }
    
    outputName = malloc(100*sizeof(char));
    sprintf(outputName,"%s/Static_%u.txt",directory,NP);
    trajName = malloc(100*sizeof(char));
    sprintf(trajName,"%s/Static_%u_%u.xyz",directory,NP,iteration);
    
    ReadInitFromFile = false;
    
    EwaldSum = true;
    
    CheckOverlap = false; //change it to false if system is too dense

    MSD_CoM = 0.0;
    
    if(MatrixCalculation) printf(">Calculating Averaged Matrix every 100 time steps\n");
    else printf(">Not Calculating Averaged Matrix\n");
    
    AVGmode = 0;
    /* Average Methods: *
     * use AVGmode = 0 if don't want to average the Matrix calculated in current iteration with any of the previous matrices
     * use AVGmode = 1 if want to averaged the current matrix with the previous matrices
     * The averaging algorism:
     ** mm'[n] = AVGfactor * mm[n] + (1 - AVGfactor) * mm'[n-1] **
     * Specify AVGfactor below or use initial value 0.5 *
     */
    
    AVGfactor = 0.25;
    
    if(AVGmode) printf("*Averaging with a factor of: %.2lf*\n",AVGfactor);
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    idum = malloc(sizeof(long));
    *idum = -1;
    ran1(idum);
    *idum = -1*initRan();

    //generateOutput();
    
    initParaEvir();
    
    initChains();
    
    if(ReadInitFromFile) readChains();
    if(EwaldSum) initM2();
    
    initMatrix();
    
    step = 0;
    ext = 0;
    
    //printInitial();
}

void generateOutput(){
    FILE *outputfile;
    outputfile = fopen(outputName, "w");
    fprintf(outputfile, "epsilon = %lf\n", EPSILON);
    fprintf(outputfile, "kappa = %lf\n", KAPPA);
    fprintf(outputfile, "N = %d\n", N);
    fprintf(outputfile, "dt = %lf\n", DT);
    fprintf(outputfile, "tmax = %d\n", TMAX);
    fprintf(outputfile, "NP = %d\n", NP);
    fprintf(outputfile, "L = %lf\n", L);
    fprintf(outputfile, "KMAX = %d\n", KMAX);
    fprintf(outputfile, "BinSize = %f\n", BIN_SIZE);
    fprintf(outputfile, "c_star = %lf\n\n",C_STAR);
    fprintf(outputfile, "SEED %ld\n",*idum);
    fclose(outputfile);
}

void readChains(){
    Chain = calloc(N*NP,sizeof(Chain_t));
    
    printf(">Reading Conformations from File...\n");
    
    FILE *trajec;
    trajec = fopen("./output/Trajectory.xyz","r");
    if(trajec==NULL){
        printf("!!Trajectory File not Founded, please use .xyz file. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i<N*NP; ++i){
        fscanf(trajec, "A %lf %lf %lf\n", &Chain[i].rx, &Chain[i].ry, &Chain[i].rz);
    }
    fclose(trajec);
    
    printf(">Read Conformation Sucessfully.\n");
}

void readIteration(){
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/iteration_%s.bin",directory,c_normalchar);
    FILE *iter = fopen(str,"rb");
    if(iter==NULL){
        printf("!!Iteration File not Founded. Exiting...\n");
        exit(EXIT_FAILURE);
    } else{
        printf(">Reading Iteration...\n");
    }
    fread(&iteration,sizeof(uint16_t),1,iter);
    fclose(iter);
    free(str);
}

void readMatrixAvg(){
    FILE *Matrix;
    char* str = malloc(100*sizeof(char));
    printf(">Reading from pre-averaged Matrix ...\n");
    if(MatrixCalculation){
        sprintf(str,"%s/Matrix_%s_%u.bin",directory,c_normalchar,iteration);
    } else{
        sprintf(str,"%s/Matrix_%s.bin",directory,c_normalchar);
    }
    
    Matrix = fopen(str, "rb");
    if(Matrix==NULL){
        printf("!!Matrix File not Founded. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    
    uint32_t Ntotal,Nc;
    uint64_t binX;
    double binsize;
    //header
    fread(&Ntotal,sizeof(uint32_t),1,Matrix);
    fread(&Nc,sizeof(uint32_t),1,Matrix);
    fread(&binX,sizeof(uint64_t),1,Matrix);
    fread(&binsize,sizeof(double),1,Matrix);
    //matrix
    uint64_t bintotal = binX*binX*binX;
    int stotal = 9*Nc*Nc;
    uint64_t mtotal = bintotal*stotal;
    
    printf(">total N:%d Nc:%d\n>total Bin:%llu total Elements:%llu\n",Ntotal,Nc,bintotal,mtotal);
    
    matrix->DiffMatrixAvg = (float*)malloc(mtotal*sizeof(float));
    matrix->selfMatrixAvg = (float*)malloc(stotal*sizeof(float));
    matrix->counter = (uint32_t*)malloc(bintotal*sizeof(uint32_t));
    
    //Original Structure
    fread(matrix->DiffMatrixAvg,sizeof(float),mtotal,Matrix);
    fread(matrix->selfMatrixAvg,sizeof(float),stotal,Matrix);
    fread(matrix->counter, sizeof(uint32_t),bintotal,Matrix);
    fclose(Matrix);
    free(str);
    
    printf(">Read Matrix Sucessfully.\n");
}

void readGeyerWinterParameters(){
    FILE* gwParameter;
    char *str = malloc(100*sizeof(char));
    sprintf(str,"%s/GeyerWinterParam_%d.bin",directory,iteration);
    gwParameter = fopen(str, "rb");
    
    gwParm->C_i = (float*)malloc(3*N*sizeof(float));
    fread(&gwParm->beta_ij,sizeof(double),1,gwParameter);
    fread(gwParm->C_i,sizeof(float),3*N,gwParameter);
    fclose(gwParameter);
    free(str);
}

void initParaEvir(){
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
    //numThreads = omp_get_num_procs();
    printf(">OpenMP Detected, Runing with %d threads\n",numThreads);
#else
    printf(">OpenMP not detected\n");
#endif
}

void initChains(){
    Chain = calloc(N*NP,sizeof(Chain_t));
    
    printf(">Initializing the Initial Conformation...\n");
    
    for(int i = 0; i<NP; ++i){
        int test = 0;
        while(test==0){
            test = 1;
            Chain[i*N].rx = ran1(idum)*L - L/2.0;
            Chain[i*N].ry = ran1(idum)*L - L/2.0;
            Chain[i*N].rz = ran1(idum)*L - L/2.0;
            
            Chain[i*N].rx -= round(Chain[i*N].rx/L)*L;
            Chain[i*N].ry -= round(Chain[i*N].ry/L)*L;
            Chain[i*N].rz -= round(Chain[i*N].rz/L)*L;
            
            for(int j = 0; j<i; ++j){
                double dx = Chain[i*N].rx - Chain[j*N].rx;
                double dy = Chain[i*N].ry - Chain[j*N].ry;
                double dz = Chain[i*N].rz - Chain[j*N].rz;
                dx -= round(dx/L)*L;
                dy -= round(dy/L)*L;
                dz -= round(dz/L)*L;
                
                if(dx*dx+dy*dy+dz*dz<2.0){
                    test = 0;
                }
            }
        }
    }
    
    for(int i = 0; i<NP; ++i){
        for(int j = 1; j<N; ++j){
            int test = 0;
            while(test==0){
                test = 1;
                double theta = ran1(idum)*2.0*3.14158;
                double phi = acos(2.0*ran1(idum)-1.0);
                Chain[i*N+j].rx = Chain[i*N+j-1].rx + 2.05*cos(theta)*sin(phi);
                Chain[i*N+j].ry = Chain[i*N+j-1].ry + 2.05*sin(theta)*sin(phi);
                Chain[i*N+j].rz = Chain[i*N+j-1].rz + 2.05*cos(phi);
                Chain[i*N+j].rx -= round(Chain[i*N+j].rx/L)*L;
                Chain[i*N+j].ry -= round(Chain[i*N+j].ry/L)*L;
                Chain[i*N+j].rz -= round(Chain[i*N+j].rz/L)*L;
                
                if(CheckOverlap){
                    int k;
                    for(k = 0; k<i*N+j; ++k){
                        double dx = Chain[i*N+j].rx-Chain[k].rx;
                        double dy = Chain[i*N+j].ry-Chain[k].ry;
                        double dz = Chain[i*N+j].rz-Chain[k].rz;
                        dx -= round(dx/L)*L;
                        dy -= round(dy/L)*L;
                        dz -= round(dz/L)*L;
                        
                        if(dx*dx+dy*dy+dz*dz<4.0){
                            test = 0;
                        }
                    }
                }
            }
        }
    }
}

void initMatrix(){
    binmax = ceil(L/BIN_SIZE); //max number of bins
    matrix->DiffMatrix = calloc(9*N*N*binmax*binmax*binmax,sizeof(float));
    matrix->selfMatrix = calloc(9*N*N,sizeof(float));
    matrix->gridCounting = calloc(binmax*binmax*binmax,sizeof(uint32_t));
}

void initM2(){
    int m2d = KMAX*2.0 + 1.0;
    m2 = calloc(m2d, sizeof(m2[0]));
    for(int i = 0; i < m2d; i++){
        m2[i] = calloc(m2d, sizeof(m2[0][0]));
        for(int j = 0; j < m2d; j++){
            m2[i][j] = calloc(m2d, sizeof(m2[0][0][0]));
            for(int k = 0; k < m2d; k++){
                m2[i][j][k] = calloc(9, sizeof(m2[0][0][0][0]));
            }
        }
    }
    
    double alpha = 6.0/L;
    
    for(int kx = -KMAX; kx<KMAX+1; ++kx){
        double rkx = 2.0*M_PI*kx/L;
        for(int ky = -KMAX; ky<KMAX+1; ++ky){
            double rky = 2.0*M_PI*ky/L;
            for(int kz = -KMAX; kz<KMAX+1; ++kz){
                double rkz = 2.0*M_PI*kz/L;
                double rkk = rkx*rkx+rky*rky+rkz*rkz;
                double mm3 = (1.0-rkk/3.0)*(1.0+rkk*pow(alpha,-2.0)/4.0+rkk*rkk*pow(alpha,-4.0)/8.0)*6.0*M_PI/rkk*exp(-rkk*pow(alpha,-2.0)/4.0);
                if(kx==0&&ky==0&kz==0) mm3 = 0;
                
                //printf("kx %d,ky %d,kz %d,mm3:%lf\n",kx,ky,kz,mm3);
                
                m2[kx+KMAX][ky+KMAX][kz+KMAX][0] += mm3*(1.0-rkx*rkx/rkk);
                m2[kx+KMAX][ky+KMAX][kz+KMAX][1] += -mm3*rkx*rky/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][2] += -mm3*rkx*rkz/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][3] += -mm3*rky*rkx/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][4] += mm3*(1.0-rky*rky/rkk);
                m2[kx+KMAX][ky+KMAX][kz+KMAX][5] += -mm3*rky*rkz/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][6] += -mm3*rkz*rkx/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][7] += -mm3*rkz*rky/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][8] += mm3*(1.0-rkz*rkz/rkk);
            }
        }
    }
}

void printInitial(){
    FILE *Trajectory;
    Trajectory = fopen(trajName, "w");
    fprintf(Trajectory, "%d\n-1\n", N*NP);
    for(int i = 0; i<N*NP; ++i){
        fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
    }
    fclose(Trajectory);
}

void printIteration(){
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/iteration_%s.bin",directory,c_normalchar);
    FILE *iterat = fopen(str,"wb");
    printf(">Generating Iteration File...\n");
    fwrite(&iteration,sizeof(uint16_t),1,iterat);
    fclose(iterat);
    free(str);
}

long initRan(){
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000+1000000000;
    
    /*time_t seconds;
    time(&seconds); //switched from time due to crashes... unsure why...
    return -1*(unsigned long)(seconds); //Dividing by 100 keeps within correct bounds? not sure why this works */
}

double ran1(long *idum){
    long j;
    long k;
    static long idum2 = 123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if(*idum <= 0){
        if(-(*idum)<1) *idum=1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for(j=NTAB+7;j>=0;--j)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum<0) *idum+=IM1;
            if(j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum<0) *idum += IM1;
    k=idum2/IQ2;
    if(*idum<0) idum2+= IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if(iy<1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}




