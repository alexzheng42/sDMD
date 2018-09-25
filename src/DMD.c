//
//  New Collision
//
//  Created by Size Zheng on 6/12/13.
//  Copyright (c) 2013 Size Zheng. All rights reserved.
//

#include "DMD.h"

//-----------------
//directory
char datadir[1024] = ".";
char savedir[1024] = ".";
//-----------------


//-----------------
//basic variables
int atomnum;
int **connectionMap;
int chargeAA[9] = {0};
long int frame = 0;
unsigned int seed;
double timestep, currenttime = 0;
double oldcurrenttime = 0;
double cutoffr;
double instTemperature, targetTemperature, oldTemperature = 0;
double TotlKinEner;
double viscosity; //dynamic u/nm/ps
double outputrate;
double outputrecord = 0;
double boxDimension[4]; //unit A
struct AtomStr *atom;
struct AAStr *aminoacid;
struct PepStr *protein;
char neworcontinue[20];
char Methodtype[20];
char temperatureType[20];
//-----------------


//-----------------
//system check
long int countReCal = 0;
long int collisioneventsum = 0;
long int HBeventsum = 0;
long int bondeventsum = 0;
long int HBNeighboreventsum = 0;
long int thermostateventsum = 0;
long int pbcandcrosseventsum = 0;
long int walleventsum = 0;
long int oldtotaleventsum = 0;
long int newtotaleventsum = 0;
//-----------------


//-----------------
//protein sequence
int numofprotein = 0;
int totalAminoAcids = 0;
//-----------------


//-----------------
//subcell property
int cellnum[4];
int *celllist;
double cellsize[4];
//-----------------


//-----------------
//binary tree
struct EventListStr CBT = {.count = 0};
int nthCheck, nthNode;
//-----------------


//-----------------
//potential pair
struct ConstraintStr potentialPairCollision[NATOMTYPE + 1][NATOMTYPE + 1];
struct ConstraintStr potentialPairHB[11][NATOMTYPE + 1][NATOMTYPE + 1];
struct HBPotentialStr HBPotential;


//-----------------
//thermostat and solvent
char thermostatType[20];
double thermoF;
//-----------------


//-----------------
//record parameters
int HBnumformed = 0;
int alphaHBformed = 0;
char timer[30];
//-----------------


//-----------------
//multi-thread variables
int threadNum;
int codeNum;
int eventToCommit;
signed int readCount = 0, writeCount = 0, activeWrite = 0;
struct ThreadStr **thread;
struct RandomListStr *threadRandomList;
struct ThreadInfoStr *thrInfo;
struct PreCalObjStr *preCalList;
struct PreCalObjStr **preCalPtr;
pthread_t *thread_t;
pthread_mutex_t mainFrameLock, rawDataLock, commitLock;
pthread_mutex_t mstThrLock, slvThrLock;
pthread_cond_t mainFrameRenew = PTHREAD_COND_INITIALIZER;
pthread_cond_t rawDataRenew = PTHREAD_COND_INITIALIZER;
pthread_cond_t readCheck = PTHREAD_COND_INITIALIZER;
pthread_cond_t writeCheck = PTHREAD_COND_INITIALIZER;
pthread_cond_t mstThrQ = PTHREAD_COND_INITIALIZER;
pthread_cond_t slvThrQ = PTHREAD_COND_INITIALIZER;


//-----------------
//REMD
struct REMDStr REMDInfo = {.flag            = 0,
                           .REMD_PortNum    = -1,
                           .REMD_T          = -1,
                           .REMD_OutputRate = 1000,
                           .REMD_ServerName = "",
                           .REMD_ExtraName  = ""};


//-----------------
//Wall
struct AtomStr *wall = NULL;
struct DynamicWallStr wallDyn = {.mark = "no",
                                 .sign = 0,
                                 .size = 0,
                                 .rate = 0,
                                 .step = 0,
                                 .touch = 0,
                                 .curTime = 0,
                                 .origBoxDim = {0}};
struct ObstructionStr obstObj = {.mark = 0,
                                 .num  = 0};
struct TunnelStr      tunlObj = {.mark = 0};
struct SphObstStr  SphObstObj = {.mark = 0};

char wallExist[20];
char wallType[20];
//-----------------


//-----------------
//Flow
struct FlowStr flow = {.mark           = 0,
                       .constV.v       = {0, 0, 0, 0},
                       .force.f        = {0, 0, 0, 0},
                       .charge.num     = 0,
                       .charge.PBCMark = 1,
                       .charge.potGrad = 0,
                       .charge.gap     = 0
};
//-----------------


//-----------------
//files for saving data
struct FileStr *fileList;
//-----------------
//To maximum speed, open and close files as less as possible.


//-----------------
//other
int visual = 0;
int tmpInt;
double tmpDouble;
//-----------------




int main(int argc, const char * argv[]) {
    clock_t st, et;
    double hours, mins, sec;
    double oldTime, oldTimeEvaluation, timeEvaluation; //based on 100 ps
    char directory[200];
    FILE *timeFile;
    
    st = clock();
    
    InputData(argc, argv); //input data from .txt files in input folder
    
    //========================================================================================
    //check parameter
    oldTime = currenttime;
    printf("\n");
    printf("\n");
    printf("--------------------------------\n");
    printf("| Simulation configurations:   |\n");
    printf("| currenttime=%16.2f |\n", currenttime);
    printf("| timestep=%19.2f |\n", timestep); //timestep -> total time step to run
    printf("| targetT=%20.2f |\n", targetTemperature);
    printf("| OutputRate=%17.2f |\n", outputrate);
    printf("| CutOffRadius=%15.2f |\n", cutoffr);
    printf("| Method=%21s |\n", Methodtype);
    printf("| Thermostat=%17s |\n", thermostatType);
    printf("| WallExist=%18s |\n", wallExist);
    printf("| WallType=%19s |\n", wallType);
    printf("| boxsize=%6.1f%7.1f%7.1f |\n", boxDimension[1], boxDimension[2], boxDimension[3]);
    printf("--------------------------------\n");
    printf("Please check the Log file for detail information!\n");
    printf("\n");
    
    sprintf(directory, "%s/timeEvaluation.dat", savedir);
    if (strncmp("new", neworcontinue, 1) == 0 || fopen(directory,"r") == NULL) {
        printf("Time evaluation is unavailable this time!\n");
        oldTimeEvaluation = 0;
    } else {
        timeFile = fopen(directory, "rb");
        fread(&oldTimeEvaluation, sizeof (double), 1, timeFile);
        timeEvaluation = (timestep - currenttime) / 100 * oldTimeEvaluation;
        hours = floor(timeEvaluation / 3600);
        mins = floor((timeEvaluation - hours * 3600) / 60);
        sec = timeEvaluation - hours * 3600 - mins * 60;
        printf("This simulation may take %.0lf hours %.0lf mins %.0lf sec\n", hours, mins, sec);
        fclose(timeFile);
    }
    //========================================================================================
    
    printf("\nSimulation starts...\n");
    ThreadProcess();
    
    //========================================================================================
    printf("\nFinished!\n");
    printf("\n");
    printf("\n");
    
    //statistic information about the finished simulation (continuous)
    newtotaleventsum = collisioneventsum   +
                       bondeventsum        +
                       HBeventsum          +
                       HBNeighboreventsum  +
                       thermostateventsum  +
                       pbcandcrosseventsum +
                       walleventsum;
    
    printf("ReCalculate times   =%li\n\n", countReCal);
    printf("collision times     =%li\n", collisioneventsum);
    printf("bond times          =%li\n", bondeventsum);
    printf("HB times            =%li\n", HBeventsum);
    printf("HB Neighbor times   =%li\n", HBNeighboreventsum);
    printf("thermostat times    =%li\n", thermostateventsum);
    printf("PBC&CC times        =%li\n", pbcandcrosseventsum);
    printf("Wall times          =%li\n", walleventsum);
    printf("percentage of thermostat=%.4f%%\n", (float) thermostateventsum / newtotaleventsum * 100);
    
    et = clock();    //record the ending time
    hours = floor((double) (et - st) / CLOCKS_PER_SEC / 3600);
    mins = floor(((double) (et - st) / CLOCKS_PER_SEC - hours * 3600) / 60);
    sec = (double) (et - st) / CLOCKS_PER_SEC - hours * 3600 - mins * 60;
    printf("\nTotal time taken: %.0f hours %.0f mins %.0f sec\n", hours, mins, sec);
    printf("Calculate event times: %li\n", newtotaleventsum - oldtotaleventsum);
    
    timeFile = fopen(directory, "wb");
    sec += hours * 3600 + mins * 60;
    timeEvaluation = (sec + oldTimeEvaluation * oldTime / 100) / timestep * 100;
    fwrite(&timeEvaluation, sizeof (double), 1, timeFile);
    fclose(timeFile);
    //========================================================================================
    
    GlobalCloseFree();
    return 0;
}

