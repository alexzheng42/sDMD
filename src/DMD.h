//
//  DMD.h
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#ifndef DMD_h
#define DMD_h

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "List.h"

//#define DEBUG_RANDOM
//#define DEBUG_PRINTF
//#define DETAILS
//#define GEL
#define HYDROGEN_BOND

#define TIME_GAP      0.5
#define INFINIT       1E8
#define ZERO          1E-8
#define INVALID      -111
#define EXTEND_RATIO  0.95 //make sure the cut off radius is less than the subcell length
#define THERMOSTAT(method) ((strncmp(method,"Andersen",1) == 0) ? (5.25) : (0.10)) //keep the frequency about 1%

#ifdef DEBUG_RANDOM
#define RANDOM_SEED 1522672160
#endif

/*
 Units of the system
 Temperature: kcal / (mol * kB) around 502 K (kB = 0.001987204118 kcal / (mol * K))
        Mass: u (unified atomic mass unit)
      Energy: kcal / mol
      Length: A
        Time: 1 time unit is approximately 48.88 fs
 */
//------------------


#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)


//connection map bit
#define BOND_CONNECT             0b00001
#define CONSTRAINT_CONNECT       0b00010
#define HB_CONNECT               0b00100
#define NEIGHBOR_CONNECT         0b01000
#define OVERLAP_CONNECT          0b10000


#define INT_CALLOC(name, size) \
name=(int *)calloc(size, sizeof(int));

#define INT_2CALLOC(name, size1, size2, NN)    \
name=(int **)calloc(size1, sizeof(int *));   \
for (NN=0; NN<size1; NN++) {    \
name[NN]=(int *)calloc(size2, sizeof(int)); \
}

#define DOUBLE_CALLOC(name, size) \
name=(double *)calloc(size, sizeof(double));

#define DOUBLE_2CALLOC(name, size1, size2)    \
name=(double **)calloc(size1, sizeof(double *));   \
for (int NN=0; NN<size1; NN++) {    \
name[NN]=(double *)calloc(size2, sizeof(double)); \
}

#define DOUBLE_3CALLOC(name, size1, size2, size3)    \
name=(double ***)calloc(size1, sizeof(double **));  \
for (int NN=0; NN<size1; NN++) {    \
name[NN]=(double **)calloc(size2, sizeof(double *));    \
for (int NN2=0; NN2<size2; NN2++) { \
name[NN][NN2]=(double *)calloc(size3, sizeof(double));  \
}   \
}

#define FREE2(name, size2)    \
for (int NN=0; NN<size2; NN++) {    \
free(name[NN]); \
}   \
free(name);

#define FREE3(name, size2, size3)   \
for (int NN=0; NN<size2; NN++) {    \
for (int NN2=0; NN2<size3; NN2++) {    \
free(name[NN][NN2]);    \
}   \
free(name[NN]); \
}   \
free(name);

#define MIN3(num1, num2, num3)    \
(((num1) < (num2)) ? (((num1) < (num3)) ? (num1) : (num3)) :    \
(((num2) < (num3)) ? (num2) : (num3)))

#define MAX3(num1, num2, num3)    \
(((num1) > (num2)) ? (((num1) > (num3)) ? (num1) : (num3)) :    \
(((num2) > (num3)) ? (num2) : (num3)))

#define MIN2(num1, num2)    \
(((num1) < (num2)) ? (num1) : (num2))

#define TRANSFER_VECTOR(recv, send)   \
recv[0]=send[0];    \
recv[1]=send[1];    \
recv[2]=send[2];    \
recv[3]=send[3];

#define POINT_TO_STRUCT(pointer, struct)  \
pointer=&struct[0];

#define ABSVALUE(num) ((num<0) ? (-num) : (num))

#define DOT_MINUS(n1, n2, out) \
out[1]=n1[1]-n2[1]; \
out[2]=n1[2]-n2[2]; \
out[3]=n1[3]-n2[3];

#define DOT_PLUS(n1, n2, out) \
out[1]=n1[1]+n2[1]; \
out[2]=n1[2]+n2[2]; \
out[3]=n1[3]+n2[3];

#define CROSS_PROD(n1, n2, result)   \
result[1]=n1[2]*n2[3]-n1[3]*n2[2];  \
result[2]=n1[3]*n2[1]-n1[1]*n2[3];  \
result[3]=n1[1]*n2[2]-n1[2]*n2[1];

#define FACTOR_PROD(factor, vector, out) \
out[1]=factor*vector[1];    \
out[2]=factor*vector[2];    \
out[3]=factor*vector[3];

#define DOT_PROD(n1, n2) \
n1[1]*n2[1]+n1[2]*n2[2]+n1[3]*n2[3];

#define FIND_DUPLICATED(list1, list2, num, flag)  \
flag=0; \
num=0;    \
while (list2[++num]) {   \
if (list1==list2[num]) {   \
flag=1; \
break;  \
}    \
}

#define FIND_DUPLICATED_IN_THE_RANGE(list1, list2, start, end, num, flag)  \
flag=0; \
for (num=start; num<=end; num++) {  \
if (list2[num]==0) {  \
break;  \
}  else if (list1==list2[num]) {   \
flag=1; \
break;  \
}   \
}   \

#define FIND_DUPLICATED_AND_ADD(list1, list2, num, flag)  \
flag=0; \
num=0;    \
while (list2[++num]) {   \
if (list1==list2[num]) {   \
flag=1; \
break;  \
}    \
}   \
if (flag==0) list2[num]=list1;

#define FIND_DUPLICATED_IN_THE_RANGE_AND_ADD(list1, list2, start, end, num, flag)  \
flag=0; \
for (num=start; num<=end; num++) {  \
if (num==end) { \
printf("List space may be not enough!\n");   \
flag=1; \
}   \
if (list2[num]==0) {  \
break;  \
}  else if (list1==list2[num]) {   \
flag=1; \
break;  \
}   \
}   \
if (flag==0 && list1!=0) list2[num]=list1;

#define READ_MUTEX_PRE(lock, writers, readersQ, readers) \
pthread_mutex_lock(&lock); \
while (!(writers == 0)) \
pthread_cond_wait(&readersQ, &lock); \
readers ++;    \
pthread_mutex_unlock(&lock);

#define READ_MUTEX_POST(lock, writersQ, readers)    \
pthread_mutex_lock(&lock);    \
if (--readers == 0)    \
pthread_cond_signal(&writersQ);    \
pthread_mutex_unlock(&lock);

#define WRITE_MUTEX_PRE(lock, readers, writersQ, writers, activeWriters)    \
pthread_mutex_lock(&lock);    \
writers ++;    \
while (!((readers == 0) && (activeWriters == 0)))    \
pthread_cond_wait(&writersQ, &lock);    \
activeWriters ++;    \
pthread_mutex_unlock(&lock);

#define WRITE_MUTEX_POST(lock, writersQ, readersQ, writers, activeWriters)    \
pthread_mutex_lock(&lock);    \
writers --;    \
activeWriters --;    \
if (writers > 0)    \
pthread_cond_signal(&writersQ);    \
else    \
pthread_cond_broadcast(&readersQ);    \
pthread_mutex_unlock(&lock);

enum EEventType {
    Coli_Event,  //collision event
    Bond_Event,  //bond event
      HB_Event,  //HB event
    HBNe_Event,  //HB neighbor event
    Ther_Event,  //thermostat event
    Wall_Event,  //wall event
    Obst_Event,  //obstruction event
    Tunl_Event,  //tunnel event
    Chrg_Event,  //charge event
    SphO_Event,  //column event
     PCC_Event,  //PBC and crossing event
    Lagv_Event,  //Langevin event
     TBD_Event,  //TBD event
    Cacl_Event,  //event needs to be cancelled
    Invd_Event   //invalid event, for initialization
};

enum EInteractionType {
    CoreColli,
    WellBounc,
    WellCaptr,
    WellEscap,
    TBDCaptr,
    TBDEscap,
    TBDHBForm,
    TBDHBBrik,
    Invalid_Inter
};

enum EHBEventType {
    NoEvent,
    HBForm,
    HBBreak,
    HBSkip
};

enum EChargeAA {
    TER,
    ARG,
    ASP,
    GLU,
    HIS,
    LYS,
    AlC,
    AlP,
    COM
};

enum FileType {
    trj,         //gro trajectory
    pot,         //potential energy
    kin,         //kinetic energy
    tem,         //temperature
    cnt,         //connection map
    HBn,         //HB number
    xyz,         //xyz trajectory
    pdb,         //pdb trajectory
    lgf,         //log file
    sysInfo,     //saved data for analyzing
    savedData,   //saved data for continuity
    RE,          //replica exchange data
    lenFileType  //length of this list
};

struct HBPotentialStr //B: Backbone, S: Sidechain
{
    double BB;
    double SS;
    double BS;
};

struct AtomSeqStr {
    int atomNum;
    int proteinNum;
    int aminoacidNum;
};

struct StepPotenStr {
    double d;
    double e;
    double accumulated;
    struct StepPotenStr *next;
};

struct ConstraintStr {
    int connection;
    double dmin;
    struct StepPotenStr *step;
    double dmax;
    struct ConstraintStr *next;
};

struct HBStr {
    int bondConnection;
    int neighbor;
    char role;
    enum EHBEventType interactionType;
};

struct HBNeighborStr{
    int neighborPartner[5];
    int partnerNum;
    int neighborStatus;
};

struct InteractionEventStr {
    int partner;
    long counter;
    long partnerCounter;
    double time;
    double potential;
    enum EEventType eventType;
    enum EInteractionType subEventType;
};

struct PropertyStr {
    char name[5];
    char extraProperty[2][10];
    char nameOfAA[5];
    int num;
    int type;
    double charge;
    double mass;
    double color[3];
    struct AtomSeqStr sequence;
    struct ConstraintStr *bond;
    struct ConstraintStr *constr;
};

struct DynamicStr {
    int cellIndex[4];
    double coordinate[4];
    double velocity[4];
    struct HBStr HB;
    struct HBNeighborStr HBNeighbor;
    struct InteractionEventStr event;
};

struct EventListStr {
    int count;
    int *node;
    int *leaf;
    double **time;
};

struct AtomStr {
    struct PropertyStr *property;   //static properties of the atom, such as mass, atom # ...
    struct DynamicStr *dynamic;    //dynamic properties of the atom, such as coordinate, velocity ...
};

struct AAStr {
    char nameOfAA[5];
    int startAtomNum;
    int endAtomNum;
    int proteinNum;
    double mass;
};

struct PepStr {
    int startAtomNum;
    int endAtomNum;
    int startAANum;
    int endAANum;
    double mass;
};

struct FileStr {
    int mark;
    char name[256];
    FILE *file;
};

struct ThreadStr {
    int tid;
    int atomNum;
    int getWork, finishWork; //0: no; 1: yes
    int hazardStatus;
    double getWorkTime;
    double finishWorkTime;
    double getWorkFrame;
    struct AtomStr **raw;       //point to the raw atom data
    struct AtomStr **listPtr;   //point to the thread atom copy
    struct AtomStr *newTarget;  //point to the obj in atomList
    struct AtomStr *newPartner; //point to the obj in atomList
    struct AtomStr oldTarget;   //store the original target data
    struct AtomStr oldPartner;  //store the original partner data
    struct FileStr *fileList;
    list atomList;  //thread atom copy
};


struct ThreadInfoStr {
    int threadID;
    int *threadRenewList;
};


struct PreCalObjStr {
    int eventStatus;
    int renewList[4];
    struct ThreadStr *data;
};


struct REMDStr {
    int flag;
    int REMD_PortNum;
    int REMD_OutputRate;
    double REMD_T;
    char REMD_ServerName[50];
    char REMD_ExtraName[50];
};


struct DynamicWallStr {
    char mark[16];
    int sign;
    int touch;
    double size;
    double rate;
    double step;
    double curTime;
    double origBoxDim[4];
};


struct ObstructionStr {
    int mark;
    int num;
    int *hole;
    double **position;
    struct AtomStr *obst;
};


struct SphObstStr {
    int mark;
    int num;
    double *radius;
    double **position;
};


struct TunnelStr {
    int mark;
    int num;
    double **position;
    double *diameter;
    struct AtomStr tunnel;
};


struct ConstVStr {
    double v[4];
};


struct ForceStr {
    double f[4];
    double *timeRec;
    double **a;
};


struct ChargeStr {
    int num;
    int PBCMark;
    double **position;
    double **velocity;
    double *potGrad;
    double *gap;
};


struct FlowStr {
    int mark;
    struct ConstVStr constV;
    struct ForceStr force;
    struct ChargeStr charge;
};


//--------------------------------------------
//Global Variables
//--------------------------------------------


//-----------------
//directory
extern char datadir[1024];
extern char savedir[1024];
//-----------------


//-----------------
//basic variables
extern int atomnum;
extern int **connectionMap;
extern int chargeAA[9];
extern long int frame;
extern unsigned int seed;
extern double timestep, currenttime;
extern double oldcurrenttime;
extern double cutoffr;
extern double instTemperature, targetTemperature, oldTemperature;
extern double TotlKinEner;
extern double viscosity; //dynamic u/nm/ps
extern double outputrate;
extern double outputrecord;
extern double boxDimension[4]; //unit A
extern struct AtomStr *atom;
extern struct AAStr *aminoacid;
extern struct PepStr *protein;
extern char neworcontinue[20];
extern char Methodtype[20];
extern char temperatureType[20];
//-----------------


//-----------------
//system check
extern long int countReCal;
extern long int collisioneventsum;
extern long int HBeventsum;
extern long int bondeventsum;
extern long int HBNeighboreventsum;
extern long int thermostateventsum;
extern long int pbcandcrosseventsum;
extern long int walleventsum;
extern long int oldtotaleventsum;
extern long int newtotaleventsum;
//-----------------


//-----------------
//protein sequence
extern int numofprotein;
extern int totalAminoAcids;
//-----------------


//-----------------
//subcell property
extern int cellnum[4];
extern int *celllist;
extern double cellsize[4];
//-----------------


//-----------------
//binary tree
extern struct EventListStr CBT;
extern int nthCheck, nthNode;
//-----------------


//-----------------
//potential pair
extern struct ConstraintStr potentialPairCollision[32][32];
extern struct ConstraintStr potentialPairHB[11][32][32];
extern struct HBPotentialStr HBPotential;


//-----------------
//thermostat and solvent
extern char thermostatType[20];
//-----------------


//-----------------
//record parameters
extern int HBnumformed;
extern int alphaHBformed;
extern char timer[30];
//-----------------


//-----------------
//multi-thread variables
extern int threadNum;
extern int codeNum;
extern int eventToCommit;
extern signed int readCount, writeCount, activeWrite;
extern struct ThreadStr **thread;
extern struct RandomListStr *threadRandomList;
extern struct ThreadInfoStr *thrInfo;
extern struct PreCalObjStr *preCalList;
extern pthread_t *thread_t;
extern pthread_mutex_t mainFrameLock, rawDataLock, commitLock;
extern pthread_mutex_t mstThrLock, slvThrLock;
extern pthread_cond_t mainFrameRenew, rawDataRenew, readCheck, writeCheck;
extern pthread_cond_t mstThrQ, slvThrQ;


//-----------------
//REMD
extern struct REMDStr REMDInfo;


//-----------------
//Wall
extern struct AtomStr *wall;
extern struct DynamicWallStr wallDyn;
extern struct ObstructionStr obstObj;
extern struct TunnelStr      tunlObj;
extern struct SphObstStr     SphObstObj;
extern char wallExist[20];
extern char wallType[20];


//-----------------
//Flow
extern struct FlowStr flow;


//-----------------
//saving data
extern struct FileStr *fileList;
//-----------------
//To maximum speed, open and close files as less as possible.


//-----------------
//other
extern int visual;
extern int tmpInt;
extern double tmpDouble;
//-----------------




//==============================================
//================FUNCTIONS=====================
//==============================================

//---------------------------------------------
//initialization
//---------------------------------------------
void InputData(int argc, const char * argv[]);
//---------------------------------------------


//---------------------------------------------
//models
//---------------------------------------------
int AAModel(char *);
int AtomModel(char *);
int HBModel(struct AtomStr *, struct AtomStr *);
int NeighborModel(char *, char *);
//---------------------------------------------


//---------------------------------------------
//TimePrediction
//---------------------------------------------
void CollisionTime(struct AtomStr *collision_i, struct ThreadStr *thisThread);
void BondTime(struct AtomStr *bond_i, struct ThreadStr *thisThread);
void HBTime(struct AtomStr *HB_i, struct ThreadStr *thisThread);
void HBNeighborTime(struct AtomStr *neighbor_i, struct ThreadStr *thisThread);
void ThermostatTime(struct AtomStr *targetAtom);
void PBCandCrossCellTime(struct AtomStr *targetAtom);
void WallTime(struct AtomStr *targetAtom);
void ObstTime(struct AtomStr *targetAtom);
void TunnelTime(struct AtomStr *targetAtom);
void ChargeTime(struct AtomStr *targetAtom);
void SphObstTime(struct AtomStr *targetAtom);
//---------------------------------------------


//---------------------------------------------
//ThreadProcess
//---------------------------------------------
int ProcessEvent(int *threadRenewList, struct ThreadStr *thisThread);
int ThreadProcess(void);
int DoEvent(struct ThreadStr *thisThread);
int HazardCheck(struct AtomStr *oldTargetAtom, struct AtomStr *oldPartner, struct AtomStr *oldTargetNeighbor, struct AtomStr *oldPartnerNeighbor, struct ThreadStr *thisThread);
void CommitEvent(struct AtomStr **destLibrary, struct AtomStr *newTargetAtom, struct AtomStr *newPartnerAtom, struct AtomStr *oldTargetAtom, struct AtomStr *oldPartnerAtom, struct AtomStr *oldTargetNeighbor, struct AtomStr *oldPartnerNeighbor);
void Predict(int *renewList, struct ThreadStr *thisThread);
void AssignJob(int *list, struct ThreadStr *thisThread);
void AssignThread(int *renewList, struct ThreadStr *thisThread);
void FirstRun(struct ThreadStr *thisThread);
void DoWallDyn(void);
struct ThreadStr* InitializeThread(int tid, struct AtomStr *atomList);
//---------------------------------------------


//---------------------------------------------
//Simulation Method
//---------------------------------------------
void MSThread(void);
void SingleThread(void);
void REMD(void);
//---------------------------------------------


//---------------------------------------------
//Event
//---------------------------------------------
void LinkList(char * type, struct AtomStr *targetAtom, struct ThreadStr *thisThread);
void PBC(char *type, struct AtomStr *targetAtom, struct ThreadStr *thisThread);
//---------------------------------------------


//---------------------------------------------
//CBT
//---------------------------------------------
void CreateCBT(void);
void UpdateCBT(int *renewList);
void InsertCBT(int *renewList);
void DeleteCBT(int *renewList);
//---------------------------------------------


//--------------------------------------------
//Tool functions
//--------------------------------------------
void SurroundingCheck(int, int);
void Rotation(double *, double *, double, char);
void TimeBenchmark(char *, char *, int);
void PrintCellList(struct ThreadStr *thisThread);
void PrintList(list *atomList);
double TotalPotentialEnergy(int collision_i, double * direction, struct ThreadStr *thisThread);
void SDEnergyMin(long int stepNum, struct ThreadStr *thisThread);
void AtomDataCpy(struct AtomStr * dest, struct AtomStr * sour, int flag);
void RenewCellList(struct AtomStr *newTargetAtom, struct AtomStr *oldTargetAtom);
void ListRefresh(int, int *, int, int);
void CreateGELCoordinate(int);
void ResetTarget(int *renewList, struct ThreadStr *thisThread);
void MakeNeighborList(int *neighborList, int targetAtom, int threadID);
void TimeForward(double time, struct ThreadStr *thisThread);
void UpdateData(double time, char *type, struct ThreadStr *thisThread);
void PBCShift(struct AtomStr *atom_i, struct AtomStr *atom_j, double *shift);
void PrintPreCalList(void);
void PrintBonds(int atomNum);
void PrintConstr(int atomNum);
void PrintStep(struct StepPotenStr *thisStep);
void PrintCollisionPotentialTable(void);
void PrintHBPotentialTable(int HBTypeNum);
void CalAccumulatedPoten(struct StepPotenStr *startPoint);
void FreeConstr(struct ConstraintStr *thisConstr);
void PrintData(int *list);
void CalCOMV(int proteinNum, double *netV);
void ChangeColor(int type, double *color);
int AtomTypeChange(int originalType, int direct);
int FindPair(struct AtomStr *atom1, struct AtomStr *atom2, char *interactionType, double direction, double distance2, double *lowerLimit, double *upperLimit, double *lowerPotential, double *upperPotential, double *accumPotential, struct AtomStr *HBPartner, int typeChange);
//HBPartner is only used during finding pair between HB target atom and its neighbor. HBPartner is the partner atom of the target atom.
//typeChange is only used during calculating the energy change before and after HB forming or breaking
int EnterTunnel(double time, struct AtomStr *targetAtom);
double RandomValue(int, int);
double Maxwell_Boltzmann_Distribution(double, double);
double RandomVelocity(double mass);
double CalKinetE(struct ThreadStr *thisThread);
double CalPotenE(struct ThreadStr *thisThread);
double CalSysTem(struct ThreadStr *thisThread);
double GetViscosity(double);
double AbsvalueVector(double * vector);
double CalDistance(double *p1, double *p2);
double FindPotWellWidth(struct AtomStr *atom_i, struct AtomStr *atom_j);
struct ConstraintStr *RightPair(int type1, int type2, int flag);
//---------------------------------------------


//---------------------------------------------
//Data parameters output
//---------------------------------------------
void DisplayTime(char *);
//---------------------------------------------


//---------------------------------------------
//DataSave
//---------------------------------------------
void SaveData(enum FileType type, struct ThreadStr *thisThread);
void SavePDB(struct FileStr *file);
void SaveGRO(struct FileStr *file);
void GlobalCloseFree(void);
struct FileStr* InitializeFiles(char *extra, struct FileStr* preList);
//---------------------------------------------


//---------------------------------------------
//Visual
//---------------------------------------------
void VisualSGThread(void);
//---------------------------------------------


#endif /* DMD_h */

