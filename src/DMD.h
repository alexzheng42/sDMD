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

//#define DEBUG_RANDOM
//#define DEBUG_PRINTF
//#define DETAILS
//#define GEL
#define HYDROGEN_BOND

#define TIME_GAP      0.5
#define INFINIT       1E8
#define ZERO          1E-8
#define INVALID      -111
#define NATOMTYPE     32
#define EXTEND_RATIO  0.90 //make sure the cut off radius is less than the subcell length

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifdef DEBUG_RANDOM
#define RANDOM_SEED 1575263273
#endif

#define FALSE 0
#define TRUE 1

/*
 Units of the system
 Temperature: kcal / (mol * kB) around 503.2 K (kB = 0.001987204118 kcal / (mol * K))
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

#define INT_2CALLOC(name, size1, size2)    \
name=(int **)calloc(size1, sizeof(int *));   \
name[0]=(int *)calloc(size1*size2, sizeof(int));  \
for (int NN=1; NN<size1; NN++) name[NN]=name[0]+NN*size2;

#define DOUBLE_CALLOC(name, size) \
name=(double *)calloc(size, sizeof(double));

#define DOUBLE_2CALLOC(name, size1, size2)    \
name=(double **)calloc(size1, sizeof(double *));   \
name[0]=(double *)calloc(size1*size2, sizeof(double));  \
for (int NN=1; NN<size1; NN++) name[NN]=name[0]+NN*size2;

#define FREE2(name) \
free(name[0]);  \
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
    CGSu_Event,  //CG surface event
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
    wallInfo,    //saved data for wall
    obstInfo,    //saved data for obstruct
    CGInfo,      //saved data for CG
    savedData,   //saved data for continuity
    RE,          //replica exchange data
    lenFileType  //length of this list
};

struct HBPotentialStr //B: Backbone, S: Sidechain
{
    double BB_r, BB_v;
    double SS_r, SS_v;
    double BS_r, BS_v;
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
    char name[8];
    char extraProperty[2][16];
    char nameofAA[8];
    int num;
    int typeofAtom;
	int typeofAA;
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
    char nameofAA[5];
    int type;
    int startAtomNum;
    int endAtomNum;
    int proteinNum;
    double mass;
};

struct CGStr {
    int start;
    int end;
    int type;
};

struct CGGovStr {
    char type[2][16];
    int mark;
    int totalNum;
    int *list;
    struct CGStr *thisCG;
};

struct CGBeadStr {
    double chi;
    double epsilon;
    double sigma;
    double delta;
};

struct CGMatrixStr {
    double theta[5]; // 1 2 3 s p
    struct CGBeadStr CGBead[NATOMTYPE];
};

struct PepStr {
    int startAtomNum;
    int endAtomNum;
    int startAANum;
    int endAANum;
    double mass;
};

struct AtomListStr {
    int targetNum;
    int x[27], y[27], z[27];
    double timeInr;
    struct AtomStr *target;
    struct AtomStr *partner;
    struct AtomStr **ptr;
};

struct FileStr {
    int mark;
    char name[256];
    FILE *file;
};


struct REMDTempStr {
    int num;
    double T;
};


struct REMDStr {
    int flag;
    int REMD_PortNum;
    int REMD_OutputRate;
    char REMD_ServerName[50];
    char REMD_ExtraName[50];
    struct REMDTempStr REMD_Temperature;
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


struct PreBondStr {
    int mark;
    int num;
    double dis;
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
extern struct CGGovStr CG;
extern struct CGMatrixStr CGMatrix;
extern struct PreBondStr preBond;
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
extern long int CGeventsum;
extern long int oldtotaleventsum;
extern long int newtotaleventsum;
extern long int warningsum;
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
extern struct ConstraintStr potentialPairCollision[NATOMTYPE + 1][NATOMTYPE + 1];
extern struct ConstraintStr potentialPairHB[12][NATOMTYPE + 1][NATOMTYPE + 1];
extern struct ConstraintStr potentialPairCG[NATOMTYPE + 1][NATOMTYPE + 1];
extern struct HBPotentialStr HBPotential;


//-----------------
//thermostat and solvent
extern char thermostatType[20];
extern double thermoF;
//-----------------


//-----------------
//record parameters
extern int HBnumformed;
extern int alphaHBformed;
extern char timer[30];
//-----------------


//-----------------
//multi-thread variables


//-----------------
//single-thread variable
extern int renewList[128];
extern struct AtomListStr atomList;


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
void InitializeCGPotentialMatrix(void);
//---------------------------------------------


//---------------------------------------------
//models
//---------------------------------------------
int FindNumInAA(char *);
int AAModel(char *);
int AtomModel(char *);
int HBModel(struct AtomStr *, struct AtomStr *);
int NeighborModel(char *, char *);
double HBBarrier(struct AtomStr *atom1, struct AtomStr *atom2);
//---------------------------------------------


//---------------------------------------------
//TimePrediction
//---------------------------------------------
void CollisionTime(struct AtomStr *collision_i, struct AtomListStr *thisList);
void BondTime(struct AtomStr *bond_i, struct AtomListStr *thisList);
void HBTime(struct AtomStr *HB_i, struct AtomListStr *thisList);
void HBNeighborTime(struct AtomStr *neighbor_i, struct AtomListStr *thisList);
void CGSurfaceTime(struct AtomStr *target);
void ThermostatTime(struct AtomStr *target);
void PBCandCrossCellTime(struct AtomStr *target);
void WallTime(struct AtomStr *target);
void ObstTime(struct AtomStr *target);
void TunnelTime(struct AtomStr *target);
void ChargeTime(struct AtomStr *target);
void SphObstTime(struct AtomStr *target);
int JobAssign(struct AtomStr* atom_i, struct AtomStr* atom_j, struct InteractionEventStr* event, enum EEventType type);
//---------------------------------------------


//---------------------------------------------
//SimulationProcess
//---------------------------------------------
void Prepare(void);
void FreeVariables(void);
void Predict(int* renewList, struct AtomListStr* thisList);
void AssignJob(struct AtomStr** target, struct AtomStr** partner, int* renewList, struct AtomListStr* thisList);
void DoEvent(struct AtomListStr *thisList);
void DoWallDyn(void);
void FixAssignRenewList(int *renewList, struct AtomListStr *thisList);
int HazardCheck(struct AtomStr *target, struct AtomStr *partner, struct AtomStr *targetNeighbor, struct AtomStr *partnerNeighbor, struct AtomListStr *thisList);
//---------------------------------------------


//---------------------------------------------
//Simulation Method
//---------------------------------------------
void SingleThread(void);
void REMD(void);
//---------------------------------------------


//---------------------------------------------
//Event
//---------------------------------------------
void LinkList(char * type, struct AtomStr *target, int *oldIdx, struct AtomListStr *thisList);
void PBC(char *type, struct AtomStr *target, struct AtomListStr *thisList);
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
void ThreadProcess(void);
void SDEnergyMin(long int stepNum, struct AtomListStr *thisList);
void AtomDataCpy(struct AtomStr * dest, struct AtomStr * sour, int flag);
void ListRefresh(int, int *, int, int);
void CreateGELCoordinate(int);
void ResetTarget(int *renewList, struct AtomListStr *thisList);
void TimeForward(double time);
void UpdateData(double time);
void PBCShift(struct AtomStr *atom_i, struct AtomStr *atom_j, double *shift);
void PrintBonds(int atomNum);
void PrintConstr(int atomNum);
void PrintStep(struct StepPotenStr *thisStep);
void PrintCollisionPotentialTable(void);
void PrintHBPotentialTable(int HBTypeNum);
void CalAccumulatedPoten(struct StepPotenStr *startPoint);
void FreeConstr(struct ConstraintStr *thisConstr);
void PrintListData(int *thisList);
void PrintHBData(void);
void CalCOMV(int proteinNum, double *netV);
void ChangeColor(int type, double *color);
int AtomTypeChange(int originalType, int direct);
int FindPair(struct AtomStr *atom1, struct AtomStr *atom2, char *interactionType, double direction, double distance2, double *lowerLimit, double *upperLimit, double *lowerPotential, double *upperPotential, double *accumPotential, struct AtomStr *HBPartner, int typeChange);
//HBPartner is only used during finding pair between HB target atom and its neighbor. HBPartner is the partner atom of the target atom.
//typeChange is only used during calculating the energy change before and after HB forming or breaking
int EnterTunnel(double time, struct AtomStr *target);
int FindElemInList(int *thisList, int value, int start, int end, int add);
double RandomValue(int, int);
double Maxwell_Boltzmann_Distribution(double, double);
double RandomVelocity(double mass);
double CalKinetE(void);
double CalSysTem(void);
double CalPotenE(void);
double CalExponential(double value, int times);
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
void SaveData(enum FileType type, struct FileStr *fileList);
void SavePDB(struct FileStr *file);
void SaveGRO(struct FileStr *file);
void GlobalCloseFree(void);
void InitializeFiles(char *extra, struct FileStr* fileList);
//---------------------------------------------


//---------------------------------------------
//Visual
//---------------------------------------------
void VisualSGThread(void);
//---------------------------------------------


#endif /* DMD_h */

