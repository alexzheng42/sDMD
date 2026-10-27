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

// size of every buffer that holds a file path (directory + file name)
#define PATH_LEN 4096

//#define DEBUG_FF
//#define DEBUG_PRINTF
//#define GEL
#define HYDROGEN_BOND

#define TIME_GAP      0.5
#define INFINIT       1E8
#define ZERO          1E-8
#define INVALID      -111
#define EXTEND_RATIO  0.90 //make sure the cut off radius is less than the subcell length

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
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

#define _CALLOC2(type, name, size1, size2)    \
name=(type **)calloc(size1, sizeof(type *));   \
name[0]=(type *)calloc(size1*size2, sizeof(type));  \
for (int NN=1; NN<size1; NN++) name[NN]=name[0]+NN*size2;

#define _FREE2(name) \
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

#define NEIGHBORSCAN(list, count) \
for (int count = 1; count <= list[0]; count ++)

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

enum EChargeRes {
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

// HB classes: energies of each HB type are scaled per class (PRO_RNA_MOL.DMDParam
// section D) and a class barrier is subtracted on formation
enum EHBClass { HB_BB, HB_BS, HB_SS, HB_CLASSNUM };

struct HBPotentialStr
{
    int typeNum;
    int *classOf;                 // class of each HB type (index 1..typeNum)
    double scale[HB_CLASSNUM];
    double barrier[HB_CLASSNUM];
};

struct AtomSeqStr {
    int atomNum;
    int molNum;
    int resNum;
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
    double dmax;
    double e;
    double accumulated;
    struct StepPotenStr *step;
    struct ConstraintStr *next;
};


// An atom holds up to maxPartner HBs at a time (1 unless its residue file raises it
// with HBMAX, e.g. a water O that accepts up to three HBs); maxPartner = 0 means 1.
#define HB_MAX_PARTNER 3

// savedData.dat starts with this marker and a format version. Files written before the
// multi-partner HB layout start directly with the atom number and are read with the legacy layout.
#define SAVED_DATA_MAGIC   (-1)
#define SAVED_DATA_VERSION 3 // 3: times relative to the time origin, timeBaseCount at the end

// SysInfo.dat (read by sAnalysis) starts with this marker and a format version; version 2 writes
// every list with an explicit count. Older files start directly with the atom number.
#define SYSINFO_MAGIC   (-1)
#define SYSINFO_VERSION 2

struct HBStr {
    int bondConnection[HB_MAX_PARTNER + 1]; //[0] records the num of HB partners, [1..] the partners
    int maxPartner;
    int neighbor[3]; //at most two neighbors, [0] records the num
    char role;
    enum EHBEventType interactionType;
};

// HB-neighbor partners per atom: a water O with three HBs accepted (two donor neighbors each)
// and two donated has 8
#define HB_NEIGHBOR_MAX 8

struct HBNeighborStr{
    int neighborPartner[HB_NEIGHBOR_MAX + 2]; //[1..HB_NEIGHBOR_MAX] used; the last stays 0 so unbounded list scans stop inside the array
    int partnerNum;
    int neighborStatus;
    int removedAtom[5];
};

struct InteractionEventStr {
    int flag;
    int dim;
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
    char extraProperty[10][16];
    char resName[8];
    char typeName[3][16]; //[0]: NULL, [1]: type1, [2]: type2
    int num;
    int typeofAtom[3]; //[0]: index, [1]: original, [2]: reacted
	int typeofRes;
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
    double lastUpdateTime; // Add this field
    struct HBStr HB;
    struct HBNeighborStr HBNeighbor;
    struct InteractionEventStr event;
};

// Per-atom dynamic data in the layout used before multi-partner HBs (one HB partner, 5-entry
// HB-neighbour list). sAnalysis reads SysInfo.dat, WallInfo.dat and ObstInfo.dat in this layout,
// and savedData.dat files from before the change use it.
struct HBStrLegacy {
    int bondConnection;
    int neighbor[3];
    char role;
    enum EHBEventType interactionType;
};

struct HBNeighborStrLegacy {
    int neighborPartner[5];
    int partnerNum;
    int neighborStatus;
    int removedAtom[5];
};

struct DynamicStrLegacy {
    int cellIndex[4];
    double coordinate[4];
    double velocity[4];
    double lastUpdateTime;
    struct HBStrLegacy HB;
    struct HBNeighborStrLegacy HBNeighbor;
    struct InteractionEventStr event;
};

void DynamicToLegacy(const struct DynamicStr *dynamic, struct DynamicStrLegacy *old);

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

struct ResidueStr {
    char name[8];
    int type;
    int startAtomNum;
    int endAtomNum;
    int molNum;
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
    struct CGBeadStr *CGBead;
};

struct MolStr {
    int startAtomNum;
    int endAtomNum;
    int startResNum;
    int endResNum;
    double mass;
};

// neighbor-cell stencil: cells of side >= cutoffr / k are scanned out to k cells
// away, keeping only cells whose nearest point is within cutoffr (k <= 3)
#define MAX_STENCIL 343

struct AtomListStr {
    int targetNum;
    int stencilNum;
    int x[MAX_STENCIL], y[MAX_STENCIL], z[MAX_STENCIL];
    double timeInr;
    struct AtomStr *target;
    struct AtomStr *partner;
    struct AtomStr **ptr;
};

struct FileStr {
    int mark;
    char name[PATH_LEN];
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


struct RegionStr {
    double xyz[4][2];
    double r;
};


struct ObstructionStr {
    int mark;
    int num;
    int *hole;
    double **position;
    struct AtomStr *obst;
    struct RegionStr *region;
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


struct FixedAtomStr {
    int mark;
    int num;
    int *list;
};


struct SurfaceStr {
    int mark;
    char type[32];
    char style[32];
};


struct FFBeadStr {
    char name[16];
    int val[4];
};


struct FFBondStr {
    char name[2][16];
    double kb, b0;
    double dmin, dmax, delta;
};


struct FFAngleStr {
    char name[2][16];
    double ktheta, theta0;
    double kub, r0;
    double minDeg, maxDeg;
    double dmin, dmax, delta;
    struct FFBeadStr interAtom;
};


struct DisPotStr {
    double val, dis, pot;
};

struct DihedralDataStr {
    int step;
    double bond[3], angle[2];
    char dataString[512];
};


struct FFDihedralStr {
    char name[2][16];
    double kphi, n, delta;
    struct FFBeadStr interAtom[2];
    struct DihedralDataStr data;
    struct DisPotStr *list;
};

struct FFImproperStr {
    char name[2][16];
    double komega, omega0;
    struct FFBeadStr interAtom[2];
    struct DihedralDataStr data;
    struct DisPotStr *list;
};

struct FFLJStr {
    char name[16];
    double mass, charge, sigma, epsilon;
};


struct FFInfoStr {
    struct FFBeadStr        *group;
    struct FFBeadStr        *bead;
    struct FFBondStr        **bond;
    struct FFAngleStr       **angle;
    struct FFDihedralStr    **dihedral;
    struct FFImproperStr    **improper;
    struct FFLJStr          *LJ;
};


//--------------------------------------------
//Global Variables
//--------------------------------------------


//-----------------
//directory
extern char datadir[PATH_LEN];
extern char savedir[PATH_LEN];
//-----------------


//-----------------
//basic variables
extern int atomnum;
extern int **connectionMap;
extern int chargeRes[9];
extern long int frame;
extern unsigned int seed;
extern double timestep, currenttime;
// Time origin: all stored times are relative to timeBaseCount * TIME_BASE_UNIT (see RebaseTime)
#define TIME_BASE_UNIT 1024.0
extern long long timeBaseCount;
#define REAL_TIME(t) ((t) + timeBaseCount * TIME_BASE_UNIT)
extern double oldcurrenttime;
extern double cutoffr;
extern double instTemperature, targetTemperature, oldTemperature;
extern double TotlKinEner;
extern double viscosity; //dynamic u/nm/ps
extern double outputrate;
extern double outputrecord;
extern double boxDimension[4]; //unit A
extern struct AtomStr *atom;
extern struct ResidueStr *residue;
extern struct MolStr *molecule;
extern struct CGGovStr CG;
extern struct CGMatrixStr CGMatrix;
extern struct PreBondStr preBond;
extern struct FixedAtomStr fixedAtomInfo;
extern char neworcontinue[32];
extern char Methodtype[32];
extern char FFLibrary[32];
extern char temperatureType[32];
//-----------------


//-----------------
//force field info
extern struct FFInfoStr FFInfo;
//-----------------


//-----------------
//system check
extern long int countReCal;
extern long int countBackup;
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
//molecule sequence
extern int numofMol;
extern int numofRes;
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
extern struct ConstraintStr **potentialPairCollision;
extern struct ConstraintStr **potentialPairCG;
extern struct ConstraintStr ***potentialPairHB;
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
extern int cellFactor;


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


//-----------------
//Surface
extern struct SurfaceStr surfInfo;


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
void InitializeCG(void);
void InitializeCGPotential(void);
//---------------------------------------------


//---------------------------------------------
//FF
//---------------------------------------------
void ReadFFInfo(char *method);
void BondDuplicate(int ResNum);
void ConstrDuplicate(int ResNum);
void RemoveConstr(struct PropertyStr *property, int atomNum, char *type);
void CreateLJData(FILE *outputFile);
void BuildSurfaceConnection(void);

int ReadAtomInfo(char *info, int resNum);
int ReadBondInfo(char *info, int resNum);
int ReadAngleInfo(char *info, int resNum);
int ReadDihedralInfo(char *infoType, char *info, int resNum);
int ReadHBInfo(char *infoType, char *info, int resNum);
int FindResIndex(char *name);
int FindAtomIndex(char *);
int FindTargetAtom(int ResNum, char *target);
int HBModel(struct AtomStr *, struct AtomStr *);
double HBBarrier(struct AtomStr *atom1, struct AtomStr *atom2);
int NeighborModel(char *, char *);
int AAList(char *name);
int SingleParticleList(char *name);
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
void InitBackupEvents(void);
void ShiftBackupEvents(double shift);
void ClearBackupEvents(int atomNum);
int TakeBackupEvent(struct AtomStr *target, struct AtomListStr *thisList);
//---------------------------------------------


//---------------------------------------------
//SimulationProcess
//---------------------------------------------
void Prepare(void);
void FreeVariables(void);
void Predict(int* renewList, struct AtomListStr* thisList);
void RebaseTime(void);
void AssignJob(struct AtomStr** target, struct AtomStr** partner, int* renewList, struct AtomListStr* thisList);
void DoEvent(struct AtomListStr *thisList);
void DoWallDyn(void);
void FixAssignRenewList(int *renewList, struct AtomListStr *thisList);
int HazardCheck(struct AtomStr *target, struct AtomStr *partner, struct AtomListStr *thisList);
//---------------------------------------------


//---------------------------------------------
//Simulation Method
//---------------------------------------------
void SingleThread(void);
void REMD(void);
void SingleThreadRun(void (*afterEvent)(void), void (*beforeFree)(void));
//---------------------------------------------


//---------------------------------------------
//Event
//---------------------------------------------
void LinkList(char * type, struct AtomStr *target, int *oldIdx, struct AtomListStr *thisList);
void PBC(char *type, struct AtomStr *target, struct AtomListStr *thisList);
void HBAddPartner(struct AtomStr *HB_i, int partner);
void HBRemovePartner(struct AtomStr *HB_i, int partner);
void AddNeighbor2HB(int *currentStatus, int *list, int *neighborList);
void AddHB2Neighbor(int *list, int HBNum);
void SwitchConnection(int atom_i, int atom_j, int *neighborList_i, int *neighborList_j);
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
void UpdateAtom(struct AtomStr *target, double targetTime);
void UpdateAllAtoms(double targetTime);
void CellPBCShift(int *cell_i, int *cell_j, double *shift);
void BuildStencil(struct AtomListStr *thisList);
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
void CalCOMV(int molNum, double *netV);
void ChangeColor(int type, double *color);
enum EInterPairType {
    InterColli,
    InterConstr,
    InterHB,
    InterNeighbor,
    InterCG,
    InterNone
};

int FindPair(struct AtomStr *atom1, struct AtomStr *atom2, enum EInterPairType interactionType, double direction, double distance2, double *lowerLimit, double *upperLimit, double *lowerPotential, double *upperPotential, double *accumPotential, struct AtomStr *HBPartner, int typeChange);
//HBPartner is only used during finding pair between HB target atom and its neighbor. HBPartner is the partner atom of the target atom.
//typeChange is only used during calculating the energy change before and after HB forming or breaking
int EnterTunnel(double time, struct AtomStr *target);
int FindElemInList(int *thisList, int value, int start, int end, int add);
// inline: called for every neighbor in the event-prediction scans
static inline int HBCapacity(const struct AtomStr *a) {
    return a->dynamic->HB.maxPartner > 0 ? a->dynamic->HB.maxPartner : 1;
}

static inline int HBPartnerNum(const struct AtomStr *a) {
    return a->dynamic->HB.bondConnection[0];
}

// TRUE if a and atom #b are HB partners
static inline int HBIsPartner(const struct AtomStr *a, int b) {
    for (int i = 1; i <= a->dynamic->HB.bondConnection[0]; i++)
        if (a->dynamic->HB.bondConnection[i] == b) return TRUE;
    return FALSE;
}

static inline int HBHasFreeSlot(const struct AtomStr *a) {
    return a->dynamic->HB.bondConnection[0] < HBCapacity(a);
}

// Neighbor pairs of an HB: host (an HB atom) - an HB.neighbor atom of one of host's
// partners. Returns that partner, or 0 if other is not such a neighbor of host.
static inline int HBNeighborHostPartner(const struct AtomStr *host, int other) {
    for (int p = 1; p <= host->dynamic->HB.bondConnection[0]; p++) {
        const int *list = atom[host->dynamic->HB.bondConnection[p]].dynamic->HB.neighbor;
        for (int n = 1; n <= list[0]; n++)
            if (list[n] == other) return host->dynamic->HB.bondConnection[p];
    }
    return 0;
}

static inline int FindStrInList(char thisList[10][16], char *str, int add) {
    for (int i = 1; i < 9; i++) {
        // compare the first char before strcmp: this is called for every neighbor scan
        if (thisList[i][0] == str[0] && strcmp(thisList[i], str) == 0) {
            return i;
        } else if (thisList[i][0] == '\0') {
            if (add)
                sprintf(thisList[i], "%s", str);
            return 0;
        }
    }
    return 0;
}
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
struct ConstraintStr *RightPair(int *atom1_typeList, int type1, int *atom2_typeList, int type2, int flag);
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

