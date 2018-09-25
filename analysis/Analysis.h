//
//  Analysis.h
//  Analysis
//
//  Created by Size Zheng on 10/30/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define TRUE      1
#define FALSE     0
#define INVALID  -111
#define BOLTZMANN 0.0019872041 //kcal / mol / K

#define bondConnect             0b00001
#define constraintConnect       0b00010
#define HBConnect               0b00100
#define neighborConnect         0b01000
#define overlapConnect          0b10000

#define findduplicatedandadd(list1, list2, num, flag)  \
flag=0; \
num=1;    \
while (list2[num]!=0) {   \
if (list1==list2[num]) {   \
flag=1; \
break;  \
} else { \
num++; \
}   \
}   \
if (flag==0) {  \
list2[num]=list1; \
}

#define min3(num1, num2, num3)    \
(((num1) < (num2)) ? (((num1) < (num3)) ? (num1) : (num3)) :    \
(((num2) < (num3)) ? (num2) : (num3)))

#define max3(num1, num2, num3)    \
(((num1) > (num2)) ? (((num1) > (num3)) ? (num1) : (num3)) :    \
(((num2) > (num3)) ? (num2) : (num3)))

#define transfer_vector(recv, send)   \
recv[1]=send[1];    \
recv[2]=send[2];    \
recv[3]=send[3];

#define dotminus(n1, n2, out) \
out[1]=n1[1]-n2[1]; \
out[2]=n1[2]-n2[2]; \
out[3]=n1[3]-n2[3];

#define dotplus(n1, n2, out) \
out[1]=n1[1]+n2[1]; \
out[2]=n1[2]+n2[2]; \
out[3]=n1[3]+n2[3];

#define crossprod(n1, n2, result)   \
result[1]=n1[2]*n2[3]-n1[3]*n2[2];  \
result[2]=n1[3]*n2[1]-n1[1]*n2[3];  \
result[3]=n1[1]*n2[2]-n1[2]*n2[1];

#define factorprod(factor, vector, out) \
out[1]=factor*vector[1];    \
out[2]=factor*vector[2];    \
out[3]=factor*vector[3];

#define dotprod(n1, n2) \
n1[1]*n2[1]+n1[2]*n2[2]+n1[3]*n2[3];

#define pointToStruct(pointer, struct)  \
pointer=&struct[0];

#define DOUBLE_2CALLOC(name, size1, size2)    \
name=(double **)calloc(size1, sizeof(double *));   \
for (int NN=0; NN<size1; NN++) {    \
name[NN]=(double *)calloc(size2, sizeof(double)); \
}

enum EOperation {
    empty,
    RemovePBC,
    AnalyzeEn,
    AnalyzeHB,
    AnalyzeAgg,
    AnalyzeRamach,
    AnalyzeConMap,
    CalculateRG,
    CalculateMSD,
    NumofElem
};

enum EFileType {
    inSys,
    inTrj,
    inCnt,
    inLog,
    inREMDTemp,
    outTrj,
    outPE,
    outKE,
    outWE,
    outTolE,
    outHBInfo,
    outAAMark,
    outTemp,
    outAgg,
    NumofFileType
};

enum EEventType {
    Coli_Event,  //collision event
    Bond_Event,  //bond event
    HB_Event,    //HB event
    HBNe_Event,  //HB neighbor event
    Ther_Event,  //thermostat event
    Wall_Event,  //wall event
    PCC_Event,   //PBC and crossing event
    Lagv_Event,  //Langevin event
    TBD_Event,   //TBD event
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
    Invalid
};

enum EHBEventType {
    NoEvent,
    HBForm,
    HBBreak,
    HBSkip
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

struct HBPotentialStr //B: Backbone, S: Sidechain
{
    double BB;
    double SS;
    double BS;
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

struct HBType {
    int helix_alpha;
    int helix_310;
    int helix_pi;
    int beta;
    int other;
    int total;
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


struct TunnelStr {
    int mark;
    double startPosition;
    double endPosition;
    double diameter;
    struct AtomStr tunnel;
};


struct WallStr {
    char wallExist[20];
    char wallType[20];
    struct AtomStr wall;
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


struct REMDStr {
    int mark;
    int numReplica;
};

struct FileStr {
    enum EFileType fileType;
    char name[256];
};

struct FileListStr {
    int count;
    char list[128][64];
};


struct SectionStr {
    long frameCount;
    double oldTime;
    double targetTime;
    double temp;
    double outputRate;
    char methodType[20];
    struct WallStr wallObj;
    struct DynamicWallStr wallDyn;
    struct ObstructionStr obstObj;
    struct TunnelStr tunlObj;
    struct FlowStr flow;
};


extern int atomnum;
extern int numofprotein;
extern int totalAminoAcids;
extern int **connectionMap;
extern int *celllist;
extern int freshStart;
extern int *analysisList;
extern int nPP;

extern double cutoffr;
extern double boxOrigDim[4];
extern double boxCurtDim[4];
extern double cellNum[4];
extern double cellSize[4];

extern char path[1024];
extern char names[NumofFileType][256];
extern char obstDir[1024];
extern char targetPeptideNum[16];

extern struct FileStr **files;
extern struct AAStr *aminoacid;
extern struct PepStr *protein;
extern struct AtomStr *atom;
extern struct HBType HBSum;
extern struct ConstraintStr potentialPairCollision[32][32];
extern struct ConstraintStr potentialPairHB[11][32][32];
extern struct HBPotentialStr HBPotential;
extern struct REMDStr RE;
extern struct FileListStr fileList;
extern struct SectionStr *sectInfo;

void InitializeFiles(int row, int column);
void SystemInformationInput(int id);
void AdjustPBC(int id);
void FreeVariables(void);
void HBInfo(int id);
void ClusterInfo(int id);
void EnergyInfo(int id);
void EstCell(void);
void LinkList(void);
void PrintProcess(long step);
void AssignFileList(int id);
void FindTargetFile(char *oldName, char *fileListName, char *newName);
int ReadLog(void);
int ReadGro(FILE *inputFile);
int ReadConnectionMap(FILE *inputFile);
int WHAM(int argc, const char * argv[]);
int CheckHBConnection(int thisAtomNum);
int AtomModel(char *type);
double absvalue_vector(double * vector);
