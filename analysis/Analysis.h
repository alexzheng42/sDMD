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
#define INFINIT   1E8
#define BOLTZMANN 0.0019872041 //kcal / mol / K

#define bondConnect             0b00001
#define constraintConnect       0b00010
#define HBConnect               0b00100
#define neighborConnect         0b01000
#define overlapConnect          0b10000

#define min3(num1, num2, num3)    \
(((num1) < (num2)) ? (((num1) < (num3)) ? (num1) : (num3)) :    \
(((num2) < (num3)) ? (num2) : (num3)))

#define max3(num1, num2, num3)    \
(((num1) > (num2)) ? (((num1) > (num3)) ? (num1) : (num3)) :    \
(((num2) > (num3)) ? (num2) : (num3)))

#define TRANSFER_VECTOR(recv, send)   \
recv[1]=send[1];    \
recv[2]=send[2];    \
recv[3]=send[3];

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

#define POINT_TO_STRUCT(pointer, struct)  \
pointer=&struct[0];

#define _CALLOC2(type, name, size1, size2)    \
name=(type **)calloc(size1, sizeof(type *));   \
name[0]=(type *)calloc(size1*size2, sizeof(type));  \
for (int NN=1; NN<size1; NN++) name[NN]=name[0]+NN*size2;

#define _FREE2(name)    \
free(name[0]);  \
free(name); \
name = NULL;



enum EOperation {
    empty,
    RemovePBC,
    AnalyzeEn,
    AnalyzeHB,
    AnalyzeAgg,
    AnalyzeRamach,
    AnalyzeConMap,
    CalculateRG,
    CalculateRG2Obst,
    CalculateDensity,
    CalculatePSD,
    CalculateMSD,
    CalculateRMSD,
    CalculatePotMap,
    NumofElem
};

enum EFileType {
    inSys,
    inWall,
    inObst,
    inCG,
    inTrj,
    inCnt,
    inLog,
    inREMDTemp,
    outTrj,
    outEne,
    outHBInfo,
    outAAMark,
    outAgg,
    outRG,
    outRMSD,
    outPotMap,
    outCG2Obst,
    outDensity,
    outPSD,
    outREMD_T,
    outPotMap_T,
    outPMFMap_T,
    NumofFileType
};

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
    double e, accumulated;
    struct StepPotenStr *step;
    struct ConstraintStr *next;
};

#define HB_CLASSNUM 3   // backbone-backbone, backbone-sidechain, sidechain-sidechain (as in sDMD)

// SysInfo.dat format marker / version written by sDMD (DMD.h)
#define SYSINFO_MAGIC   (-1)
#define SYSINFO_VERSION 2

struct HBPotentialStr
{
    int typeNum;
    int *classOf;                 // class of each HB type (index 1..typeNum)
    double scale[HB_CLASSNUM];
    double barrier[HB_CLASSNUM];
};

struct HBStr {
    int bondConnection;
    int neighbor[3]; //at most two neighbors, [0] records the num
    char role;
    enum EHBEventType interactionType;
};

struct HBNeighborStr{
    int neighborPartner[5];
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

// same layout as sDMD's DynamicStrLegacy (the layout of WallInfo.dat / ObstInfo.dat)
struct DynamicStr {
    int cellIndex[4];
    double coordinate[4];
    double velocity[4];
    double lastUpdateTime;
    struct HBStr HB;
    struct HBNeighborStr HBNeighbor;
    struct InteractionEventStr event;
};

struct AtomStr {
    struct PropertyStr *property;   //static properties of the atom, such as mass, atom # ...
    struct DynamicStr *dynamic;    //dynamic properties of the atom, such as coordinate, velocity ...
};

struct ResidueStr {
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

struct DisStr {
    char type[256];
    int mark;
};


struct DensityStr {
    int mark;
    int dim;
    int pointNum;
};


struct PSDStr {
    int mark;
    int axis[4];
};


struct MolStr {
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
    int SS;
    int other;
    int total;
    int interBB;
    int interBS;
    int interSS;
    int intraBB;
    int intraBS;
    int intraSS;
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
    double kphi, n, phi0;
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
    struct FFLJStr          *LJ;
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
    char methodType[20];
    long frameCount;
    double oldTime;
    double targetTime;
    double temp;
    double outputRate;
    struct WallStr wallObj;
    struct DynamicWallStr wallDyn;
    struct ObstructionStr obstObj;
    struct TunnelStr tunlObj;
    struct FlowStr flow;
};

struct EnergyReadStr {
    double step;
    double T;
    double KE;
    double PE;
    double WE;
    double TE;
};

struct HBReadStr {
    double step;
    struct HBType type;
};

struct RMSDReadStr {
    double step;
    double vRMSD[2];
};

struct GroupStr {
    int molNum;
    int *list_Mol;
    struct MolStr **molList;
    int atomNum;
    int *list_Atom;
    struct AtomStr **atomList;
};


extern int atomnum;
extern int numofMol;
extern int numofRes;
extern int **connectionMap;
extern int *celllist;
extern int freshStart;
extern int *analysisList;

extern double cutoffr;
extern double boxOrigDim[4];
extern double boxCurtDim[4];
extern double cellNum[4];
extern double cellSize[4];

extern char path[1024];
extern char names[NumofFileType][256];

extern struct FileStr **files;
extern struct ResidueStr *residue;
extern struct MolStr *molecule;
extern struct CGGovStr CG;
extern struct DisStr disInfo;
extern struct DensityStr densInfo;
extern struct PSDStr psdInfo;
extern struct AtomStr *atom;
extern struct HBType HBSum;
extern struct ConstraintStr **potentialPairCollision;
extern struct ConstraintStr **potentialPairCG;
extern struct ConstraintStr ***potentialPairHB;
extern struct HBPotentialStr HBPotential;
extern struct REMDStr RE;
extern struct FileListStr *fileList;
extern struct FileListStr RMSDFile;
extern struct SectionStr *sectInfo;
extern struct EnergyReadStr energy;
extern struct HBReadStr hb;
extern struct RMSDReadStr rmsd;
extern struct FFInfoStr FFInfo;
extern struct GroupStr group;

void InitializeFiles(int row, int column);
void SystemInformationInput(int id);
void ReadObst(int id, int sectNum);
void ReadWall(int id, int sectNum);
void ReadCGPot(int id);
void AdjustPBC(int id);
void FreeVariables(void);
void HBInfo(int id);
void ClusterInfo(int id);
void EnergyInfo(int id);
void RGInfo(int id);
void DisInfo(int id);
void DensityInfo(int id);
void PSDInfo(int id);
void RMSDInfo(int id);
void PESurfaceInfo(int id);
void PESurfaceInfoPerT(void);
void EstCell(void);
void LinkList(void);
void PrintProcess(long step);
void AssignFileList(int id);
void FindTargetFile(char *oldName, char *fileListName, char *newName);
void REMDPESurfaceCombine(void);
void PrintStep(struct StepPotenStr *thisStep);
void CreateGroup(void);
int ReadLog(int id);
int ReadGro(FILE *inputFile);
int ReadConnectionMap(FILE *inputFile);
int WHAM(int argc, const char * argv[]);
int CheckHBConnection(int thisAtomNum);
int AtomModel(char *type);
int AAModel(char *type);
int HBModel(struct AtomStr *atom1, struct AtomStr *atom2);
int HBNeighborHost(int host, int other);
int ReadEnergyFile(FILE *EnergyInputFile, struct EnergyReadStr *thisE);
int ReadHBFile(FILE *HBInputFile, struct HBReadStr *thisHB);
int ReadRMSDFile(FILE *RMSDInputFile, struct RMSDReadStr *thisRMSD);
int ReadREMDTempFile(FILE *REMDTempFile);
int FindElemInList(int *list, int value, int start, int end, int add);
int FindStrInList(char thisList[10][16], char *str, int add);
double absvalue_vector(double * vector);
