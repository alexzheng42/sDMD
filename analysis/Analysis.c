//
//  main.c
//  HBAnalysis
//
//  Created by Size Zheng on 4/14/15.
//  Copyright (c) 2015 Size Zheng. All rights reserved.
//

#include "Analysis.h"

int atomnum;
int numofprotein;
int totalAminoAcids;
int *celllist = NULL;
int freshStart = 0;
int *analysisList;
int nPP = 0;

double cutoffr;
double boxOrigDim[4] = {0};
double boxCurtDim[4] = {0};
double cellNum[4] = {0};
double cellSize[4] = {0};

char path[1024];
char names[NumofFileType][256] = {"SysInfo.dat", "", "", "out_log.txt", "out_REMD.txt", "rPBCGRO.gro", "Energy.txt", "HBInfo.txt", "AAMark.txt", "AggNum.txt", "RG.txt", "RMSD.txt", "PotMap.txt", "PotMap_T.txt", "PMFMap_T.txt"};
char obstDir[1024] = "NUll";
char targetPeptideNum[16] = "";

struct FileStr **files;
struct AAStr *aminoacid;
struct PepStr *protein;
struct AtomStr *atom;
struct ConstraintStr potentialPairCollision[NATOMTYPE + 1][NATOMTYPE + 1];
struct ConstraintStr potentialPairHB[12][NATOMTYPE + 1][NATOMTYPE + 1];
struct HBPotentialStr HBPotential;
struct REMDStr RE = {.mark = 0, .numReplica = -1};
struct FileListStr fileList = {.count = 0};
struct FileListStr RMSDFile = {.count = 0};
struct SectionStr *sectInfo;


static void Initialization(int argc, const char *argv[]);
static void DoAnalysis(int argc, const char *argv[]);


int main(int argc, const char * argv[]) {
    
    Initialization(argc, argv);
    
    printf("\n==============================================================\n");
    printf("Analyzing starts:\n");
    printf("\n");
    
    DoAnalysis(argc, argv);
    FreeVariables();

    printf("Finished!\n");
    printf("==============================================================\n");
    printf("\n");

    return 0;
}


void Initialization(int argc, const char *argv[]) {
    analysisList = calloc(NumofElem, sizeof(int));
    files = calloc(NumofFileType, sizeof(struct FileStr*));
    
    for (int i = 1; i < argc; ++i) {
        if (strncmp(argv[i], "-", 1) == 0) {
            if (strcmp(argv[i], "-path") == 0) {
                sprintf(path, "%s", argv[i + 1]);
            } else if (strcmp(argv[i], "-rPBC") == 0) {
                analysisList[RemovePBC] = 1;
            } else if (strcmp(argv[i], "-HB") == 0) {
                analysisList[AnalyzeHB] = 1;
            } else if (strcmp(argv[i], "-En") == 0) {
                analysisList[AnalyzeEn] = 1;
            } else if (strcmp(argv[i], "-Ag") == 0) {
                analysisList[AnalyzeAgg] = 1;
            } else if (strcmp(argv[i], "-RG") == 0) {
                analysisList[CalculateRG] = 1;
            } else if (strcmp(argv[i], "-MSD") == 0) {
                analysisList[CalculateMSD] = 1;
            } else if (strcmp(argv[i], "-RMSD") == 0) {
                analysisList[CalculateRMSD] = 1;
            } else if (strcmp(argv[i], "-PotMap") == 0) {
                analysisList[CalculatePotMap] = 1;
            } else if (strcmp(argv[i], "-Ramach") == 0) {
                analysisList[AnalyzeRamach] = 1;
            } else if (strcmp(argv[i], "-ConMap") == 0) {
                analysisList[AnalyzeConMap] = 1;
            } else if (strcmp(argv[i], "-trj") == 0) {
                sprintf(names[inTrj], "%s", argv[i + 1]);
            } else if (strcmp(argv[i], "-cnt") == 0) {
                sprintf(names[inCnt], "%s", argv[i + 1]);
            } else if (strcmp(argv[i], "-sys") == 0) {
                sprintf(names[inSys], "%s", argv[i + 1]);
            } else if (strcmp(argv[i], "-log") == 0) {
                sprintf(names[inLog], "%s", argv[i + 1]);
            } else if (strcmp(argv[i], "-ref") == 0) {
                while (strncmp(argv[i + RMSDFile.count + 1], "-", 1) != 0) {
                    sprintf(RMSDFile.list[RMSDFile.count], "%s", argv[i + (++RMSDFile.count)]);
                }
            } else if (strcmp(argv[i], "-obs") == 0) {
                sprintf(obstDir, "%s", argv[i + 1]);
            } else if (strcmp(argv[i], "-nPP") == 0) {
                sprintf(targetPeptideNum, "%s", argv[i + 1]);
                nPP = atoi(argv[i + 1]);
            } else if (strcmp(argv[i], "-sum") == 0 ) {
                freshStart = 1;
            } else if (strcmp(argv[i], "-REMD") == 0) {
                RE.mark = 1;
            } else if (strcmp(argv[i], "-reNo") == 0) {
                if (i + 1 >= argc || argv[i + 1][0] == '-') {
                    printf("!!ERROR!!: Please provide a valid replica number for this REMD set!\n\n");
                    goto help;
                } else {
                    RE.numReplica = atoi(argv[i + 1]);
                }
            } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
            help:
                printf("Data Analysis Program:\n");
                printf("-path       exact path of data folder\n");
                printf("-trj        trajectory input file\n");
                printf("-cnt        connect map input file\n");
                printf("-sys        (optional) system info input file\n");
                printf("-log        (optional) log input file\n");
                printf("-ref        (optional) reference structure input file for RMSD\n");
                printf("-obs        (optional) exact path of obstruction info input file\n");
                printf("-nPP        (optional) only analyze this specific number of peptide\n");
                printf("-sum        (optional) analyze all the files in the \"path\"\n");
                printf("-rPBC       remove PBC, call at the beginning. require -trj input file\n");
                printf("-HB         analyze HB info. require -cnt input file\n");
                printf("-En         analyze energy info\n");
                printf("-Ag         analyze aggregation info\n");
                printf("-RG         calculate radius of gyration\n");
                printf("-MSD        calculate MSD\n");
                printf("-RMSD       (require -ref) calculate root-mean-square deviation\n");
                printf("-PotMap     (require -ref) calculate PE surface WRT RMSD (need 2x refs)\n");
                printf("-Ramach     (require -HB)  plot Ramachandran plots for each amino acid\n");
                printf("-ConMap     (require -HB)  plot contact maps for atoms and amino acids\n");
                printf("-REMD       analyze REMD data. follow the flags below:\n");
                printf("  -reNo     number of replicas\n");
                printf("  -bin      (optional) number of bins in histogram. default 10\n");
                printf("  -aHB      (optional) use alpha HB number as the reaction coordinate\n");
                printf("  -bHB      (optional) use beta HB number as the reaction coordinate\n");
                printf("  -tHB      (optional) use total HB number as the reaction coordinate\n");
                printf("  -rmsd     (optional) use RMSD as the reaction coordinate (require -RMSD & -ref)\n");
                printf("  -maxE     (optional) max value of potential energy (most negative). default -150\n");
                printf("  -minE     (optional) min value of potential energy (least negative). default 0\n");
                printf("  -maxRC    (optional) max value of reaction coordinate. default 10\n");
                printf("  -minRC    (optional) min value of reaction coordinate. default 0\n");
                printf("  -rate     (optional) data dumping rate, unit is time. default 10\n");
                printf("  -temp     (optional) target temperature. default 300\n");
                printf("  -st       (optional) start frame of analyzing. default 0\n");
                printf("  -et       (optional) end frame of analyzing. default -1 (the final frame)\n");
                exit(EXIT_SUCCESS);
            } else {
                if (!(strcmp(argv[i], "-bin")   == 0 ||
                      strcmp(argv[i], "-maxE")  == 0 ||
                      strcmp(argv[i], "-minE")  == 0 ||
                      strcmp(argv[i], "-rate")  == 0 ||
                      strcmp(argv[i], "-temp")  == 0 ||
                      strcmp(argv[i], "-st")    == 0 ||
                      strcmp(argv[i], "-et")    == 0 ||
                      strcmp(argv[i], "-aHB")   == 0 ||
                      strcmp(argv[i], "-bHB")   == 0 ||
                      strcmp(argv[i], "-tHB")   == 0 ||
                      strcmp(argv[i], "-rmsd")  == 0 ||
                      strcmp(argv[i], "-minRC") == 0 ||
                      strcmp(argv[i], "-maxRC") == 0)) {
                    printf ("!!ERROR!!: invalid flag: %s!\n %s:%i", argv[i], __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    
    if (RE.mark == 1 && RE.numReplica < 0) {
        printf("!!ERROR!!: Please provide a valid replica number for this REMD set!\n\n");
        goto help;
    }
    
    for (int i = 0; i < NumofFileType; i ++) {
        if (RE.mark) {
            files[i] = (struct FileStr*)calloc(RE.numReplica, sizeof(struct FileStr));
        } else {
            files[i] = (struct FileStr*)calloc(1, sizeof(struct FileStr));
        }
    }
    InitializeFiles(NumofFileType, ((RE.mark) ? RE.numReplica : 1));
    
    return;
}


void DoAnalysis(int argc, const char *argv[]) {
    for (int i = 0; i < ((RE.mark) ? RE.numReplica : 1); i ++) {
        if (RE.mark) {
            printf("Replica #%i:\n", i);
        }
        
        SystemInformationInput(i);
        AssignFileList(i);
        ReadLog();
        EstCell();
        
        for (int n = 1; n < NumofElem; n ++) {
            if (analysisList[n]) {
                switch (n) {
                    case RemovePBC:
                        AdjustPBC(i);
                        break;
                        
                    case AnalyzeHB:
                        HBInfo(i);
                        break;
                        
                    case AnalyzeEn:
                        EnergyInfo(i);
                        break;
                        
                    case AnalyzeAgg:
                        ClusterInfo(i);
                        break;
                        
                    case CalculateRG:
                        RGInfo(i);
                        break;
                        
                    case CalculateMSD:
                        break;
                        
                    case CalculateRMSD:
                        RMSDInfo(i);
                        break;
                        
                    case CalculatePotMap:
                        PESurfaceInfo(i);
                        break;
                        
                    default:
                        break;
                }
            }
        }
        printf("\n");
    }
    
    if (RE.mark) {
        WHAM(argc, argv);
        
        if (analysisList[CalculatePotMap]) {
            REMDPESurfaceCombine();
        }
    }
}

