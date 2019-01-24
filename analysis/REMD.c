//
//  REMD.c
//  Analysis
//
//  Created by Size Zheng on 11/18/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"

struct RCOptionStr { //reaction coordiate option
    int aHB, bHB, tHB;
    int rmsd;
};

struct REMDAnaStr {
    int replicaNo;
    int exRate, binNum;
    long st, et;
    double dumpRate;
    double targetT;
    double min_E, max_E;
    double min_RC, max_RC;
    double gap_E, gap_RC;
    struct RCOptionStr option;
};

static long **count; //[replica number/temperature number][bin number]
static long *totCount;
static double *temperature;
static double **probability; //[replica number/temperature number][bin number]
static double *potenE; //[bin number]
static double *DoS; //[bin number]
static double *fValue;
static struct REMDAnaStr REMD = {
    .binNum = 10,
    .exRate = 10,
    .option.aHB  = 0,
    .option.bHB  = 0,
    .option.tHB  = 0,
    .option.rmsd = 0,
    .st = 0,
    .et = -1,
    .dumpRate = 10,
    .targetT = 300,
    .min_E = 0,
    .max_E = 200,
    .min_RC = 0,
    .max_RC = 10
};


static int FindRightElem(double value_i, double *list, int listSize);
static int CheckCount(long **count, long *total, int repNo, int binNum);
static int CountBin(FILE *inEnergFile, FILE *inReCorFile, char *type, int rate, long startTime, long endTime, long *thisCount, long *totCount, double *pE);
static int CountPosition(double value, double min, double gap);
static int CheckValidity(double);
static double FindCurT(FILE *inputFile, long *step, long startTime, long endTime);
static double CalBeta(double T);
static double CalAvePoten(double tgtT, double *dos, double *pE, int binNum, int dimension);
static double CalCv(double tgtT);
static double FindOffSet(double *dos, int binNum);
static void CountReplica(long startTime, long endTime);
static void ReadConfigFile(void);
static void CalDoS(double *dos, long **n, long *totN, double *pE, int repNo, int binNo);
static void CalProbability(double *dos, double *f, double *pE, int replicaNo, int binNum);
static void CalFreeEnergy(double *tgtT, int repNo, int binNum, double *dos, double *pE);
static void CalEntropy(double *dos, int binNum);
static void CalCvPotT(int bin, double minT, double maxT);
static void SaveProbability(int repNo, int binNum, double *tgtT);
static void FreeVariable(void);
static void InitializeVariables(int argc, const char *argv[]);
static FILE *OpenInputFile(char *name);
void PrintCountinBin(int elem); //for debug
void PrintDoS(void);   //for debug
void PrintProbability(int id);  //for debug


int WHAM(int argc, const char *argv[]) {
    
    printf("Executing WHAM analysis...");
    fflush(stdout);
    
    InitializeVariables(argc, argv);
    ReadConfigFile();
    
    /*
     These two functions, CountReplica() and CalDoS(), have to be called
     to calculate counts in bins and density of state.
     Here bins can be related to any reaction coordinates.
     Right now the reaction coordinates can be hydrogen number and RMSD.
     */
    CountReplica(REMD.st, REMD.et);
    if (CheckCount(count, totCount, REMD.replicaNo, REMD.binNum) < 0) {
        printf("!!ERROR!!: count results have something wrong! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    CalDoS(DoS, count, totCount, potenE, REMD.replicaNo, REMD.binNum);
    
    CalFreeEnergy(temperature, REMD.replicaNo, REMD.binNum, DoS, potenE);
    CalEntropy(DoS, REMD.binNum);
    CalProbability(DoS, fValue, potenE, REMD.replicaNo, REMD.binNum);

    if (!(REMD.option.aHB ||
          REMD.option.bHB ||
          REMD.option.tHB ||
          REMD.option.rmsd)) {
        CalCvPotT(REMD.binNum, temperature[0], temperature[REMD.replicaNo - 1]);
    }

    FreeVariable();
    printf("Done!\n");
    return 0;
}

static void InitializeVariables(int argc, const char *argv[]) {
    int checkFlag = 1;
    
    REMD.replicaNo = RE.numReplica;
    for (int i = 1; i < argc; ++i) {
        if (checkFlag) {
            if (strcmp(argv[i], "-REMD") == 0) {
                checkFlag = 0;
            }
            continue;
        }
        
        if (strncmp(argv[i], "-", 1) == 0) {
            if (strcmp(argv[i], "-bin") == 0) {
                REMD.binNum = atoi(argv[i + 1]); //bin number
            } else if (strcmp(argv[i], "-aHB") == 0) {
                REMD.option.aHB = 1; //use alpha HB as the reaction coordinate
            } else if (strcmp(argv[i], "-bHB") == 0) {
                REMD.option.bHB = 1; //use beta HB as the reaction coordinate
            } else if (strcmp(argv[i], "-tHB") == 0) {
                REMD.option.tHB = 1; //use total HB as the reaction coordinate
            } else if (strcmp(argv[i], "-rmsd") == 0) {
                REMD.option.rmsd = 1; //use RMSD as the reaction coordinate
                if (argv[i + 1][0] != '-') {
                    REMD.option.rmsd = atoi(argv[i + 1]);
                }
            } else if (strcmp(argv[i], "-maxE") == 0) {
                REMD.max_E = atof(argv[i + 1]);
                if (REMD.max_E > 0) {
                    REMD.max_E *= -1;
                }
            } else if (strcmp(argv[i], "-minE") == 0) {
                REMD.min_E = atof(argv[i + 1]);
                if (REMD.min_E > 0) {
                    REMD.min_E *= -1;
                }
            } else if (strcmp(argv[i], "-minRC") == 0) {
                REMD.min_RC = atof(argv[i + 1]);
            } else if (strcmp(argv[i], "-maxRC") == 0) {
                REMD.max_RC = atof(argv[i + 1]);
            } else if (strcmp(argv[i], "-rate") == 0) {
                REMD.dumpRate = atof(argv[i + 1]);
            } else if (strcmp(argv[i], "-temp") == 0) {
                REMD.targetT = atof(argv[i + 1]);
            } else if (strcmp(argv[i], "-st") == 0) {
                REMD.st = atol(argv[i + 1]);
            } else if (strcmp(argv[i], "-et") == 0) {
                REMD.et = atol(argv[i + 1]);
            }
        }
    }
    
    if (REMD.et > 0 && REMD.et <= REMD.st) {
        printf("!!ERROR!!: end time should be larger than the start time! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (REMD.option.aHB || REMD.option.bHB || REMD.option.tHB) {
        REMD.binNum = REMD.max_RC - REMD.min_RC;
    }
    
    REMD.gap_E  = fabs((REMD.max_E  - REMD.min_E)  / REMD.binNum);
    REMD.gap_RC = fabs((REMD.max_RC - REMD.min_RC) / REMD.binNum);
      LONG_2CALLOC(count,       REMD.replicaNo, REMD.binNum);
    DOUBLE_2CALLOC(probability, REMD.replicaNo, REMD.binNum);
    
    potenE   = (double *)calloc(REMD.binNum,    sizeof(double));
    DoS      = (double *)calloc(REMD.binNum,    sizeof(double));
    totCount = (long   *)calloc(REMD.replicaNo, sizeof(long));
}


FILE *OpenInputFile(char *name) {
    char directory[1024], buffer[1024], letter;
    FILE *inputFile;
    
    sprintf(directory, "%s%s", path, name);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: the file %s does not exist! please check your directory!\n", directory);
        printf("           you may need to implement some other analyses before doing REMD analysis!\n");
        exit(EXIT_FAILURE);
    }
    while (fscanf(inputFile, "%c", &letter) != EOF && letter == '@') {
        fscanf(inputFile, "%[^\n]\n", buffer);
    }
    fseek(inputFile, -1, SEEK_CUR);
    
    return inputFile;
}


void CountReplica(long startTime, long endTime) {
    char type[256];
    FILE *inEnergyFile, *inTempFile, *inReCorFile = NULL;

    for (int id = 0; id < REMD.replicaNo; id ++) {
        int elem = -1;
        long step = startTime;
        double curT;
        
        inEnergyFile = OpenInputFile(files[outEne    ][id].name);
        inTempFile   = OpenInputFile(files[inREMDTemp][id].name);
        
        //set the reaction coordinate file and type
        if (REMD.option.aHB ||
            REMD.option.bHB ||
            REMD.option.tHB) {
            inReCorFile = OpenInputFile(files[outHBInfo][id].name);
            sprintf(type, "hb");
        } else if (REMD.option.rmsd) {
            inReCorFile = OpenInputFile(files[outRMSD  ][id].name);
            sprintf(type, "rmsd");
        } else {
            sprintf(type, "energy");
        }

        while (step <= endTime || endTime < 0) {
            curT = FindCurT(inTempFile, &step, startTime, endTime);
            if (curT < 0) {
                break; //reach the end of the file
            }
            
            if ((elem = FindRightElem(curT, temperature, REMD.replicaNo)) < 0) {
                printf("!!ERROR!!: cannot find the right temperature! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            if (CountBin(inEnergyFile, inReCorFile, type, REMD.exRate, step, endTime, count[elem], &totCount[elem], potenE))
                break;
            
            step += REMD.exRate * REMD.dumpRate;
        }
        
        fclose(inEnergyFile);
        fclose(inTempFile);
    }
    
    if (inReCorFile) { //the reaction coordinate is assigned
        for (int l = 0; l < REMD.binNum; l ++) {
            long sum = 0;
            for (int i = 0; i < REMD.replicaNo; i ++) {
                sum += count[i][l];
            }
            if (sum) {
                potenE[l] /= sum;
            }
        }
    } else {
        for (int l = 0; l < REMD.binNum; l ++) {
            potenE[l] = REMD.min_E - (l + 0.5) * REMD.gap_E;
        }
    }
    
    return;
}


int CountBin(FILE *inEnergFile, FILE *inReCorFile, char *type, int rate, long startTime, long endTime, long *thisCount, long *totCount, double *pE) {
    int check = 0;
    double *curFrame_E  = &energy.step;
    double *curPotenE = &energy.PE;
    
    check = ReadEnergyFile(inEnergFile, &energy);
    while ((long)(*curFrame_E) < startTime && (check = ReadEnergyFile(inEnergFile, &energy)) == 0);
    if (check < 0 || (endTime > 0 && (long)(*curFrame_E) > endTime)) {
        return -1; // end of the file
    }
    
    if (inReCorFile) {
        if (strcmp(type, "hb") == 0) {
            int HBNum = 0;
            double *curFrame_HB  = &hb.step;
            
            check = ReadHBFile(inReCorFile, &hb);
            while ((long)(*curFrame_HB) < startTime && (check = ReadHBFile(inReCorFile, &hb)) == 0);
            if (check < 0 || (endTime > 0 && (long)(*curFrame_HB) > endTime)) {
                return -1; // end of the file
            }
            
            if (REMD.option.aHB) HBNum = hb.alpha;
            else if (REMD.option.bHB) HBNum = hb.beta;
            else if (REMD.option.tHB) HBNum = hb.totl;
            else {
                printf("!!ERROR!!: invalid reaction coordinate option! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            pE[HBNum] += *curPotenE;
            thisCount[HBNum] ++;
            (*totCount) ++;
            
            for (int i = 1; i < rate; i ++) {
                if (ReadEnergyFile(inEnergFile, &energy) < 0)
                    return -1; // the file could be end at the middle of calculations
                if (ReadHBFile(inReCorFile, &hb) < 0)
                    return -1; // the file could be end at the middle of calculations
                
                if (REMD.option.aHB) HBNum = hb.alpha;
                else if (REMD.option.bHB) HBNum = hb.beta;
                else if (REMD.option.tHB) HBNum = hb.totl;
                else {
                    printf("!!ERROR!!: invalid reaction coordinate option! %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                
                pE[HBNum] += *curPotenE;
                thisCount[HBNum] ++;
                (*totCount) ++;
            }
        } else if (strcmp(type, "rmsd") == 0) {
            double *curFrame_rmsd = &rmsd.step;
            
            check = ReadRMSDFile(inReCorFile, &rmsd, REMD.option.rmsd);
            while ((long)(*curFrame_rmsd) < startTime && (check = ReadRMSDFile(inReCorFile, &rmsd, REMD.option.rmsd)) == 0);
            if (check < 0 || (endTime > 0 && (long)(*curFrame_rmsd) > endTime)) {
                return -1;
            }
            
            int position = CountPosition(rmsd.vRMSD, REMD.min_RC, REMD.gap_RC);
            pE[position] += *curPotenE;
            thisCount[position] ++;
            (*totCount) ++;
            
            for (int i = 1; i < rate ; i ++) {
                if (ReadEnergyFile(inEnergFile, &energy) < 0)
                    return -1;
                if (ReadRMSDFile(inReCorFile, &rmsd, REMD.option.rmsd) < 0)
                    return -1;
                
                position = CountPosition(rmsd.vRMSD, REMD.min_RC, REMD.gap_RC);
                pE[position] += *curPotenE;
                thisCount[position] ++;
                (*totCount) ++;
            }
        }
    } else if (!inReCorFile && strcmp(type, "energy") == 0) {
        int position = CountPosition(*curPotenE, REMD.min_E, REMD.gap_E);
        thisCount[position] ++;
        (*totCount) ++;
        
        for (int i = 1; i < rate; i ++) {
            if (ReadEnergyFile(inEnergFile, &energy) < 0) return -1; // the file could be end at the middle of calculations
            
            position = CountPosition(*curPotenE, REMD.min_E, REMD.gap_E);
            thisCount[position] ++;
            (*totCount) ++;
        }
    } else {
        printf("!!ERROR!!: bin count option is invalid! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return 0;
}


void CalDoS(double *dos, long **n, long *totN, double *pE, int repNo, int binNo) {
    double *newFValue;
    double numerator, denominator;
    double tolerance = 10E-7;
    
    fValue = (double *)calloc(repNo, sizeof(double));
    newFValue = (double *)calloc(repNo, sizeof(double));
    
repeat:
    for (int l = 0; l < binNo; l ++) {
        if (pE[l] == 0) {
            continue;
        }
        numerator = denominator = 0;
        for (int i = 0; i < repNo; i ++) {
            numerator += n[i][l];
            denominator += totN[i] / exp(fValue[i])
                         * exp(-1 * CalBeta(temperature[i]) * pE[l] / 2)
                         * exp(-1 * CalBeta(temperature[i]) * pE[l] / 2); //avoid INF or NaN
        }
        
        if (CheckValidity(denominator)) {
            dos[l] = numerator / denominator;
        }
    }
    
    for (int i = 0; i < repNo; i ++) {
        denominator = 0;
        for (int l = 0; l < binNo; l ++) {
            if (pE[l] == 0) {
                continue;
            }
            denominator += dos[l]
                         * exp(-1 * CalBeta(temperature[i]) * pE[l] / 2)
                         * exp(-1 * CalBeta(temperature[i]) * pE[l] / 2); //avoid INF or NaN
        }
        
        newFValue[i] = log(denominator);
    }
    
    double error[2] = {0};
    for (int i = 0; i < repNo; i ++) {
        error[1] = fabs((newFValue[i] - fValue[i]) / fValue[i]);
        error[0] += error[1] * error[1];
    }
    
    if (error[0] > tolerance) {
        for (int n = 0; n < repNo; n ++) {
            fValue[n] = newFValue[n];
        }
        goto repeat;
    }
    
    free(newFValue); newFValue = NULL;
    return;
}


void CalProbability(double *dos, double *f, double *pE, int reNo, int binNum) {
    for (int i = 0; i < reNo; i ++) {
        for (int l = 0; l < binNum; l ++) {
            probability[i][l] = dos[l] / exp(f[i])
                              * exp(-1 * CalBeta(temperature[i]) * pE[l] / 2)
                              * exp(-1 * CalBeta(temperature[i]) * pE[l] / 2);
        }
    }
    
    SaveProbability(reNo, binNum, temperature);
    return;
}


double CalBeta(double T) {
    return 1 / (T * BOLTZMANN);
}


int FindRightElem(double value_i, double *list, int listSize) {
    double deviation = 5;
    
    for (int i = 0; i < listSize; i ++) {
        if (fabs(value_i - list[i]) < deviation) {
            return i;
        }
    }
    
    return -1;
}


double FindCurT(FILE *inputFile, long *step, long startTime, long endTime) {
    int check = 0;
    double curFrame = 0;
    double thisT;
    
    check = fscanf(inputFile, "%lf%lf", &curFrame, &thisT);
    while ((long)curFrame < startTime && (check = fscanf(inputFile, "%lf%lf", &curFrame, &thisT)) != EOF);
    
    if (check < 0 || (endTime > 0 && (long)curFrame > endTime) || thisT == 0) {
        return -1; // end of the file
    }
    
    if (*step != curFrame) {
        *step = curFrame;
    }
    
    return thisT / BOLTZMANN;
}

int CountPosition(double value, double min, double gap) {
    int num = 0;
    while (fabs(value) < min + num * gap ||
           fabs(value) > min + (num + 1) * gap) {
        num ++;
    }
    
    if (num > REMD.binNum) {
        printf("!!ERROR!!: potential from input file excesses the assigned max potential! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return num;
}


void ReadConfigFile() {
    int num;
    char directory[1024], buffer[1024];
    FILE *inConfigFile;
    
    sprintf(directory, "%sinput/REMDConfig.txt", path);
    inConfigFile = fopen(directory, "r");
    if (inConfigFile == NULL) {
        printf("!!ERROR!!: the configuration file %s does not exist! please check your directory!\n", directory);
        exit(EXIT_FAILURE);
    }
    
    fscanf(inConfigFile, "%[^\n]\n", buffer);
    fscanf(inConfigFile, "%s%s%i", buffer, buffer, &num);
    
    fscanf(inConfigFile, "%s%s", buffer, buffer);
    temperature = (double *)calloc(num, sizeof(double));
    for (int i = 0; i < num; i ++) {
        fscanf(inConfigFile, "%lf", &temperature[i]);
    }
    
    fscanf(inConfigFile, "\n");
    fscanf(inConfigFile, "%[^\n]\n", buffer);
    fscanf(inConfigFile, "%s%s%i", buffer, buffer, &REMD.exRate);
    REMD.exRate /= REMD.dumpRate;
    
    fclose(inConfigFile);
    
    if (REMD.st > 0) {
        REMD.st = (REMD.st % REMD.exRate == 0) ? REMD.st : (REMD.st / REMD.exRate + 1) * REMD.exRate;
    }
    if (REMD.et > 0) {
        REMD.et = (REMD.et % REMD.exRate == 0) ? REMD.et : (REMD.et / REMD.exRate + 1) * REMD.exRate;
    }
    return;
}

void PrintCountinBin(int elem) {
    printf("Elem #%i:\n", elem);
    for (int i = 0; i < REMD.binNum; i ++) {
        printf("bin #%i: %li\n", i, count[elem][i]);
    }
    
    return;
}

void PrintDoS(void) {
    for (int i = 0; i < REMD.binNum; i ++) {
        printf("Bin #%i: %e\n", i, DoS[i]);
    }
    
    return;
}

void PrintProbability(int id) {
    for (int i = 0; i < REMD.binNum; i ++) {
        printf("Bin #%i: %lf\n", i, probability[id][i]);
    }
    
    return;
}

void SaveProbability(int repNo, int binNum, double *tgtT) {
    char directory[1024];
    double x, y;
    FILE *outputFile;
    
    sprintf(directory, "%sProbability.txt", path);
    outputFile = fopen(directory, "w");
    
    fprintf(outputFile, "# Probability\n");
    fprintf(outputFile, "# %8s ", "binNum");
    for (int l = 0; l < repNo; l ++) {
        fprintf(outputFile, "%20.2lf (K) ", tgtT[l]);
    }
    fprintf(outputFile, "\n");
    
    for (int l = 0; l < binNum; l ++) {
        if (REMD.option.aHB ||
            REMD.option.bHB ||
            REMD.option.tHB ||
            REMD.option.rmsd) {
            x = l * REMD.gap_RC;
        } else {
            x = potenE[l];
        }
        fprintf(outputFile, "%10.4lf ", x);
        
        for (int i = 0; i < repNo; i ++) {
            y = probability[i][l];
            if (CheckValidity(y) && probability[i][l] > 1E-4) {
                fprintf(outputFile, "%20.4e ", y);
            } else {
                fprintf(outputFile, "%20.4e ", 0.0);
            }
        }
        fprintf(outputFile, "\n");
    }
    
    fclose(outputFile);
    return;
}

int CheckCount(long **count, long *total, int repNo, int binNum) {
    long sum;
    for (int i = 0; i < repNo; i ++) {
        sum = 0;
        for (int l = 0; l < binNum; l ++) {
            sum += count[i][l];
        }
        
        if (sum != total[i]) {
            return -1;
        }
    }
    
    return 0;
}

void CalFreeEnergy(double *tgtT, int repNo, int binNum, double *dos, double *pE) {
    double x, y;
    char directory[1024];
    FILE *outputFile;
    
    sprintf(directory, "%sFreeEnergy.txt", path);
    outputFile = fopen(directory, "w");

    fprintf(outputFile, "# Free Energy\n");
    fprintf(outputFile, "# %8s ", "binNum");
    for (int l = 0; l < repNo; l ++) {
        fprintf(outputFile, "%20.2lf (K) ", tgtT[l]);
    }
    fprintf(outputFile, "\n");
    
    for (int l = 0; l < binNum; l ++) {
        if (REMD.option.aHB ||
            REMD.option.bHB ||
            REMD.option.tHB ||
            REMD.option.rmsd) {
            x = l * REMD.gap_RC;
        } else {
            x = potenE[l];
        }
        fprintf(outputFile, "%10.2lf ", x);
        
        for (int i = 0; i < repNo; i ++) {
            if (dos[l] != 0 && pE[l] != 0) {
                y = pE[l] - (log(dos[l]) - FindOffSet(dos, binNum)) / CalBeta(tgtT[i]);
                fprintf(outputFile, "%20.4e ", y);
            } else {
                fprintf(outputFile, "%20.4e ", 0.0);
            }
        }
        fprintf(outputFile, "\n");
    }
    
    fclose(outputFile);
    return;
}


void CalEntropy(double *dos, int binNum) {
    char directory[1024];
    double x, y;
    FILE *outputFile;
    
    sprintf(directory, "%sEntropy.txt", path);
    outputFile = fopen(directory, "w");
    
    for (int l = 0; l < binNum; l ++) {
        if (REMD.option.aHB ||
            REMD.option.bHB ||
            REMD.option.tHB ||
            REMD.option.rmsd) {
            x = l * REMD.gap_RC;
        } else {
            x = potenE[l];
        }
        
        y = log(dos[l]) - FindOffSet(dos, binNum); //unit of kB
        
        if (CheckValidity(y)) {
            fprintf(outputFile, "%10.4lf%20.4e\n", x, y);
        }
    }
    
    fclose(outputFile);
}


double FindOffSet(double *dos, int binNum) {
    int l = -1;
    double offSet = 0;
    
    while(dos[++l] == 0);
    offSet = dos[l];
    for (; l < binNum - 1; l ++) {
        if (dos[l + 1] > 0 && dos[l] > dos[l + 1]) {
            offSet = dos[l + 1];
        }
    }
    
    return log(offSet);
}


double CalAvePoten(double tgtT, double *dos, double *pE, int binNum, int dimension) {
    double numerator, denominator;
    double PEn;
    
    numerator = denominator = 0;
    for (int l = 0; l < binNum; l ++) {
        PEn = pE[l];
        for (int i = 1; i < dimension; i ++) {
            PEn *= pE[l];
        }
        numerator   += PEn * dos[l] * exp(-1 * CalBeta(tgtT) * pE[l] / 2) * exp(-1 * CalBeta(tgtT) * pE[l] / 2); //avoid INF or NaN
        denominator +=       dos[l] * exp(-1 * CalBeta(tgtT) * pE[l] / 2) * exp(-1 * CalBeta(tgtT) * pE[l] / 2);
    }
    return numerator / denominator;
}


double CalCv(double tgtT) {
    double U, U2;
    double denominator;
    
     U = CalAvePoten(tgtT, DoS, potenE, REMD.binNum, 1);
    U2 = CalAvePoten(tgtT, DoS, potenE, REMD.binNum, 2);
    
    denominator = BOLTZMANN * tgtT;
    denominator *= denominator;
    
    if (CheckValidity(U2) && CheckValidity(U) && CheckValidity(denominator)) {
        return (U2 - U * U) / denominator;
    } else {
        return 0;
    }
}


//To calculate Cv & U vs T, it has to count bins of potential energy
void CalCvPotT(int bin, double minT, double maxT) {
    double gapT, tgtT;
    char directory[1024];
    FILE *outputFile;
    
    sprintf(directory, "%sCvPotT.txt", path);
    outputFile = fopen(directory, "w");
    fprintf(outputFile, "# Cv and PE at different T\n");
    fprintf(outputFile, "# %8s %20s %20s\n", "binNum", "U", "Cv");
    
    gapT = (maxT - minT) / bin;
    
    for (int i = 0; i <= bin; i ++) {
        tgtT = minT + i * gapT;
        fprintf(outputFile, "%10.2lf %20.4lf %20.4lf\n", minT + i * gapT, CalAvePoten(tgtT, DoS, potenE, bin, 1), CalCv(tgtT));
    }
    
    fclose(outputFile);
}


void FreeVariable(void) {
    _2FREE(count);
    _2FREE(probability);
    
    free(potenE);       potenE      = NULL;
    free(DoS);          DoS         = NULL;
    free(totCount);     totCount    = NULL;
    free(fValue);       fValue      = NULL;
    free(temperature);  temperature = NULL;
    
    return;
}


int CheckValidity(double value) {
    if (isinf(value) || value != value) {
        return 0;
    }
    
    return 1;
}


