//
//  REMD.c
//  Analysis
//
//  Created by Size Zheng on 11/18/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"


int replicaNo;
int exRate, binNum = 10;
int aHB = 0, bHB = 0, tHB = 0;
long **count; //[replica number/temperature number][bin number]
long *totCount;
long st = 0, et = -1;
double dumpRate = 10.0;
double *temperature;
double minPotential = 0, maxPotential = -200;
double potenInterval;
double targetT = 300;
double **probability; //[replica number/temperature number][bin number]
double *potenE; //[bin number]
double *DoS; //[bin number]
double *fValue;


int FindRightElem(double value_i, double *list, int listSize);
int CheckCount(long **count, long *total, int repNo, int binNum);
int CountPotenE(FILE *inputFile, FILE *inReCorFile, int rate, long startTime, long endTime, long *thisCount, long *totCount, double *pE);
int CountPotenEPosition(double curPotenE);
int CheckValidity(double);
double FindCurT(FILE *inputFile, long *step, long startTime, long endTime);
double CalBeta(double T);
double CalAvePoten(double tgtT, double *dos, double *pE, int binNum, int dimension);
double CalCv(double tgtT);
double FindOffSet(double *dos, int binNum);
void CountReplica(long startTime, long endTime);
void ReadConfigFile(void);
void PrintCountinBin(int elem);
void PrintProbability0(void);
void PrintProbability(int id);
void CalDoS(double *dos, long **n, long *totN, double *pE, int repNo, int binNo);
void CalProbability(double *dos, double *f, double *pE, int replicaNo, int binNum);
void CalFreeEnergy(double *tgtT, int repNo, int binNum, double *dos, double *pE);
void CalEntropy(double *dos, int binNum);
void SaveCvT(int bin, double minT, double maxT);
void SavePotenT(int bin, double minT, double maxT);
void SaveProbability(int reNo, int binNum);
void FreeVariable(void);
static void InitializeVariables(int argc, const char *argv[]);
FILE *OpenInputFile(char *name);


int WHAM(int argc, const char *argv[]) {
    
    printf("Executing WHAM analysis...");
    fflush(stdout);
    
    InitializeVariables(argc, argv);
    ReadConfigFile();
    
    CountReplica(st, et);
    if (CheckCount(count, totCount, replicaNo, binNum) < 0) {
        printf("!!ERROR!!: count results have something wrong! %s:%i\n", __FILE__, __LINE__);
    }
    
    CalDoS(DoS, count, totCount, potenE, replicaNo, binNum);
    CalFreeEnergy(temperature, replicaNo, binNum, DoS, potenE);
    CalEntropy(DoS, binNum);
    if (!(aHB || bHB || tHB)) {
        CalProbability(DoS, fValue, potenE, replicaNo, binNum);
        SaveCvT(binNum, temperature[0], temperature[replicaNo - 1]);
        SavePotenT(binNum, temperature[0], temperature[replicaNo - 1]);
    }

    FreeVariable();
    printf("Done!\n");
    return 0;
}

static void InitializeVariables(int argc, const char *argv[]) {
    int checkFlag = 1;
    
    replicaNo = RE.numReplica;
    for (int i = 1; i < argc; ++i) {
        if (checkFlag) {
            if (strcmp(argv[i], "-REMD") == 0) {
                checkFlag = 0;
            }
            continue;
        }
        
        if (strncmp(argv[i], "-", 1) == 0) {
            if (strcmp(argv[i], "-bin") == 0) {
                binNum = atoi(argv[i + 1]); //bin number
            } else if (strcmp(argv[i], "-aHB") == 0) {
                aHB = 1; //use alpha HB as the reaction coordinate
            } else if (strcmp(argv[i], "-bHB") == 0) {
                bHB = 1; //use beta HB as the reaction coordinate
            } else if (strcmp(argv[i], "-tHB") == 0) {
                tHB = 1;
            } else if (strcmp(argv[i], "-maxE") == 0) {
                maxPotential = atof(argv[i + 1]);
                if (maxPotential > 0) {
                    maxPotential *= -1;
                }
            } else if (strcmp(argv[i], "-minE") == 0) {
                minPotential = atof(argv[i + 1]);
                if (minPotential > 0) {
                    minPotential *= -1;
                }
            } else if (strcmp(argv[i], "-rate") == 0) {
                dumpRate = atof(argv[i + 1]);
            } else if (strcmp(argv[i], "-temp") == 0) {
                targetT = atof(argv[i + 1]);
            } else if (strcmp(argv[i], "-st") == 0) {
                st = atol(argv[i + 1]);
            } else if (strcmp(argv[i], "-et") == 0) {
                et = atol(argv[i + 1]);
            }
        }
    }
    
    if (et > 0 && et <= st) {
        printf("!!ERROR!!: end time should be larger than the start time! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    potenInterval = fabs((maxPotential - minPotential) / binNum);
    count       = (long   **)calloc(replicaNo, sizeof(long   *));
    probability = (double **)calloc(replicaNo, sizeof(double *));
    for (int i = 0; i < replicaNo; i ++) {
        count[i]       = (long   *)calloc(binNum, sizeof(long  ));
        probability[i] = (double *)calloc(binNum, sizeof(double));
    }
    
    potenE   = (double *)calloc(binNum,    sizeof(double));
    DoS      = (double *)calloc(binNum,    sizeof(double));
    totCount = (long   *)calloc(replicaNo, sizeof(long));
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
    FILE *inEnergyFile, *inTempFile, *inReCorFile = NULL;

    for (int id = 0; id < replicaNo; id ++) {
        int elem = -1;
        long step = startTime;
        double curT;
        
        inEnergyFile = OpenInputFile(files[outPE   ][id].name);
        inTempFile   = OpenInputFile(files[inREMDTemp][id].name);
        if (aHB || bHB || tHB) {
            inReCorFile = OpenInputFile(files[outHBInfo][id].name);
        }

        while (step <= endTime || endTime < 0) {
            curT = FindCurT(inTempFile, &step, startTime, endTime);
            if (curT < 0) {
                break; //reach the end of the file
            }
            
            if ((elem = FindRightElem(curT, temperature, replicaNo)) < 0) {
                printf("!!ERROR!!: cannot find the right temperature! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            if (CountPotenE(inEnergyFile, inReCorFile, exRate, step, endTime, count[elem], &totCount[elem], potenE) < 0)
                break; //reach the end of the file
            step += exRate * dumpRate;
        }
        
        fclose(inEnergyFile);
        fclose(inTempFile);
    }
    
    if (aHB || bHB || tHB) {
        for (int l = 0; l < binNum; l ++) {
            long sum = 0;
            for (int i = 0; i < replicaNo; i ++) {
                sum += count[i][l];
            }
            if (sum) {
                potenE[l] /= sum;
            }
        }
    } else {
        for (int l = 0; l < binNum; l ++) {
            potenE[l] = minPotential - (l + 0.5) * potenInterval;
        }
    }
    
    return;
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
    
    free(newFValue);
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
    
    SaveProbability(reNo, binNum);
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


int CountPotenE(FILE *inputFile, FILE *inReCorFile, int rate, long startTime, long endTime, long *thisCount, long *totCount, double *pE) {
    char buffer[256];
    int check = 0;
    double curFrame = 0;
    double curPotenE;
    
    check = fscanf(inputFile, "%lf%lf", &curFrame, &curPotenE);
    while ((long)curFrame < startTime && (check = fscanf(inputFile, "%lf%lf", &curFrame, &curPotenE)) != EOF);
    
    if (check < 0 || (endTime > 0 && (long)curFrame > endTime)) {
        return -1; // end of the file
    }
    
    if (inReCorFile) {
        check = fscanf(inReCorFile, "%lf%[^\n]\n", &curFrame, buffer);
        while ((long)curFrame < startTime && (check = fscanf(inReCorFile, "%lf%[^\n]\n", &curFrame, buffer)) != EOF);
        
        if (check < 0 || (endTime > 0 && (long)curFrame > endTime)) {
            return -1; // end of the file
        }
        
        int HBNum = 0;
        if (aHB) {
            sscanf(buffer, "%i", &HBNum);
        } else if (bHB) {
            sscanf(buffer, "%i%i%i%i", &HBNum, &HBNum, &HBNum, &HBNum);
        } else if (tHB) {
            sscanf(buffer, "%i%i%i%i%i", &HBNum, &HBNum, &HBNum, &HBNum, &HBNum);
        } else {
            printf("!!ERROR!!: invalid reaction coordinate option! %s:%i\n", __FILE__, __LINE__);
        }
        pE[HBNum] += curPotenE;
        thisCount[HBNum] ++;
        (*totCount) ++;
        
        for (int i = 1; i < rate; i ++) {
            if (fscanf(inputFile, "%s%lf", buffer, &curPotenE) == EOF)
                return -1; // the file could be end at the middle of calculations
            if (fscanf(inReCorFile, "%lf%[^\n]\n", &curFrame, buffer) == EOF)
                return -1; // the file could be end at the middle of calculations
            
            if (aHB) {
                sscanf(buffer, "%i", &HBNum);
            } else if (bHB) {
                sscanf(buffer, "%i%i%i%i", &HBNum, &HBNum, &HBNum, &HBNum);
            } else if (tHB) {
                sscanf(buffer, "%i%i%i%i%i", &HBNum, &HBNum, &HBNum, &HBNum, &HBNum);
            } else {
                printf("!!ERROR!!: invalid reaction coordinate option! %s:%i\n", __FILE__, __LINE__);
            }
            pE[HBNum] += curPotenE;
            thisCount[HBNum] ++;
            (*totCount) ++;
        }
        
    } else {
        int position = CountPotenEPosition(curPotenE);
        thisCount[position] ++;
        (*totCount) ++;
        
        for (int i = 1; i < rate; i ++) {
            if (fscanf(inputFile, "%s%lf", buffer, &curPotenE) == EOF)
                return -1; // the file could be end at the middle of calculations
            position = CountPotenEPosition(curPotenE);
            thisCount[position] ++;
            (*totCount) ++;
        }
    }

    return 0;
}


int CountPotenEPosition(double curPotenE) {
    int num = 1;
    while (curPotenE < minPotential - num * potenInterval) {
        num ++;
    }
    
    if (num > binNum) {
        printf("!!ERROR!!: potential from input file excesses the assigned max potential! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return num - 1;
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
    fscanf(inConfigFile, "%s%s%i", buffer, buffer, &exRate);
    exRate /= dumpRate;
    
    fclose(inConfigFile);
    
    if (st > 0) {
        st = (st % exRate == 0) ? st : (st / exRate + 1) * exRate;
    }
    if (et > 0) {
        et = (et % exRate == 0) ? et : (et / exRate + 1) * exRate;
    }
    return;
}

void PrintCountinBin(int elem) {
    printf("Elem #%i:\n", elem);
    for (int i = 0; i < binNum; i ++) {
        printf("bin #%i: %li\n", i, count[elem][i]);
    }
    
    return;
}

void PrintProbability0() {
    for (int i = 0; i < binNum; i ++) {
        printf("Bin #%i: %e\n", i, DoS[i]);
    }
    
    return;
}

void PrintProbability(int id) {
    for (int i = 0; i < binNum; i ++) {
        printf("Bin #%i: %lf\n", i, probability[id][i]);
    }
    
    return;
}

void SaveProbability(int reNo, int binNum) {
    char directory[1024];
    double x, y;
    FILE *outputFile;
    
    for (int i = 0; i < reNo; i ++) {
        sprintf(directory, "%sProbability%i.txt", path, i);
        outputFile = fopen(directory, "w");
        
        for (int l = 0; l < binNum; l ++) {
            x = minPotential - (l + 0.5) * potenInterval;
            y = probability[i][l];
            
            if (CheckValidity(y) && probability[i][l] > 1E-4) {
                fprintf(outputFile, "%10.4lf%20.4e\n", x, y);
            }
        }
        
        fclose(outputFile);
    }
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
    
    for (int i = 0; i < repNo; i ++) {
        sprintf(directory, "%sFreeEnergy_%.0lf.txt", path, tgtT[i]);
        outputFile = fopen(directory, "w");
        
        for (int l = 0; l < binNum; l ++) {
            if (dos[l] != 0 && pE[l] != 0) {
                x = l;
                y = pE[l] - (log(dos[l]) - FindOffSet(dos, binNum)) / CalBeta(tgtT[i]);
                fprintf(outputFile, "%10.4lf%20.4e\n", x, y);
            }
        }
        
        fclose(outputFile);
    }
    return;
}


void CalEntropy(double *dos, int binNum) {
    char directory[1024];
    double x, y;
    FILE *outputFile;
    
    sprintf(directory, "%sEntropy.txt", path);
    outputFile = fopen(directory, "w");
    
    for (int l = 0; l < binNum; l ++) {
        if (aHB || bHB || tHB) {
            x = l;
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
    
     U = CalAvePoten(tgtT, DoS, potenE, binNum, 1);
    U2 = CalAvePoten(tgtT, DoS, potenE, binNum, 2);
    
    denominator = BOLTZMANN * tgtT;
    denominator *= denominator;
    
    if (CheckValidity(U2) && CheckValidity(U) && CheckValidity(denominator)) {
        return (U2 - U * U) / denominator;
    } else {
        return 0;
    }
}

void SaveCvT(int bin, double minT, double maxT) {
    double gapT, tgtT;
    char directory[1024];
    FILE *outputFile;
    
    sprintf(directory, "%sCvT.txt", path);
    outputFile = fopen(directory, "w");
    
    gapT = (maxT - minT) / bin;
    
    for (int i = 0; i <= bin; i ++) {
        tgtT = minT + i * gapT;
        fprintf(outputFile, "%10.2lf%20.4lf\n", minT + i * gapT, CalCv(tgtT));
    }
    
    fclose(outputFile);
}


void SavePotenT(int bin, double minT, double maxT) {
    double gapT, tgtT;
    char directory[1024];
    FILE *outputFile;
    
    sprintf(directory, "%sUT.txt", path);
    outputFile = fopen(directory, "w");
    
    gapT = (maxT - minT) / bin;
    
    for (int i = 0; i <= bin; i ++) {
        tgtT = minT + i * gapT;
        fprintf(outputFile, "%10.2lf%20.4lf\n", minT + i * gapT,
                CalAvePoten(tgtT, DoS, potenE, bin, 1));
    }
    
    fclose(outputFile);
}

void FreeVariable(void) {
    for (int i = 0; i < replicaNo; i ++) {
        free(count[i]);
        free(probability[i]);
    }
    free(count);
    free(probability);
    free(potenE);
    free(DoS);
    free(totCount);
    free(fValue);
    free(temperature);
    
    return;
}


int CheckValidity(double value) {
    if (isinf(value) || value != value) {
        return 0;
    }
    
    return 1;
}
