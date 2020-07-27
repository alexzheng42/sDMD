//
//  RMSD.c
//  Analysis
//
//  Created by Size on 2018/12/20.
//  Copyright © 2018 Size Zheng. All rights reserved.
//

#include "Analysis.h"

#define DIM 3
#define RMSDRG 1.0
#define RMSDGRID 20
#define CPROD(a, b, c) {    \
    c[0] = a[1] * b[2] - a[2] * b[1]; \
    c[1] = a[2] * b[0] - a[0] * b[2]; \
    c[2] = a[0] * b[1] - a[1] * b[0]; \
}

#define MATRIXZERO(matrix, row, col)    \
memset(matrix[0], 0, sizeof(double) * row * col);

#ifndef SQRT2
#define SQRT2 sqrt(2.0)
#endif

struct PESurStr {
    int count;
    double PE;
};

struct RMSDReadStr rmsd;

static int nRow, nCol;
static int beginP, endP, offSet, thisAtomNum;
static double **curCoor, ***refCoor;
static double **omega, **om, **vh, **vk, **u;
static double *backboneMass;
struct PESurStr **EnSurGrid;

static int CheckBackbone(int num);
static void InitializeVariables(void);
static void FreeVar(void);
static void ReadRef(char *fileName, int ndim, int natoms, double **ref);
static void DoTranslation(int ndim, int natoms, double **coor);
static void DoRotation(int ndim, int natoms, double **curCoordinate, double **refCoordinate);
static void DoTranspose(double **matrix, int row, int col);
static void DoMultiply(double **LMatrix, double **RMatrix, double **result, int LRow, int LCol, int RRow, int RCol);
static void Coor2Matrix(double **matrix);
static void CalFitRotMatrix(int ndim, int natoms, double **refCoor, double **curCoor, double **R);
static void Jacobi(double **a, int n, double d[], double **v);
static void Rotate(double **a, int i, int j, int k, int l, double tau, double s);
static double CalRMSD(int ndim, int natoms, double **ref, double **cur);

void SaveMatrix(double **matrix, char *fileName);
void PrintMatrix(double **matrix, int row, int col);


void RMSDInfo(int id) {
    long step = 0, totalFrame;
    char directory[1024];
    struct SectionStr *sect;
    FILE *inputTrjFile;
    FILE *outputRMSDFile;

    /*
    if (RE.mark) {
        printf("\n");
        printf("!!ERROR!!: sAnalysis could not perform RMSD calculation during REMD analysis. Please do RMSD independently.\n");
        exit(EXIT_FAILURE);
    }
     */
    
    printf("Analyzing RMSD... ");
    fflush(stdout);
    
    InitializeVariables();
    
    if (RMSDFile.count == 0) {
        printf("!!ERROR!! no reference structure, which is required by RMSD calculation. please provide a valid reference structure by -ref! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    //there would be multiple reference structures
    for (int i = 0; i < RMSDFile.count; i ++) {
        ReadRef(RMSDFile.list[i], DIM, nRow, refCoor[i]);
    }
    
    //start to read the coordinates from each frame and calculate RMSD
    sprintf(directory, "%s%s", path, files[outRMSD][id].name);
    outputRMSDFile = fopen(directory, "w");
    fprintf(outputRMSDFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(outputRMSDFile, "%s\n", "@    yaxis  label \"RMSD (nm)\"");
    
    sprintf(directory, "%s%s", path, files[outTrj][id].name);
    inputTrjFile = fopen(directory, "r");
    if (inputTrjFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outTrj][id].name, path, __FILE__, __LINE__);
        printf("           You may need to execute -rPBC first to generate the required file!\n");
        exit(EXIT_FAILURE);
    }
    
    for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
        totalFrame = (sectInfo[sectNum].frameCount + sectInfo[sectNum].oldTime) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
            if (ReadGro(inputTrjFile)) break;
            
            PrintProcess(step);
            sect = &sectInfo[sectNum];
            
            fprintf(outputRMSDFile,   "%15.2lf ", step * sect->outputRate);
            for (int i = 0; i < RMSDFile.count; i ++) {
                fprintf(outputRMSDFile, "%15.4lf ", CalRMSD(DIM, nRow, refCoor[i], curCoor) / 10);
            }
            fprintf(outputRMSDFile, "\n");
        }
    }
    
    fclose(inputTrjFile);
    fclose(outputRMSDFile);    
    FreeVar();
    
    printf("\bDone\n");
    fflush(stdout);
    
    return;
}

void PESurfaceInfo(int id) {
    int place[2];
    long step = 0, totalFrame;
    double thisRMSD[2], thisPE;
    double gap = RMSDRG / RMSDGRID;
    char directory[512], buffer[512];
    struct SectionStr *sect;
    FILE *inputTrjFile, *inputCntFile, *inputPotFile;
    FILE *outputEnSurFile;
    
    printf("Analyzing potential energy surface... ");
    fflush(stdout);
    
    if (RMSDFile.count == 1) {
        printf("!!ERROR!!: calculate energy surface would need TWO reference structures. please provide the second one! ");
        printf("%s:%i \n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    InitializeVariables();
    
    //there would be multiple reference structures
    for (int i = 0; i < RMSDFile.count; i ++) {
        ReadRef(RMSDFile.list[i], DIM, nRow, refCoor[i]);
    }
    
    //start to read the coordinates from each frame and calculate RMSD
    sprintf(directory, "%s%s", path, files[outPotMap][id].name);
    outputEnSurFile = fopen(directory, "w");
    
    sprintf(directory, "%s%s", path, files[outTrj][id].name);
    inputTrjFile = fopen(directory, "r");
    if (inputTrjFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outTrj][id].name, path, __FILE__, __LINE__);
        printf("           You may need to execute -rPBC first to generate the required file!\n");
        exit(EXIT_FAILURE);
    }
    
    sprintf(directory, "%s%s", path, files[outEne][id].name);
    inputPotFile = fopen(directory, "r");
    if (inputPotFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outEne][id].name, path, __FILE__, __LINE__);
        printf("           You may need to execute -En first to generate the required file!\n");
        exit(EXIT_FAILURE);
    }
    
    for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
        
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inCnt][id].name, fileList[id].list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        inputCntFile = fopen(directory, "r");
        if (inputCntFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", buffer, path, __FILE__, __LINE__);
            printf("           You may need to specify the connection information file by using -cnt flag!\n");
            exit(EXIT_FAILURE);
        }
        
        totalFrame = (sectInfo[sectNum].frameCount + sectInfo[sectNum].oldTime) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
            if (ReadGro(inputTrjFile) || ReadConnectionMap(inputCntFile) || ReadEnergyFile(inputPotFile, &energy)) break;
            
            thisPE = energy.PE;
            PrintProcess(step);
            sect = &sectInfo[sectNum];
            
            thisRMSD[0] = CalRMSD(DIM, nRow, refCoor[0], curCoor) / 10;
            thisRMSD[1] = CalRMSD(DIM, nRow, refCoor[1], curCoor) / 10;
            
            for (int n = 0; n < RMSDGRID; n ++) {
                if ((thisRMSD[0] > n * gap) &&
                    (thisRMSD[0] <= (n + 1) * gap)) {
                    place[0] = n;
                    break;
                }
            }
            
            for (int n = 0; n < RMSDGRID; n ++) {
                if ((thisRMSD[1] > n * gap) &&
                    (thisRMSD[1] <= (n + 1) * gap)) {
                    place[1] = n;
                    break;
                }
            }
            
            EnSurGrid[place[0]][place[1]].count ++;
            EnSurGrid[place[0]][place[1]].PE += thisPE;
        }
        
        fclose(inputCntFile);
    }
    
    for (int i = 0; i < RMSDGRID; i ++) {
        for (int n = 0; n < RMSDGRID; n ++) {
            if (EnSurGrid[i][n].count)
                EnSurGrid[i][n].PE /= EnSurGrid[i][n].count;
            fprintf(outputEnSurFile, "%15.4lf %15.4lf %15.4lf %15i\n",
                    (i + 0.5) * gap,
                    (n + 0.5) * gap,
                    EnSurGrid[i][n].PE,
                    EnSurGrid[i][n].count);
        }
    }
    
    fclose(inputTrjFile);
    fclose(inputPotFile);
    fclose(outputEnSurFile);
    FreeVar();
    
    printf("\bDone\n");
    fflush(stdout);
    return;
}

void PESurfaceInfoPerT(void) {
    int num;
    int place[2];
    long step = 0, totalFrame;
    double thisRMSD[2], thisPE;
    double gap = RMSDRG / RMSDGRID;
    char directory[512], buffer[512];
    struct SectionStr *sect;
    struct PESurStr ***gridPerT;
    FILE *inputTrjFile, *inputCntFile, *inputPotFile, *inputREMDTFile;
    FILE *outputEnSurFile;
    
    printf("Analyzing potential energy surface per temperature... ");
    fflush(stdout);
    
    if (RMSDFile.count == 1) {
        printf("!!ERROR!!: calculate energy surface would need TWO reference structures. please provide the second one! ");
        printf("%s:%i \n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    gridPerT = (struct PESurStr ***)calloc(RE.numReplica, sizeof(struct PESurStr **));
    for (int m = 0; m < RE.numReplica; m ++) {
        gridPerT[m] = (struct PESurStr **)calloc(RMSDGRID, sizeof(struct PESurStr *));
        gridPerT[m][0] = (struct PESurStr *)calloc(RMSDGRID * RMSDGRID, sizeof(struct PESurStr));
        for (int n = 1; n < RMSDGRID; n ++) gridPerT[m][n] = gridPerT[m][0] + n * RMSDGRID;
        for (int i = 0; i < RMSDGRID; i ++) {
            for (int j = 0; j < RMSDGRID; j ++) {
                gridPerT[m][i][j].count = 0;
                gridPerT[m][i][j].PE = 0;
            }
        }
    }
    
    for (int id = 0; id < RE.numReplica; id ++) {
        
        InitializeVariables();
        
        //there would be multiple reference structures
        for (int i = 0; i < RMSDFile.count; i ++) {
            ReadRef(RMSDFile.list[i], DIM, nRow, refCoor[i]);
        }
        
        //start to read the coordinates from each frame and calculate RMSD
        sprintf(directory, "%s%s", path, files[outTrj][id].name);
        inputTrjFile = fopen(directory, "r");
        if (inputTrjFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outTrj][id].name, path, __FILE__, __LINE__);
            printf("           You may need to execute -rPBC first to generate the required file!\n");
            exit(EXIT_FAILURE);
        }
        
        sprintf(directory, "%s%s", path, files[outEne][id].name);
        inputPotFile = fopen(directory, "r");
        if (inputPotFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outEne][id].name, path, __FILE__, __LINE__);
            printf("           You may need to execute -En first to generate the required file!\n");
            exit(EXIT_FAILURE);
        }
        
        sprintf(directory, "%s%s", path, files[inREMDTemp][id].name);
        inputREMDTFile = fopen(directory, "r");
        if (inputREMDTFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inREMDTemp][id].name, path, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
            
            memset(buffer, '\0', sizeof(buffer));
            FindTargetFile(files[inCnt][id].name, fileList[id].list[sectNum + 1], buffer);
            
            sprintf(directory, "%s%s", path, buffer);
            inputCntFile = fopen(directory, "r");
            if (inputCntFile == NULL) {
                printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", buffer, path, __FILE__, __LINE__);
                printf("           You may need to specify the connection information file by using -cnt flag!\n");
                exit(EXIT_FAILURE);
            }
            
            totalFrame = (sectInfo[sectNum].frameCount + sectInfo[sectNum].oldTime) / sectInfo[sectNum].outputRate + 1;
            for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
                if (ReadGro(inputTrjFile) || ReadConnectionMap(inputCntFile) || ReadEnergyFile(inputPotFile, &energy) ||
                    (num = ReadREMDTempFile(inputREMDTFile) - 1) < 0)
                    break;
                
                thisPE = energy.PE;
                PrintProcess(step);
                sect = &sectInfo[sectNum];
                
                thisRMSD[0] = CalRMSD(DIM, nRow, refCoor[0], curCoor) / 10;
                thisRMSD[1] = CalRMSD(DIM, nRow, refCoor[1], curCoor) / 10;
                
                for (int n = 0; n < RMSDGRID; n ++) {
                    if ((thisRMSD[0] > n * gap) &&
                        (thisRMSD[0] <= (n + 1) * gap)) {
                        place[0] = n;
                        break;
                    }
                }
                
                for (int n = 0; n < RMSDGRID; n ++) {
                    if ((thisRMSD[1] > n * gap) &&
                        (thisRMSD[1] <= (n + 1) * gap)) {
                        place[1] = n;
                        break;
                    }
                }
                
                gridPerT[num][place[0]][place[1]].count ++;
                gridPerT[num][place[0]][place[1]].PE += thisPE;
            }
            
            fclose(inputCntFile);
        }
        
        fclose(inputTrjFile);
        fclose(inputPotFile);
        FreeVar();
    }
    
    for (int m = 0; m < RE.numReplica; m ++) {
        sprintf(directory, "%s%s", path, files[outPotMap_T][m].name);
        outputEnSurFile = fopen(directory, "w");
        
        for (int i = 0; i < RMSDGRID; i ++) {
            for (int n = 0; n < RMSDGRID; n ++) {
                if (gridPerT[m][i][n].count)
                    gridPerT[m][i][n].PE /= gridPerT[m][i][n].count;
                fprintf(outputEnSurFile, "%15.4lf %15.4lf %15.4lf %15i\n",
                        (i + 0.5) * gap,
                        (n + 0.5) * gap,
                        gridPerT[m][i][n].PE,
                        gridPerT[m][i][n].count);
            }
        }
        
        fclose(outputEnSurFile);
    }
    
    for (int i = 0; i < RE.numReplica; i ++) {
        free(gridPerT[i][0]);
        free(gridPerT[i]);
    }
    free(gridPerT);
    
    printf("\bDone\n");
    fflush(stdout);
    return;
}

void REMDPESurfaceCombine(void) {
    int thisCount;
    double thisPE;
    char directory[1024];
    FILE *inputFile, *outputFile;
    
    printf("Combine potential energy surface data... ");
    fflush(stdout);
    
    EnSurGrid = (struct PESurStr **)calloc(RMSDGRID, sizeof(struct PESurStr *));
    EnSurGrid[0] = (struct PESurStr *)calloc(RMSDGRID * RMSDGRID, sizeof(struct PESurStr));
    for (int n = 1; n < RMSDGRID; n ++) EnSurGrid[n] = EnSurGrid[0] + n * RMSDGRID;
    for (int i = 0; i < RMSDGRID; i ++) {
        for (int j = 0; j < RMSDGRID; j ++) {
            EnSurGrid[i][j].count = 0;
            EnSurGrid[i][j].PE = 0;
        }
    }
    
    sprintf(directory, "%s%s", path, "CombinedPotMap.txt");
    outputFile = fopen(directory, "w");
    
    for (int i = 0; i < RE.numReplica; i ++) {
        sprintf(directory, "%s%s", path, files[outPotMap][i].name);
        inputFile = fopen(directory, "r");
        if (inputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outPotMap][i].name, path, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        for (int m = 0; m < RMSDGRID; m ++) {
            for (int n = 0; n < RMSDGRID; n ++) {
                fscanf(inputFile, "%lf%lf%lf%i", &thisPE, &thisPE, &thisPE, &thisCount);
                if (thisCount) {
                    EnSurGrid[m][n].PE += thisPE * thisCount;
                    EnSurGrid[m][n].count += thisCount;
                }
            }
        }
        
        fclose(inputFile);
    }
    
    for (int i = 0; i < RMSDGRID; i ++) {
        for (int j = 0; j < RMSDGRID; j ++) {
            if (EnSurGrid[i][j].count) {
                EnSurGrid[i][j].PE /= EnSurGrid[i][j].count;
            }
            fprintf(outputFile, "%15.4lf %15.4lf %15.4lf %15i\n", /*RMSDRG - */(i + 0.5) / RMSDGRID, /*RMSDRG - */(j + 0.5) / RMSDGRID, EnSurGrid[i][j].PE, EnSurGrid[i][j].count);
            //fprintf(outputFile, "%15.4lf ", EnSurGrid[i][j].PE); //matrix form
        }
        //fprintf(outputFile, "\n");
    }

    fclose(outputFile);
    free(EnSurGrid[0]);
    free(EnSurGrid);    EnSurGrid = NULL;
    
    printf("\bDone!\n");
    fflush(stdout);
    
    return;
}

void InitializeVariables(void) {
    beginP = 1;
    endP = numofprotein;
    thisAtomNum = atomnum;
    offSet = 0;
    
    if (nPP) {
        beginP = endP = nPP;
        thisAtomNum = protein[nPP].endAtomNum - protein[nPP].startAtomNum + 1;
        offSet = protein[nPP].startAtomNum - 1;
    }
    
    nRow = 0;
    nCol = DIM;
    for (int i = 1 + offSet; i <= offSet + thisAtomNum; i ++) {
        if (CheckBackbone(i)) {
            nRow ++;
        }
    }
    
    backboneMass = (double *)calloc(nRow, sizeof(double));
    DOUBLE_2CALLOC(curCoor, nRow, nCol);
    DOUBLE_2CALLOC(omega, 2 * DIM, 2 * DIM);
    DOUBLE_2CALLOC(om,    2 * DIM, 2 * DIM);
    DOUBLE_2CALLOC(vh,    2 * DIM, 2 * DIM);
    DOUBLE_2CALLOC(vk,    2 * DIM, 2 * DIM);
    DOUBLE_2CALLOC(u,     2 * DIM, 2 * DIM);

    refCoor = (double ***)calloc(RMSDFile.count, sizeof(double **));
    for (int i = 0; i < RMSDFile.count; i ++) {
        DOUBLE_2CALLOC(refCoor[i], nRow, nCol);
    }
    
    if (analysisList[CalculatePotMap]) {
        EnSurGrid = (struct PESurStr **)calloc(RMSDGRID, sizeof(struct PESurStr *));
        EnSurGrid[0] = (struct PESurStr *)calloc(RMSDGRID * RMSDGRID, sizeof(struct PESurStr));
        for (int n = 1; n < RMSDGRID; n ++) EnSurGrid[n] = EnSurGrid[0] + n * RMSDGRID;
        for (int i = 0; i < RMSDGRID; i ++) {
            for (int j = 0; j < RMSDGRID; j ++) {
                EnSurGrid[i][j].count = 0;
                EnSurGrid[i][j].PE = 0;
            }
        }
    }
    
    return;
}

static void FreeVar(void) {
    free(backboneMass); backboneMass = NULL;
    _2FREE(curCoor);
    _2FREE(omega);
    _2FREE(om);
    _2FREE(vh);
    _2FREE(vk);
    _2FREE(u);
    
    for (int i = 0; i < RMSDFile.count; i ++) {
        _2FREE(refCoor[i]);
    }
    free(refCoor); refCoor = NULL;
    
    if (analysisList[CalculatePotMap]) {
        free(EnSurGrid[0]);
        free(EnSurGrid); EnSurGrid = NULL;
    }
}

void ReadRef(char *fileName, int ndim, int natoms, double **ref) {
    char directory[1024];
    FILE *inputRefFile;
    
    sprintf(directory, "%s%s", path, fileName);
    inputRefFile = fopen(directory, "r");
    if (inputRefFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", fileName, path, __FILE__, __LINE__);
        printf("           You may need to use flag -ref to specify the file name of the reference structure.\n");
        exit(EXIT_FAILURE);
    }
    
    ReadGro(inputRefFile);
    Coor2Matrix(ref);
    DoTranslation(ndim, natoms, ref);
    fclose(inputRefFile);
    
    return;
}

double CalRMSD(int ndim, int natoms, double **ref, double **cur) {
    double RMSD = 0;
    
    Coor2Matrix(cur);
    DoTranslation(ndim, natoms, cur);
    DoRotation(3, natoms, cur, ref);
    
    for (int i = 0; i < natoms; i ++) {
        for (int n = 0; n < ndim; n ++) {
            RMSD += (ref[i][n] - cur[i][n]) * (ref[i][n] - cur[i][n]);
        }
    }
    RMSD /= natoms;
    RMSD = sqrt(RMSD);
    
    return RMSD;
}

void DoTranslation(int ndim, int natoms, double **coor) {
    double *accuCoor = calloc(ndim, sizeof(double));
    
    for (int i = 0; i < natoms; i ++) {
        for (int n = 0; n < ndim; n ++) {
            accuCoor[n] += coor[i][n];
        }
    }
    
    for (int n = 0; n < ndim; n ++) {
        accuCoor[n] /= natoms;
        for (int i = 0; i < natoms; i ++) {
            coor[i][n] -= accuCoor[n];
        }
    }
    
    free(accuCoor); accuCoor = NULL;
    return;
}

void DoRotation(int ndim, int natoms, double **curCoordinate, double **refCoordinate) {
    double **clnCoor, **R;
    DOUBLE_2CALLOC(clnCoor, natoms, ndim);
    DOUBLE_2CALLOC(R, 3, 3);
    
    CalFitRotMatrix(3, natoms, refCoordinate, curCoordinate, R);
    DoTranspose(R, 3, 3); //the original coordinate matrix has been transposed:
                          //it is natoms x 3 instead of 3 x natoms
                          //since AB = (B^T)(A^T), so here R is transposed
    memcpy(clnCoor[0], curCoordinate[0], sizeof(double) * natoms * ndim);
    DoMultiply(clnCoor, R, curCoordinate, natoms, ndim, 3, 3);
    
    _2FREE(clnCoor);
    _2FREE(R);
    return;
}

void Coor2Matrix(double **matrix) {
    int m = 0;
    
    for (int i = 1 + offSet; i <= thisAtomNum + offSet; i ++) {
        if (CheckBackbone(i)) {
            for (int n = 1; n <= 3 ; n ++) {
                matrix[m][n - 1] = atom[i].dynamic->coordinate[n];
                backboneMass[m] = atom[i].property->mass;
            }
            m ++;
        } //only consider the backbone atoms
    }
    
    return;
}

void SaveMatrix(double **matrix, char *fileName) {
    int m = 0;
    char directory[1024];
    FILE *saveFile;
    
    sprintf(directory, "%s%s", path, fileName);
    saveFile = fopen(directory, "w");
    
    fprintf(saveFile, "Protein %i\n", nPP);
    fprintf(saveFile, "%5i\n", nRow);
    
    for (int i = 1 + offSet; i <= thisAtomNum + offSet; i ++) {
        if (!CheckBackbone(i)) continue; //only consider the backbone atoms
        fprintf(saveFile, "%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n",
                atom[i].property->sequence.aminoacidNum, atom[i].property->nameofAA, atom[i].property->name, m + 1,
                matrix[m][0] / 10, matrix[m][1] / 10, matrix[m][2] / 10);
        m ++;
    }
    
    fprintf(saveFile, "%10.5f%10.5f%10.5f\n",
            boxCurtDim[1] / 10,
            boxCurtDim[2] / 10,
            boxCurtDim[3] / 10);

    fclose(saveFile);
    return;
}

static void DoTranspose(double **matrix, int row, int col) {
    double swap;
    
    for (int i = 0; i < row; i ++) {
        for (int n = i + 1; n < col; n ++) {
            swap = matrix[i][n];
            matrix[i][n] = matrix[n][i];
            matrix[n][i] = swap;
        }
    }
    
    return;
}

static void DoMultiply(double **LMatrix, double **RMatrix, double **result, int LRow, int LCol, int RRow, int RCol) {
    if (LCol != RRow) {
        printf("!!ERROR!!: the column number of the 1st matrix is not equal to the row number of the 2nd matrix: they cannot be multiplied! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < LRow; i ++) {
        for (int j = 0; j < RCol; j ++) {
            result[i][j] = 0;
            for (int k = 0; k < LCol; k ++) {
                result[i][j] += LMatrix[i][k] * RMatrix[k][j];
            }
        }
    }
    
    return;
}

int CheckBackbone(int num) {
    if (strcmp(atom[num].property->name, "N") == 0 ||
        strcmp(atom[num].property->name, "CA") == 0 ||
        strcmp(atom[num].property->name, "C") == 0) {
        return 1;
    }
    
    return 0;
}

void PrintMatrix(double **matrix, int row, int col) {
    for (int i = 0; i < row; i ++) {
        for (int n = 0; n < col; n ++) {
            printf("%10.4lf ", matrix[i][n]);
        }
        printf("\n");
    }
}

void CalFitRotMatrix(int ndim, int natoms, double **refCoor, double **curCoor, double **R) {
    int    c, r, n, j, i, s, index;
    double *d = calloc(2 * ndim, sizeof(double));
    double xnr, xpc;
    double max_d;
    
    MATRIXZERO(omega, 2 * DIM, 2 * DIM);
    MATRIXZERO(om,    2 * DIM, 2 * DIM);
    MATRIXZERO(vh,    2 * DIM, 2 * DIM);
    MATRIXZERO(vk,    2 * DIM, 2 * DIM);
    MATRIXZERO(u,     2 * DIM, 2 * DIM);

    /*calculate the matrix U*/
    for (n = 0; n < nRow; n ++) {
        for (c = 0; c < ndim; c ++) {
            xpc = refCoor[n][c];
            for (r = 0; r < ndim; r ++) {
                xnr      = curCoor[n][r];
                u[c][r] += xnr * xpc * backboneMass[n];
            }
        }
    }
    
    /*construct omega*/
    /*omega is symmetric -> omega==omega' */
    for (r = 0; r < 2 * ndim; r ++) {
        for (c = 0; c <= r; c ++) {
            if (r >= ndim && c < ndim) {
                omega[r][c] = u[r - ndim][c];
                omega[c][r] = u[r - ndim][c];
            } else {
                omega[r][c] = 0;
                omega[c][r] = 0;
            }
        }
    }
    
    /*determine h and k*/
    Jacobi(omega, 2 * ndim, d, om);
    /* real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
     * int     natoms = number of rows and columns
     * real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
     * real       **v = v[0..n-1][0..n-1] contains the vectors in columns
     */
    
    index = 0; /* For the compiler only */
    
    /* Copy only the first ndim-1 eigenvectors */
    for (j = 0; j < ndim - 1; j ++) {
        max_d = -1000;
        for (i = 0; i < 2 * ndim; i ++) {
            if (d[i] > max_d) {
                max_d = d[i];
                index = i;
            }
        }
        d[index] = -10000;
        for (i = 0; i < ndim; i ++) {
            vh[j][i] = SQRT2 * om[i       ][index];
            vk[j][i] = SQRT2 * om[i + ndim][index];
        }
    }

    if (ndim == 3) {
        /* Calculate the last eigenvector as the outer-product of the first two.
         * This insures that the conformation is not mirrored and
         * prevents problems with completely flat reference structures.
         */
        CPROD(vh[0], vh[1], vh[2]);
        CPROD(vk[0], vk[1], vk[2]);
    } else if (ndim == 2) {
        /* Calculate the last eigenvector from the first one */
        vh[1][0] = -vh[0][1];
        vh[1][1] =  vh[0][0];
        vk[1][0] = -vk[0][1];
        vk[1][1] =  vk[0][0];
    }
    
    /* determine R */
    MATRIXZERO(R, 3, 3);
    for (r = 0; r < ndim; r ++) {
        for (c = 0; c < ndim; c ++) {
            for (s = 0; s < ndim; s ++) {
                R[r][c] += vk[s][r] * vh[s][c];
            }
        }
    }

    for (r = ndim; r < 3; r ++) R[r][r] = 1;
    free(d); d = NULL;
    
    return;
}

void Rotate(double **a, int i, int j, int k, int l, double tau, double s) {
    double g, h;
    g       = a[i][j];
    h       = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

void Jacobi(double **a, int n, double d[], double **v) {
    int    j, i;
    int    iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;
    
    b = (double *)calloc(n, sizeof(double));
    z = (double *)calloc(n, sizeof(double));

    for (ip = 0; ip < n; ip ++) {
        for (iq = 0; iq < n; iq ++) {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }

    for (ip = 0; ip < n; ip ++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }

    for (i = 1; i <= 50; i ++) {
        sm = 0.0;
        for (ip = 0; ip < n - 1; ip ++) {
            for (iq = ip + 1; iq < n; iq ++) {
                sm += fabs(a[ip][iq]);
            }
        }

        if (sm == 0.0) {
            free(z); z = NULL;
            free(b); b = NULL;
            return;
        }

        if (i < 4) {
            tresh = 0.2 * sm / (n * n);
        } else {
            tresh = 0.0;
        }

        for (ip = 0; ip < n - 1; ip ++) {
            for (iq = ip + 1; iq < n; iq ++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && fabs(d[ip]) + g == fabs(d[ip])
                    && fabs(d[iq]) + g == fabs(d[iq])) {
                    a[ip][iq] = 0.0;
                } else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (fabs(h) + g == fabs(h)) {
                        t = (a[ip][iq]) / h;
                    } else {
                        theta = 0.5 * h / (a[ip][iq]);
                        t     = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }

                    c         = 1.0 / sqrt(1 + t * t);
                    s         = t * c;
                    tau       = s / (1.0 + c);
                    h         = t * a[ip][iq];
                    z[ip]    -= h;
                    z[iq]    += h;
                    d[ip]    -= h;
                    d[iq]    += h;
                    a[ip][iq] = 0.0;

                    for (j = 0; j < ip; j ++) {
                        Rotate(a, j, ip, j, iq, tau, s);
                    }

                    for (j = ip + 1; j < iq; j ++) {
                        Rotate(a, ip, j, j, iq, tau, s);
                    }

                    for (j = iq + 1; j < n; j ++) {
                        Rotate(a, ip, j, iq, j, tau, s);
                    }

                    for (j = 0; j < n; j ++) {
                        Rotate(v, j, ip, j, iq, tau, s);
                    }
                }
            }
        }

        for (ip = 0; ip < n; ip ++) {
            b[ip] +=  z[ip];
            d[ip]  =  b[ip];
            z[ip]  =  0.0;
        }
    }
}


int ReadRMSDFile(FILE *RMSDInputFile, struct RMSDReadStr *thisRMSD) {
    char buffer[1024];
    int cur = 0, pos = 0;
    
repeat:
    if (fgets(buffer, sizeof(buffer), RMSDInputFile) == NULL)
        return -1;
    
    if (buffer[0] == '@' || buffer[0] == '#') {
        goto repeat;
    }
    
    sscanf(buffer, "%lf%n",
           &thisRMSD->step,
           &pos);
    
    cur += pos;
    int count = 0;
    while (sscanf(buffer + cur, "%lf%n", &thisRMSD->vRMSD[count ++], &pos) != EOF && count < 2) {
        cur += pos;
    }
    
    return 0;
}
