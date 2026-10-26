//
//  PSD.c
//  Analysis
//
//  Created by Size on 2021/7/31.
//  Copyright © 2021 Size Zheng. All rights reserved.
//

#include "Analysis.h"

#define HISSIZE 0.2 //unit A, histogram width

static int hisNum;
static long totalCount = 0, *count;
static double *perc, *prob, *psd;
static double distributionSize = HISSIZE;

static void CalPSD(void);
static int FindPositionIndex(double value, double gap, int maxN);

void PSDInfo(int id) {
    long step = 0, totalFrame;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE *PSDOutputFile, *TrjInputFile;
    
    for (int i = 1; i <= 3; i ++) {
        if (psdInfo.axis[i] < cellNum[i]) {
            psdInfo.axis[i] = cellNum[i];
        }
    }
    
    hisNum = 5.0 / HISSIZE; //max distance is 5.0 A
    count = (  long *)calloc(hisNum, sizeof(long  ));
    perc  = (double *)calloc(hisNum, sizeof(double));
    prob  = (double *)calloc(hisNum, sizeof(double));
    psd   = (double *)calloc(hisNum, sizeof(double));
    
    printf("Analyzing PSD information... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outPSD][id].name);
    PSDOutputFile = fopen(directory, "w");
    fprintf(PSDOutputFile, "# %17s %17s %17s %17s %17s\n",
            "Diameter",
            "Count",
            "Probability",
            "Accumulated",
            "Distribution");

    for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inTrj][id].name, fileList[id].list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        TrjInputFile = fopen(directory, "r");
        if (TrjInputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", buffer, path, __FILE__, __LINE__);
            printf("           You may need to execute -rPBC first to generate the required file!\n");
            exit(EXIT_FAILURE);
        }

        totalFrame = (sectInfo[sectNum].oldTime + sectInfo[sectNum].frameCount) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
            if (ReadGro(TrjInputFile)) break;
            
            PrintProcess(step);

            sect = &sectInfo[sectNum];
            fprintf(PSDOutputFile, "#frame = %.2lf\n", step * sect->outputRate);
            CalPSD();
            for (int i = 0; i < hisNum; i ++) {
                fprintf(PSDOutputFile, "%12.2lf %12li %12.4lf %12.4lf %12.4lf\n",
                        distributionSize * (i + 0.5),
                        count[i],
                        perc[i],
                        prob[i],
                        psd[i]);
            }
        }
        fclose(TrjInputFile);
    }
    fclose(PSDOutputFile);
    
    printf("\bDone!\n");
    fflush(stdout);
    return;
}


void CalPSD(void) {
    double gridSize[4] = {0,
        boxOrigDim[1] / psdInfo.axis[1],
        boxOrigDim[2] / psdInfo.axis[2],
        boxOrigDim[3] / psdInfo.axis[3]
    };
    
    int *gridNum = psdInfo.axis;
    int cellindex, pair, occupied;
    int curCell[4], cell_neighbor[4];
    double disCor[4], dis, minDis;
    double norm = -1;
    double *boxsize;
    double position[4] = {0};
    double positionshift[4] = {0};
    signed int x[27]={0,1,0,0,0,1,1,1,-1,0,0,0,-1,-1,-1,0,0,1,-1,1,-1,-1,-1,1,1,1,-1};
    signed int y[27]={0,0,1,0,1,0,1,1,0,-1,0,-1,0,-1,-1,-1,1,0,0,-1,1,-1,1,-1,1,-1,1};
    signed int z[27]={0,0,0,1,1,1,0,1,0,0,-1,-1,-1,0,-1,1,-1,-1,1,0,0,1,-1,-1,-1,1,1};
    
    LinkList();
    POINT_TO_STRUCT(boxsize, boxOrigDim);

    for (int nx = 0; nx < gridNum[1]; nx ++) {
        for (int ny = 0; ny < gridNum[2]; ny ++) {
            for (int nz = 0; nz < gridNum[3]; nz ++) {
                position[1] = (nx + 0.5) * gridSize[1];
                position[2] = (ny + 0.5) * gridSize[2];
                position[3] = (nz + 0.5) * gridSize[3];
                
                for (int i = 1; i <= 3; i ++) {
                    curCell[i] = position[i] / cellSize[i];
                    if (curCell[i] == cellNum[i])
                        curCell[i] --;
                }

                minDis = INFINIT * INFINIT;
                occupied = 0;
                for (int n = 0; n <= 26 && !occupied; n ++) {
                    cell_neighbor[3] = curCell[3] + z[n];
                    cell_neighbor[2] = curCell[2] + y[n];
                    cell_neighbor[1] = curCell[1] + x[n]; //scan the neighborhood 27 subcells, include the target subcell itself
                    
                    for (int j = 1; j <= 3; j ++) {
                        positionshift[j] = 0; //PBC position shift
                        if (cell_neighbor[j] < 0) {
                            cell_neighbor[j] = cellNum[j] - 1;
                            positionshift[j] = -1 * boxsize[j];
                        } else if (cell_neighbor[j] >= cellNum[j]) {
                            cell_neighbor[j] = 0;
                            positionshift[j] = boxsize[j];
                        }
                    }
                    
                    cellindex = cell_neighbor[3] * cellNum[1] * cellNum[2] +
                                cell_neighbor[2] * cellNum[1] +
                                cell_neighbor[1] + 1;
                    
                    pair = celllist[atomnum + cellindex];
                    while (pair != 0) {
                        DOT_PLUS(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);
                        DOT_MINUS(atom[pair].dynamic->coordinate, position, disCor);
                        DOT_MINUS(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);

                        dis = DOT_PROD(disCor, disCor);
                        dis = sqrt(dis); // - FFInfo.LJ[atom[pair].property->typeofAtom[0]].sigma / 2;
                        
                        if (dis <= 0) { //if this grid is occupied?
                            occupied = 1;
                            break;
                        }
                        if (minDis > dis) minDis = dis;
                                                
                        pair = celllist[pair];
                    }
                }
                
                if (!occupied) {
                    int pos = FindPositionIndex(minDis, distributionSize, hisNum);
                    if (pos >= 0) {
                        count[pos] ++;
                        totalCount ++;
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < hisNum; i ++) {
        for (int n = i; n < hisNum; n ++) {
            prob[i] += count[n];
        }
        prob[i] /= totalCount;
        perc[i] = (double)count[i] / totalCount;
    }
    
    for (int i = 1; i < hisNum; i ++) {
        psd[i] = -1 * (prob[i] - prob[i - 1]) / distributionSize;
        if (norm < psd[i]) norm = psd[i];
    }
    
    //normalize
    for (int i = 0; i < hisNum; i ++) {
        psd[i] /= norm;
    }
    
    return;
}


int FindPositionIndex(double value, double gap, int maxN) {
    for (int i = 0; i < maxN; i ++) {
        if ((value >= gap * i) &&
            (value <  gap * (i + 1))) {
            return i;
        }
    }
    
    return -1;
}
