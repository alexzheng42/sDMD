//
//  Cluster.c
//  Analysis
//
//  Created by Size Zheng on 11/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"

#define gapTime 200

static int *proteinList;
static int **cluster;
static long **interaction;

static int CheckInteraction (struct AtomStr *atom1, struct AtomStr *atom2, int **connection);
static int CheckListandAdd(int *list, int listLen, int num);
static void CalCluster(int id);
static void SetListZero(int *list, int listLen);
static void SetClusterZero(int **cluster);
static void InitializeVariables(void);
static void FreeVar(void);
void PrintInteraction(long **inter);
void PrintCluster(int **cluster);


void ClusterInfo(int id) {
    printf("Analyzing aggregation information... ");
    fflush(stdout);
    
    InitializeVariables();
    CalCluster(id);
    FreeVar();
    
    printf("\bDone\n");
    fflush(stdout);
}

static void InitializeVariables(void) {
    proteinList = (int *)calloc(numofprotein + 1, sizeof(int));
    LONG_2CALLOC(interaction, (atomnum + 1), (atomnum + 1));
    INT_2CALLOC(cluster, (numofprotein + 1), (numofprotein + 1));
}

static void FreeVar() {
    free(proteinList); proteinList = NULL;
    _2FREE(interaction);
    _2FREE(cluster);
}

void CalCluster(int id) {
    long step = 0, totalFrame;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE *inputTrjFile, *inputCntFile;
    FILE *outputAggNumFile;
    
    sprintf(directory, "%s%s", path, files[outAgg][id].name);
    outputAggNumFile = fopen(directory, "w");
    
    sprintf(directory, "%s%s", path, files[outTrj][id].name);
    inputTrjFile = fopen(directory, "r");
    if (inputTrjFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outTrj][id].name, path, __FILE__, __LINE__);
        printf("           You may need to execute -rPBC first to generate the required file!\n");
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
            if (ReadGro(inputTrjFile) || ReadConnectionMap(inputCntFile)) break;
            
            PrintProcess(step);
            
            sect = &sectInfo[sectNum];
            for (int i1 = 1; i1 <= numofprotein; i1 ++) {
                for (int n1 = protein[i1].startAtomNum; n1 <= protein[i1].endAtomNum; n1 ++) {
                    
                    for (int i2 = i1 + 1; i2 <= numofprotein; i2 ++) {
                        for (int n2 = protein[i2].startAtomNum; n2 <= protein[i2].endAtomNum; n2 ++) {
                            if (CheckInteraction(&atom[n1], &atom[n2], connectionMap)) {
                                if (interaction[n1][n2] == -1) {
                                    interaction[n1][n2] = step;
                                }
                            } else {
                                interaction[n1][n2] = -1;
                            }
                        }
                    }
                }
            }
            
            if (step * sect->outputRate > gapTime) {
                SetClusterZero(cluster);
                SetListZero(proteinList, numofprotein);
                cluster[0][0] = numofprotein;
                for (int i1 = 1; i1 <= numofprotein; i1 ++) {
                    for (int n1 = protein[i1].startAtomNum; n1 <= protein[i1].endAtomNum; n1 ++) {
                        
                        for (int i2 = i1 + 1; i2 <= numofprotein; i2 ++) {
                            if (cluster[i1][i2] == 0) {
                                for (int n2 = protein[i2].startAtomNum; n2 <= protein[i2].endAtomNum; n2 ++) {
                                    if (interaction[n1][n2] > 0 && (step - interaction[n1][n2] > gapTime / sect->outputRate)) {
                                        if (!CheckListandAdd(proteinList, numofprotein, i2)) {
                                            cluster[0][0] --;
                                            cluster[i1][i2] = 1;
                                            cluster[i2][i1] = 1;
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                fprintf(outputAggNumFile, "%10.2lf%20i\n", (double)(step * sect->outputRate), cluster[0][0]);
            }
        }
        fclose(inputCntFile);
    }
    
    fclose(inputTrjFile);
    fclose(outputAggNumFile);
    
    return;
}

int CheckInteraction (struct AtomStr *atom1, struct AtomStr *atom2, int **connection) {
    int target_i = atom1->property->num;
    int target_j = atom2->property->num;
    double position_i[4], position_j[4];
    double r_ij[4], distance2;
    double range = cutoffr * cutoffr;
    
    TRANSFER_VECTOR(position_i, atom1->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, atom2->dynamic->coordinate);
    DOT_MINUS(position_i, position_j, r_ij);
    distance2 = DOT_PROD(r_ij, r_ij);
    
    if (distance2 < range || (connection[target_i][target_j] & HBConnect)) {
        return 1;
    }
    
    return 0;
}


void SetListZero(int *list, int listLen) {
    for (int i = 0; i < listLen; i ++) {
        list[i] = 0;
    }
}


void SetClusterZero(int **cluster) {
    for (int i = 0; i <= numofprotein; i ++) {
        for (int n = 0; n <= numofprotein; n ++) {
            cluster[i][n] = 0;
        }
    }
}


int CheckListandAdd(int *list, int listLen, int num) {
    int i;
    
    for (i = 0; i < listLen; i ++) {
        if (list[i] == 0) {
            break;
        }
        
        if (list[i] == num) {
            return 1;
        }
    }
    
    if (i == listLen) {
        printf("!!ERROR!!: protein list is overflow! %s:%i\n", __FILE__, __LINE__);
    }
    
    list[i] = num;
    return 0;
}

void PrintInteraction(long **inter) {
    for (int i = 1; i <= atomnum; i ++) {
        if (i <= atomnum / 2) {
            for (int n = 1; n <= atomnum; n ++) {
                if (inter[i][n] > 0) {
                    printf("atom %i : %i, protein %i : %i\n", i, n, atom[i].property->sequence.proteinNum, atom[n].property->sequence.proteinNum);
                }
            }
        }
    }
    
    return;
}

void PrintCluster(int **cluster) {
    for (int i = 1; i <= numofprotein; i ++) {
        for (int n = 1; n <= numofprotein; n ++) {
            printf("%i ", cluster[i][n]);
        }
        printf("\n");
    }
}

