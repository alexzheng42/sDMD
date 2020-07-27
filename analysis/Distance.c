//
//  Distance.c
//  Analysis
//
//  Created by Size Zheng on 05/25/20.
//  Copyright © 2020 Size Zheng. All rights reserved.
//

#include "Analysis.h"

static int beginP, endP, offSet, thisAtomNum;

void DisInfo(int id) {
    long step = 0, totalFrame;
    char directory[1024], buffer[1024];
    double **dist;
    struct SectionStr *sect;
    FILE *TrjInputFile;
    FILE *CG2ObstOutputFile;
    
    beginP = 1;
    endP = numofprotein;
    thisAtomNum = atomnum;
    offSet = 0;

    if (nPP) {
        beginP = endP = nPP;
        thisAtomNum = protein[nPP].endAtomNum - protein[nPP].startAtomNum + 1;
        offSet = protein[nPP].startAtomNum - 1;
    }

    if (strcmp(disInfo.type, "CG2Wall") == 0) {
        sprintf(directory, "%s%s", path, files[outCG2Obst][id].name);
        CG2ObstOutputFile = fopen(directory, "w");
        fprintf(CG2ObstOutputFile, "# %s\t", "Step");
        for (int i = 1; i <= totalAminoAcids; i ++) {
            fprintf(CG2ObstOutputFile, "\tAA%i", i);
        }
        fprintf(CG2ObstOutputFile, "\n");
    }
    
    printf("Analyzing distancing... ");
    fflush(stdout);
    
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

        sect = &sectInfo[sectNum];
        DOUBLE_2CALLOC(dist, (totalAminoAcids + 1), sect->obstObj.num);
        
        totalFrame = (sectInfo[sectNum].oldTime + sectInfo[sectNum].frameCount) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
            if (ReadGro(TrjInputFile)) break;
            
            PrintProcess(step);
            if (strcmp(disInfo.type, "CG2Wall") == 0) {
                CalParticle2Wall("alphaC", sect, dist);

                fprintf(CG2ObstOutputFile, "%10.2lf ", step * sect->outputRate);
                for (int AANum = 1; AANum <= totalAminoAcids; AANum ++) {
                    for (int obstNum = 0; obstNum < sect->obstObj.num; obstNum ++) {
                        fprintf(CG2ObstOutputFile, "%8.2lf ", dist[AANum][obstNum]);
                    }
                    fprintf(CG2ObstOutputFile, "  ");
                }
                fprintf(CG2ObstOutputFile, "\n");
            }
        }
        fclose(TrjInputFile);
        _2FREE(dist);
    }
    fclose(CG2ObstOutputFile);

    printf("\bDone!\n");
    fflush(stdout);
    return;
}

void CalParticle2Wall(char *style, struct SectionStr *sect, double **dist) {
    int currentAA = 0;
    double position_i[4] = {0}, position_j[4] = {0};
    double r_ij[4];

    if (strcmp(style, "alphaC") == 0) {
        for (int n = 0; n < sect->obstObj.num; n ++) {
            struct AtomStr *thisSurface = &sect->obstObj.obst[n];
            for (int dim = 1; dim <= 3; dim ++) { //only for surface on xy, xz or yz
                if (sect->obstObj.position[n][dim] < 0) {
                    continue;
                }

                TRANSFER_VECTOR(position_j, sect->obstObj.position[n]);
                for (int i = 1; i <= 3; i ++) {
                    if (i != dim) {
                        position_j[i] = 0;
                    }
                }
            
                for (int num = 1; num <= atomnum; num ++) {
                    if (strcmp(atom[num].property->extraProperty[0], "CG") == 0 &&
		            	strcmp(atom[num].property->name,             "CA") == 0) {

                        struct AtomStr *target = &atom[num];
                        double distance2 = 0;

                        TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
                        for (int i = 1; i <= 3; i ++) {
                            if (i != dim) position_i[i] = 0;
                        }

                        DOT_MINUS(position_i, position_j, r_ij);
                        distance2 = DOT_PROD(r_ij, r_ij);

                        currentAA = protein[target->property->sequence.proteinNum - 1].endAANum +
                                    target->property->sequence.aminoacidNum;
                        dist[currentAA][n] = sqrt(distance2);
                    }
                }
            }
        }
    }

    return;
}