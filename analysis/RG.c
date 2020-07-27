//
//  RG.c
//  Analysis
//
//  Created by Size on 2018/12/19.
//  Copyright © 2018 Size Zheng. All rights reserved.
//

#include "Analysis.h"

static int beginP, endP, offSet, thisAtomNum;
static double mass;

static void InitializeVariables(void);
static void CalRGMono(int id);
static void FindCOMSingle(double *COM);
static void FreeVar(void);

void RGInfo(int id) {
    printf("Analyzing radius of gyration... ");
    fflush(stdout);
    
    InitializeVariables();
    CalRGMono(id);
    FreeVar();
    
    printf("\bDone\n");
    fflush(stdout);
}

void InitializeVariables(void) {
    beginP = 1;
    endP = numofprotein;
    thisAtomNum = atomnum;
    offSet = 0;
    mass = 0;
    
    if (nPP) {
        beginP = endP = nPP;
        thisAtomNum = protein[nPP].endAtomNum - protein[nPP].startAtomNum + 1;
        offSet = protein[nPP].startAtomNum - 1;
    }
    
    for (int n = 1 + offSet; n <= thisAtomNum + offSet; n ++) {
        mass += atom[n].property->mass;
    }
}

void FreeVar(void) {
}

void CalRGMono(int id) //IUPAC Definition
{
    double RG, COM[4] = {0};
    long step = 0, totalFrame;
    char directory[1024];
    struct SectionStr *sect;
    FILE *inputTrjFile;
    FILE *outputRGFile;
    
    sprintf(directory, "%s%s", path, files[outRG][id].name);
    outputRGFile = fopen(directory, "w");
    fprintf(outputRGFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(outputRGFile, "%s\n", "@    yaxis  label \"Rg (nm)\"\n");
    
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

            RG = 0;
            sect = &sectInfo[sectNum];
            FindCOMSingle(COM);
            for (int n = 1 + offSet; n <= thisAtomNum + offSet; n ++) {
                for (int i = 1; i <= 3; i ++) {
                    RG += (COM[i] - atom[n].dynamic->coordinate[i]) * (COM[i] - atom[n].dynamic->coordinate[i]) * atom[n].property->mass;
                }
            }

            RG /= mass;
            RG = sqrt(RG);

            fprintf(outputRGFile, "%15.2lf %15.4lf\n", step * sect->outputRate, RG / 10);
        }
    }

    fclose(inputTrjFile);
    fclose(outputRGFile);
    return;
}

void FindCOMSingle(double *COM) {
    COM[1] = COM[2] = COM[3] = 0;

    for (int n = 1 + offSet; n <= thisAtomNum + offSet; n ++) {
        for (int i = 1; i <= 3; i ++) {
            COM[i] += atom[n].dynamic->coordinate[i] * atom[n].property->mass;
        }
    }

    for (int i = 1; i <= 3; i ++) {
        COM[i] /= mass;
    }

    return;
}
