//
//  Density.c
//  Analysis
//
//  Created by Size Zheng on 07/22/21.
//  Copyright © 2021 Size Zheng. All rights reserved.
//

#include "Analysis.h"

static double CalDensitySlice(int dim, double start, double end);

void DensityInfo(int id) {
    int dim = densInfo.dim;
    long step = 0, totalFrame;
    double pos, gap = boxOrigDim[densInfo.dim] / densInfo.pointNum;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE *densityOutputFile, *TrjInputFile;

    printf("Analyzing density information... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outDensity][id].name);
    densityOutputFile = fopen(directory, "w");

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
            fprintf(densityOutputFile, "#frame = %.2lf\n", step * sect->outputRate);
            pos = 0;
            while (pos + gap <= boxOrigDim[dim]) {
                fprintf(densityOutputFile, "%10.2lf %10.2lf\n",
                        pos + gap * 0.5,
                        CalDensitySlice(dim, pos, pos + gap));
                pos += gap;
            }
        }
        fclose(TrjInputFile);
    }
    fclose(densityOutputFile);

    printf("\bDone!\n");
    fflush(stdout);
    return;
}


double CalDensitySlice(int dim, double start, double end) {
    double mass = 0;
    double vol;

    for (int i = 1; i <= atomnum; i ++) {
        if (atom[i].dynamic->coordinate[dim] >= start &&
            atom[i].dynamic->coordinate[dim] <  end) {
            mass += atom[i].property->mass;
        }
    }

    vol = end - start;
    for (int i = 1; i <= 3; i ++) {
        if (i != dim) {
            vol *= boxOrigDim[i];
        }
    }

    return mass / vol;
}
