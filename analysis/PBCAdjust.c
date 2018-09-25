//
//  PBCAdjust.c
//  DataAnalysis
//
//  Created by Size Zheng on 4/16/15.
//  Copyright (c) 2015 Size Zheng. All rights reserved.
//

#include "Analysis.h"

void SaveGRO(FILE *saveFile);
void SaveGROMono(FILE *saveFile);


void AdjustPBC(int id)
{
    int dim = 3, step = 1;
    long frameCount, totalFrame;
    double gap[4];
    char directory[1024], buffer[1024];
    FILE *positionInput;
    FILE *positionOutput;
	FILE *positionOutputMono = NULL;
    
    printf("Checking trajectory file... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outTrj][id].name);
    positionOutput = fopen(directory, "w");

	if (nPP) {
		sprintf(directory, "%srPBCGRO_%i.gro", path, nPP);
		positionOutputMono = fopen(directory, "w");
	}
    
    for (int sectNum = 0; sectNum < fileList.count; sectNum ++) {
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inTrj][id].name, fileList.list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        positionInput = fopen(directory, "r");
        if (positionInput == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inTrj][id].name, path, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        if (strcmp(sectInfo[sectNum].wallObj.wallExist, "no")) {
            if (strncmp(sectInfo[sectNum].wallObj.wallType, "sphere", 1) == 0)
                dim = 0;
            if (strncmp(sectInfo[sectNum].wallObj.wallType, "parallel", 1) == 0)
                step = 2;
            else if (strncmp(sectInfo[sectNum].wallObj.wallType, "cylinder", 1) == 0)
                dim = 1;
        }
        
        totalFrame = (sectInfo[sectNum].oldTime + sectInfo[sectNum].frameCount) / sectInfo[sectNum].outputRate;
        for (frameCount = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; frameCount < totalFrame; frameCount ++) {
            if (ReadGro(positionInput)) break;
            
            PrintProcess(frameCount);
            
            gap[1] = boxCurtDim[1] / 2;
            gap[2] = boxCurtDim[2] / 2;
            gap[3] = boxCurtDim[3] / 2;
            
            for (int j = 1; j <= numofprotein; j ++) {
                int target = (protein[j].startAtomNum + protein[j].endAtomNum) * 0.5;
                
                for (int i = target + 1; i <= protein[j].endAtomNum; i ++) {
                    for (int n = 1; n <= dim; n += step) {
                        if (atom[i].dynamic->coordinate[n] - atom[i - 1].dynamic->coordinate[n] > gap[n]) {
                            atom[i].dynamic->coordinate[n] -= boxCurtDim[n];
                        } else if (atom[i - 1].dynamic->coordinate[n] - atom[i].dynamic->coordinate[n] > gap[n]) {
                            atom[i].dynamic->coordinate[n] += boxCurtDim[n];
                        }
                    }
                }
                
                for (int i = target - 1; i >= protein[j].startAtomNum; i --) {
                    for (int n = 1; n <= dim; n += step) {
                        if (atom[i].dynamic->coordinate[n] - atom[i + 1].dynamic->coordinate[n] > gap[n]) {
                            atom[i].dynamic->coordinate[n] -= boxCurtDim[n];
                        } else if (atom[i + 1].dynamic->coordinate[n] - atom[i].dynamic->coordinate[n] > gap[n]) {
                            atom[i].dynamic->coordinate[n] += boxCurtDim[n];
                        }
                    }
                }
            }
            SaveGRO(positionOutput);

			if (nPP) {
				SaveGROMono(positionOutputMono);
			}
        }
        fclose(positionInput);
    }
    
    fclose(positionOutput);
    
    printf("\bDone!\n");
    fflush(stdout);
    
    return;
}


int ReadGro(FILE *inputFile) { //per frame
    int num;
    char buffer[1024];
    
    if (fscanf(inputFile, "%[^\n]\n", buffer) == EOF)
        return -1;
    
    fscanf(inputFile, "%i", &num);
    
    if (num != atomnum) {
        printf("!!ERROR!!: atom number does not match! from info file: %i, from trajectory: %i! %s:%i\n", atomnum, num, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    for (int i = 1; i <= num; i ++) {
        fscanf(inputFile, "%s%s%s%lf%lf%lf%lf%lf%lf", buffer, buffer, buffer,
               &atom[i].dynamic->coordinate[1],
               &atom[i].dynamic->coordinate[2],
               &atom[i].dynamic->coordinate[3],
               &atom[i].dynamic->velocity[1],
               &atom[i].dynamic->velocity[2],
               &atom[i].dynamic->velocity[3]);
        atom[i].dynamic->coordinate[1] *= 10;
        atom[i].dynamic->coordinate[2] *= 10;
        atom[i].dynamic->coordinate[3] *= 10;
    }
    fscanf(inputFile, "\n");
    fscanf(inputFile, "%lf%lf%lf\n",
           &boxCurtDim[1],
           &boxCurtDim[2],
           &boxCurtDim[3]);
    boxCurtDim[1] *= 10;
    boxCurtDim[2] *= 10;
    boxCurtDim[3] *= 10;
    
    return 0;
}


void SaveGROMono(FILE *saveFile) {
    fprintf(saveFile, "Protein %i\n", nPP);
    fprintf(saveFile, "%5i\n", protein[nPP].endAtomNum - protein[nPP].startAtomNum + 1);
    
    for (int i = protein[nPP].startAtomNum; i <= protein[nPP].endAtomNum; i++) {
        fprintf(saveFile, "%5i%-5s%5s%5i%8.3f%8.3f%8.3f%10.4f%10.4f%10.4f\n",
                atom[i].property->sequence.aminoacidNum, atom[i].property->nameOfAA, atom[i].property->name, i,
                atom[i].dynamic->coordinate[1] / 10,
                atom[i].dynamic->coordinate[2] / 10,
                atom[i].dynamic->coordinate[3] / 10,
                atom[i].dynamic->velocity[1],
                atom[i].dynamic->velocity[2],
                atom[i].dynamic->velocity[3]);
    }
    
    fprintf(saveFile, "%10.5f%10.5f%10.5f\n",
            boxCurtDim[1] / 10,
            boxCurtDim[2] / 10,
            boxCurtDim[3] / 10);
    fflush(saveFile);
}


void SaveGRO(FILE *saveFile) {
	fprintf(saveFile, "model\n");
	fprintf(saveFile, "%5i\n", atomnum);

	for (int i = 1; i <= atomnum; i++) {
		fprintf(saveFile, "%5i%-5s%5s%5i%8.3f%8.3f%8.3f%10.4f%10.4f%10.4f\n",
			atom[i].property->sequence.aminoacidNum, atom[i].property->nameOfAA, atom[i].property->name, i,
			atom[i].dynamic->coordinate[1] / 10,
			atom[i].dynamic->coordinate[2] / 10,
			atom[i].dynamic->coordinate[3] / 10,
			atom[i].dynamic->velocity[1],
			atom[i].dynamic->velocity[2],
			atom[i].dynamic->velocity[3]);
	}

	fprintf(saveFile, "%10.5f%10.5f%10.5f\n",
		boxCurtDim[1] / 10,
		boxCurtDim[2] / 10,
		boxCurtDim[3] / 10);
	fflush(saveFile);
}

