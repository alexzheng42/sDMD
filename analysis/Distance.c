//
//  Distance.c
//  Analysis
//
//  Created by Size Zheng on 05/25/20.
//  Copyright © 2020 Size Zheng. All rights reserved.
//

#include "Analysis.h"
#include <stdio.h>

struct DistStr {
    double dist;
    double eng;
};

static void CalParticle2Wall(char *style, struct SectionStr *sect, struct DistStr **dist, struct GroupStr group);
static void CalCGSurEnergy(struct SectionStr *sect, struct DistStr **dist, struct GroupStr group);


void DisInfo(int id) {
    long step = 0, totalFrame;
    char directory[1024], buffer[1024];
    struct DistStr **dist;
    struct SectionStr *sect;
    FILE *TrjInputFile = NULL;
    FILE **CG2ObstOutputFile = NULL;
    
    if (strcmp(disInfo.type, "CG2Wall") == 0) {
        CG2ObstOutputFile = calloc(sectInfo[0].obstObj.num, sizeof(FILE*));
        if (sectInfo[0].obstObj.num == 1) {
            sprintf(directory, "%s%s", path, files[outCG2Obst][id].name);
            CG2ObstOutputFile[0] = fopen(directory, "w");
            fprintf(CG2ObstOutputFile[0], "# TotalAA\t%i\n", numofRes);
            fprintf(CG2ObstOutputFile[0], "# TotalMol\t%i\n", numofMol);
            fprintf(CG2ObstOutputFile[0], "# AASeq\t");
            for (int i = 1; i <= numofMol; i ++) {
                fprintf(CG2ObstOutputFile[0], "%i\t", molecule[i].endAANum - molecule[i].startAANum + 1);
            }
            fprintf(CG2ObstOutputFile[0], "\n");
            fprintf(CG2ObstOutputFile[0], "# %s\t", "Step");
            for (int i = 1; i <= numofRes; i ++) {
                fprintf(CG2ObstOutputFile[0], "\tAA%i", i);
            }
            fprintf(CG2ObstOutputFile[0], "\n");
            fprintf(CG2ObstOutputFile[0], "# %s\t", "Step");
            for (int i = 1; i <= numofRes; i ++) {
                fprintf(CG2ObstOutputFile[0], "\t%s", residue[i].nameofAA);
            }
            fprintf(CG2ObstOutputFile[0], "\n");
        } else {
            for (int i = 0; i < sectInfo[0].obstObj.num; i ++) {
                sprintf(directory, "%s%s_%i", path, files[outCG2Obst][id].name, (i + 1));
                CG2ObstOutputFile[i] = fopen(directory, "w");
                fprintf(CG2ObstOutputFile[0], "# TotalAA %i\n", numofRes);
                fprintf(CG2ObstOutputFile[0], "# TotalMol\t%i\n", numofMol);
                fprintf(CG2ObstOutputFile[0], "# AASeq \t");
                for (int i = 1; i <= numofMol; i ++) {
                    fprintf(CG2ObstOutputFile[0], "%i\t", molecule[i].endAANum - molecule[i].startAANum + 1);
                }
                fprintf(CG2ObstOutputFile[0], "\n");
                fprintf(CG2ObstOutputFile[i], "# %s\t", "Step");
                for (int n = 1; n <= numofRes; n ++) {
                    fprintf(CG2ObstOutputFile[i], "\tAA%i", n);
                }
                fprintf(CG2ObstOutputFile[i], "\n");
                fprintf(CG2ObstOutputFile[0], "# %s\t", "Step");
                for (int i = 1; i <= numofRes; i ++) {
                    fprintf(CG2ObstOutputFile[0], "\t%s", residue[i].nameofAA);
                }
                fprintf(CG2ObstOutputFile[0], "\n");
            }
        }
        
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
        dist = (struct DistStr **)calloc(numofRes + 1, sizeof(struct DistStr *));
        for (int i = 0; i <= numofRes; i ++) {
            dist[i] = (struct DistStr *)calloc(sect->obstObj.num, sizeof(struct DistStr));
        }

        totalFrame = (sectInfo[sectNum].oldTime + sectInfo[sectNum].frameCount) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
            if (ReadGro(TrjInputFile)) break;
            
            PrintProcess(step);
            if (strcmp(disInfo.type, "CG2Wall") == 0) {
                CalParticle2Wall("alphaC", sect, dist, group);
                CalCGSurEnergy(sect, dist, group);

                for (int obstNum = 0; obstNum < sect->obstObj.num; obstNum ++) {
                    fprintf(CG2ObstOutputFile[obstNum], "%10.2lf ", step * sect->outputRate);
                    for (int AANum = 1; AANum <= numofRes; AANum ++) {
                        fprintf(CG2ObstOutputFile[obstNum], "%8.2lf %8.2lf", dist[AANum][obstNum].dist, dist[AANum][obstNum].eng);
                    }
                    fprintf(CG2ObstOutputFile[obstNum], "\n");
                }
            }
        }
        fclose(TrjInputFile);
        for (int i = 0; i <= numofRes; i ++) {
            free(dist[i]);
        }
        free(dist);
    }

    for (int i = 0; i < sectInfo[0].obstObj.num; i ++) {
        fclose(CG2ObstOutputFile[i]);
    }

    printf("\bDone!\n");
    fflush(stdout);
    return;
}


void CalParticle2Wall(char *style, struct SectionStr *sect, struct DistStr **dist, struct GroupStr group) {
    int currentAA = 0;
    double position_i[4] = {0}, position_j[4] = {0};
    double r_ij[4];

    if (strcmp(style, "alphaC") == 0) {
        for (int n = 0; n < sect->obstObj.num; n ++) {
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
            
                for (int m = 1; m <= group.atomNum; m ++) {
                    int num = group.list_Atom[m];
                    if (FindStrInList(atom[num].property->extraProperty, "CG", 0) &&
		            	strcmp(atom[num].property->name, "CA") == 0) {

                        struct AtomStr *target = &atom[num];
                        double distance2 = 0;

                        TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
                        for (int i = 1; i <= 3; i ++) {
                            if (i != dim) position_i[i] = 0;
                        }

                        DOT_MINUS(position_i, position_j, r_ij);
                        distance2 = DOT_PROD(r_ij, r_ij);

                        currentAA = molecule[target->property->sequence.molNum - 1].endAANum +
                                    target->property->sequence.resNum;
                        dist[currentAA][n].dist = sqrt(distance2);
                    }
                }
            }
        }
    }

    return;
}


void CalCGSurEnergy(struct SectionStr *sect, struct DistStr **dist, struct GroupStr group) {
    int atomType, surfType;
    int currentAA = 0;
	double tolerance, floatzero = 10E-8;
    double position_i[4], position_j[4];
	double speed_i[4], speed_j[4] = { 0 };
    double r_ij[4], v_ij[4], direction;
    
    for (int n = 0; n < sect->obstObj.num; n ++) {
        struct AtomStr *thisSurface = &sect->obstObj.obst[n];
        for (int dim = 1; dim <= 3; dim ++) {
            if (sect->obstObj.position[n][dim] < 0) {
                continue;
            }
            
            TRANSFER_VECTOR(position_j, sect->obstObj.position[n]);
            for (int i = 1; i <= 3; i ++) {
                if (i != dim) {
                    position_j[i] = 0;
                }
            }
            
            for (int m = 1; m <= group.atomNum; m ++) {
                int num = group.list_Atom[m];
                if (strcmp(atom[num].property->extraProperty[0], "CG") == 0 &&
					strcmp(atom[num].property->name,             "CA") == 0) {
                    
                    struct AtomStr *target = &atom[num];
                    double distance2 = 0;
                    double TCGSurE = 0;
                    
                    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
					TRANSFER_VECTOR(   speed_i, target->dynamic->velocity);
                    for (int i = 1; i <= 3; i ++) {
                        if (i != dim) {
                            position_i[i] = 0;
							   speed_i[i] = 0;
                        }
                    }
                    
                    DOT_MINUS(position_i, position_j, r_ij);
					DOT_MINUS(   speed_i,    speed_j, v_ij);
                    distance2 = DOT_PROD(r_ij, r_ij);
					direction = DOT_PROD(r_ij, v_ij);

					if (direction < 0)
						direction = 1;  // -> <-
					else
						direction = -1; // <- ->
					tolerance = direction * floatzero;
                    
                    atomType = target->property->typeofRes;
                    surfType = AAModel(thisSurface->property->resName);
                    struct ConstraintStr *thisConstr = &potentialPairCG[surfType][atomType];
                    
                    struct StepPotenStr *thisStep = thisConstr->step;                    
                    while (thisStep != NULL && thisStep->d != 0) {
                        if (distance2 < thisStep->d + tolerance) {
                            TCGSurE += thisStep->accumulated;
                            break;
                        }
                        thisStep = thisStep->next;
                    }
                    
                    currentAA = molecule[target->property->sequence.molNum - 1].endAANum +
                                target->property->sequence.resNum;
                    dist[currentAA][n].eng = TCGSurE;
                }
            }
        }
    }
    return;
}