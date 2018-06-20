//
//  SystemInformation.c
//  HBAnalysis
//
//  Created by Size Zheng on 4/14/15.
//  Copyright (c) 2015 Size Zheng. All rights reserved.
//

#include "Analysis.h"

void ReadObst(struct AtomStr *obst, char *type);


void SystemInformationInput(int id) {
    char directory[256];
    FILE *SysInfoFile;
    
    printf("Collecting system information...");
    fflush(stdout);

    sprintf(directory, "%s%s", path, files[inSys][id].name);
    SysInfoFile = fopen(directory,"rb");
    if (SysInfoFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inSys][id].name, path, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    fread(&atomnum, sizeof (int), 1, SysInfoFile);
    fread(&numofprotein, sizeof (int), 1, SysInfoFile);
    fread(&totalAminoAcids, sizeof (int), 1, SysInfoFile);
    
    connectionMap = (int **)calloc(atomnum + 1, sizeof(int *));
    for (int i = 0; i <= atomnum; i ++) {
        connectionMap[i] = (int *)calloc(atomnum + 1, sizeof(int));
    }
    
    aminoacid = (struct AAStr *)calloc(totalAminoAcids + 1, sizeof(struct AAStr));
    fread(aminoacid, sizeof (struct AAStr), (totalAminoAcids + 1), SysInfoFile);
    
    protein = (struct PepStr *)calloc(numofprotein + 1, sizeof(struct PepStr));
    fread(protein, sizeof (struct PepStr), (numofprotein + 1), SysInfoFile);
    
    atom = (struct AtomStr *)calloc(atomnum + 1, sizeof(struct AtomStr));
    for (int i = 0; i <= atomnum; i++) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));
        
        atom[i].property = property;
        atom[i].dynamic = dynamic;
        
        fread(atom[i].property, sizeof(struct PropertyStr), 1, SysInfoFile);
        struct ConstraintStr *thisConstr = NULL;
        struct ConstraintStr *nextConstr = atom[i].property->constr;
        
        int flag = 0;
        while (nextConstr) {
            struct ConstraintStr *newConstr = calloc(1, sizeof(struct ConstraintStr));
            fread(newConstr, sizeof(struct ConstraintStr), 1, SysInfoFile);
            
            struct StepPotenStr *thisStep = NULL;
            struct StepPotenStr *nextStep = newConstr->step;
            
            int flag2 = 0;
            while (nextStep) {
                struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
                fread(newStep, sizeof(struct StepPotenStr), 1, SysInfoFile);
                
                if (flag2 == 0) {
                    newConstr->step = newStep;
                    thisStep = newStep;
                    flag2 = 1;
                } else {
                    thisStep->next = newStep;
                    thisStep = thisStep->next;
                }
                
                nextStep = newStep->next;
            }
            
            if (flag == 0) {
                atom[i].property->constr = newConstr;
                thisConstr = newConstr;
                flag = 1;
            } else {
                thisConstr->next = newConstr;
                thisConstr = thisConstr->next;
            }
            
            nextConstr = newConstr->next;
        }
    }
    
    fread(&HBPotential, sizeof(struct HBPotentialStr), 1, SysInfoFile);
    
    for (int i = 0; i < 32; i ++) {
        for (int n = 0; n < 32; n ++) {
            fread(&potentialPairCollision[i][n], sizeof(struct ConstraintStr), 1, SysInfoFile);
            struct StepPotenStr *thisStep = NULL;
            struct StepPotenStr *nextStep = potentialPairCollision[i][n].step;
            
            int flag2 = 0;
            while (nextStep) {
                struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
                fread(newStep, sizeof(struct StepPotenStr), 1, SysInfoFile);
                
                if (flag2 == 0) {
                    potentialPairCollision[i][n].step = newStep;
                    thisStep = newStep;
                    flag2 = 1;
                } else {
                    thisStep->next = newStep;
                    thisStep = thisStep->next;
                }
                
                nextStep = newStep->next;
            }
        }
    }
    
    for (int type = 0; type < 11; type ++) {
        for (int i = 0; i < 32; i ++) {
            for (int n = 0; n < 32; n ++) {
                fread(&potentialPairHB[type][i][n], sizeof(struct ConstraintStr), 1, SysInfoFile);
                struct StepPotenStr *thisStep = NULL;
                struct StepPotenStr *nextStep = potentialPairHB[type][i][n].step;
                
                int flag2 = 0;
                while (nextStep) {
                    struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
                    fread(newStep, sizeof(struct StepPotenStr), 1, SysInfoFile);
                    
                    if (flag2 == 0) {
                        potentialPairHB[type][i][n].step = newStep;
                        thisStep = newStep;
                        flag2 = 1;
                    } else {
                        thisStep->next = newStep;
                        thisStep = thisStep->next;
                    }
                    
                    nextStep = newStep->next;
                }
            }
        }
    }
    fclose(SysInfoFile);
        
    printf("Done!\n");
    fflush(stdout);

    return;
}


int ReadLog(void) {
    char buffer[1024], sign[4];
    FILE *inputFile;
    
    for (int sectNum = 0; sectNum < fileList.count; sectNum ++) {
        sprintf(buffer, "%s%s", path, fileList.list[sectNum + 1]);
        inputFile = fopen(buffer, "r");
        if (inputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", fileList.list[sectNum + 1], path, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        struct SectionStr *sect = &sectInfo[sectNum];
        
        while (fgets(buffer, sizeof(buffer), inputFile) != NULL) {
            if (buffer[0] == '#') {
                fgets(buffer, sizeof(buffer), inputFile);   //System
                fgets(buffer, sizeof(buffer), inputFile);   //Protein Num
                fgets(buffer, sizeof(buffer), inputFile);   //Protein Seq
                
                for (int i = 0; i < numofprotein; i ++) {
                    fgets(buffer, sizeof(buffer), inputFile);
                }
                
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->oldTime);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->targetTime);

                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->temp);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->outputRate);
                
                if (sectNum == 0) {
                    fscanf(inputFile, "%s%s%lf\n", buffer, sign, &cutoffr); //cut off radius cannot change between sections
                } else {
                    fgets(buffer, sizeof(buffer), inputFile);
                }
                
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->methodType);
                
                fgets(buffer, sizeof(buffer), inputFile);   //ThermostatType
                fgets(buffer, sizeof(buffer), inputFile);   //ThermostatParameter
                fgets(buffer, sizeof(buffer), inputFile);   //DMDMethod
                fgets(buffer, sizeof(buffer), inputFile);   //ThreadNo
                
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->wallObj.wallExist);
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->wallObj.wallType);
                
                if (strcmp(sect->wallObj.wallExist, "smooth") == 0) {
                    ReadObst(&sect->wallObj.wall, "Wall");
                }
                
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->wallDyn.mark);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->wallDyn.size);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->wallDyn.rate);
                
                if (sectNum == 0) { //original box size cannot change between sections
                    if (strcmp(sect->wallDyn.mark, "no")) {
                        fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                               &boxOrigDim[1],
                               &boxOrigDim[2],
                               &boxOrigDim[3]);
                        fgets(buffer, sizeof(buffer), inputFile);
                    } else {
                        fgets(buffer, sizeof(buffer), inputFile);
                        fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                               &boxOrigDim[1],
                               &boxOrigDim[2],
                               &boxOrigDim[3]);
                    }
                } else {
                    fgets(buffer, sizeof(buffer), inputFile);
                    fgets(buffer, sizeof(buffer), inputFile);
                }
                
                fscanf(inputFile, "%s%s%i\n",  buffer, sign, &sect->flow.mark);
                if (sect->flow.mark == 1) {
                    fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                           &sect->flow.constV.v[1],
                           &sect->flow.constV.v[2],
                           &sect->flow.constV.v[3]);
                } else if (sect->flow.mark == 2) {
                    fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                           &sect->flow.force.f[1],
                           &sect->flow.force.f[2],
                           &sect->flow.force.f[3]);
                } else if (sect->flow.mark == 3) {
                    fscanf(inputFile, "%s%s%i", buffer, sign, &sect->flow.charge.num);
                    
                    DOUBLE_2CALLOC(sect->flow.charge.position, sect->flow.charge.num, 4);
                    DOUBLE_2CALLOC(sect->flow.charge.velocity, sect->flow.charge.num, 4);
                    sect->flow.charge.potGrad = (double *)calloc(sect->flow.charge.num, sizeof(double));
                    sect->flow.charge.gap     = (double *)calloc(sect->flow.charge.num, sizeof(double));
                    
                    for (int i = 0; i < sect->flow.charge.num; i ++) {
                        fscanf(inputFile, "%s%lf%lf%lf%lf%lf%s",
                               sign,
                               &sect->flow.charge.position[i][1],
                               &sect->flow.charge.position[i][2],
                               &sect->flow.charge.position[i][3],
                               &sect->flow.charge.potGrad[i],
                               &sect->flow.charge.gap[i],
                               sign);
                    }
                    fgets(buffer, sizeof(buffer), inputFile);   //the last \n
                } else {
                    fgets(buffer, sizeof(buffer), inputFile);
                }
                
                fscanf(inputFile, "%s%s%i\n", buffer, sign, &sect->obstObj.num);
                if (sect->obstObj.num) {
                    DOUBLE_2CALLOC(sect->obstObj.position, sect->obstObj.num, 4);
                    sect->obstObj.obst = (struct AtomStr *)calloc(sect->obstObj.num, sizeof(struct AtomStr));
                    fscanf(inputFile, "%s%s", buffer, sign);
                    
                    for (int i = 0; i < sect->obstObj.num; i ++) {
                        fscanf(inputFile, "%lf%lf%lf",
                               &sect->obstObj.position[i][1],
                               &sect->obstObj.position[i][2],
                               &sect->obstObj.position[i][3]);
                        if (i < sect->obstObj.num - 1) {
                            fscanf(inputFile, "%s", sign);
                        }
                        
                        ReadObst(&sect->obstObj.obst[i], "Obst");
                    }
                    fgets(buffer, sizeof(buffer), inputFile);   //the last \n
                } else {
                    fgets(buffer, sizeof(buffer), inputFile);
                }
                
                fscanf(inputFile, "%s%s%i\n", buffer, sign, &sect->tunlObj.mark);
                if (sect->tunlObj.mark) {
                    fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                           &sect->tunlObj.startPosition,
                           &sect->tunlObj.endPosition,
                           &sect->tunlObj.diameter);
                    
                    ReadObst(&sect->tunlObj.tunnel, "Tunl");                    
                } else {
                    fgets(buffer, sizeof(buffer), inputFile);
                }
            }
        }
        
        if (strcmp(sect->wallObj.wallExist, "smooth") == 0) {
            sect->wallObj.wall.dynamic->coordinate[1] = 0.5 * boxOrigDim[1];
            sect->wallObj.wall.dynamic->coordinate[2] = 0.5 * boxOrigDim[2];
            sect->wallObj.wall.dynamic->coordinate[3] = 0.5 * boxOrigDim[3];
            
            if (strcmp(sect->wallObj.wallType, "parallel") == 0) {
                sect->wallObj.wall.dynamic->coordinate[1] =
                sect->wallObj.wall.dynamic->coordinate[3] = 0;
            } else if (strcmp(sect->wallObj.wallType, "cylinder") == 0) {
                sect->wallObj.wall.dynamic->coordinate[1] = 0;
            }
        }
        
        if (sect->tunlObj.mark) {
            sect->tunlObj.tunnel.dynamic->coordinate[1] = 0;
            sect->tunlObj.tunnel.dynamic->coordinate[2] = 0.5 * boxOrigDim[2];
            sect->tunlObj.tunnel.dynamic->coordinate[3] = 0.5 * boxOrigDim[3];
        }
        
        fclose(inputFile);
    }
    
    for (int sectNum = 0; sectNum < fileList.count; sectNum ++) {
        if (sectNum < fileList.count - 1) {
            sectInfo[sectNum].frameCount = sectInfo[sectNum + 1].oldTime - sectInfo[sectNum].oldTime;
        } else {
            sectInfo[sectNum].frameCount = sectInfo[sectNum].targetTime - sectInfo[sectNum].oldTime;
        }
    }
    
    return 1;
}


void ReadObst(struct AtomStr *obst, char *type) {
    char rName[16];
    char directory[1024], buffer[1024];
    FILE *inputFile;
    
    sprintf(directory, "%s", obstDir);
    inputFile = fopen(directory, "r");
    
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find the obstruction info file at %s!\n", directory);
        printf("           you may need to use flag -obs to provide the path to the obstruction info file!\n");
        exit(EXIT_FAILURE);
    }
    
    if (obst->property == NULL) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        obst->property = property;
    }
    
    if (obst->dynamic == NULL) {
        struct DynamicStr  *dynamic  = calloc(1, sizeof(struct DynamicStr));
        obst->dynamic  = dynamic;
    }
    
    fgets(buffer, sizeof(buffer), inputFile);
    fscanf(inputFile, "%s%s%s%lf%s%s", buffer, obst->property->name, rName, &obst->property->mass,
           obst->property->extraProperty[0],
           obst->property->extraProperty[1]);
    
    obst->property->type = AtomModel(rName);
    sprintf(obst->property->nameOfAA, "%s", type);
    
    obst->dynamic->event.partner      = INVALID;
    obst->dynamic->event.subEventType = Invalid;
    obst->dynamic->event.eventType    = Invd_Event;
    
    fclose(inputFile);
    return;
}


void FreeVariables(void) {
    free(aminoacid);
    free(protein);
    
    for (int i = 0; i <= atomnum; i ++) {
        free(connectionMap[i]);
        free(atom[i].dynamic);
        free(atom[i].property);
    }
    free(connectionMap);
    free(atom);
    free(celllist);

    for (int i = 0; i < NumofFileType; i++) {
        free(files[i]);
    }
    free(files);
    free(analysisList);
    
    return;
}
