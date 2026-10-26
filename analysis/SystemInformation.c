//
//  SystemInformation.c
//  HBAnalysis
//
//  Created by Size Zheng on 4/14/15.
//  Copyright (c) 2015 Size Zheng. All rights reserved.
//

#include "Analysis.h"


// SysInfo.dat written with raw structs (sDMD before the counted format): the stored list
// pointers only tell whether another constraint / step follows
static void ReadSysInfoLegacy(FILE *SysInfoFile, int nType) {
    int num[2] = {0, nType};
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
        
        fread(&atom[i].dynamic->HB, sizeof(struct HBStr), 1, SysInfoFile);
    }
    
    struct { int typeNum; double scale; double barrier; } oldHB;   // one scale and barrier for all HBs
    fread(&oldHB, sizeof(oldHB), 1, SysInfoFile);
    HBPotential.typeNum = oldHB.typeNum;
    HBPotential.classOf = (int *)calloc(HBPotential.typeNum + 1, sizeof(int));
    for (int c = 0; c < HB_CLASSNUM; c++) {
        HBPotential.scale[c] = oldHB.scale;
        HBPotential.barrier[c] = oldHB.barrier;
    }
    
    _CALLOC2(struct ConstraintStr, potentialPairCollision, num[1], num[1]);
    for (int i = 0; i < num[1]; i ++) {
        for (int n = 0; n < num[1]; n ++) {
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
    
    potentialPairHB = (struct ConstraintStr ***)calloc(HBPotential.typeNum + 1, sizeof(struct ConstraintStr **));
    for (int i = 0; i <= HBPotential.typeNum; i ++) {
        _CALLOC2(struct ConstraintStr, potentialPairHB[i], num[1], num[1]);
    }
    for (int type = 0; type <= HBPotential.typeNum; type ++) {
        for (int i = 0; i < num[1]; i ++) {
            for (int n = 0; n < num[1]; n ++) {
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

}


// one constraint / potential table in the counted format: fields, number of steps, steps
static void ReadConstraint(struct ConstraintStr *c, FILE *in) {
    fread(&c->connection, sizeof(int), 1, in);
    fread(&c->dmin, sizeof(double), 1, in);
    fread(&c->dmax, sizeof(double), 1, in);
    fread(&c->e, sizeof(double), 1, in);
    fread(&c->accumulated, sizeof(double), 1, in);
    c->step = NULL;
    c->next = NULL;
    int nStep = 0;
    fread(&nStep, sizeof(int), 1, in);
    struct StepPotenStr **tail = &c->step;
    for (int k = 0; k < nStep; k++) {
        struct StepPotenStr *s = calloc(1, sizeof(struct StepPotenStr));
        fread(&s->d, sizeof(double), 1, in);
        fread(&s->e, sizeof(double), 1, in);
        fread(&s->accumulated, sizeof(double), 1, in);
        *tail = s;
        tail = &s->next;
    }
}


// SysInfo.dat in the counted format (sDMD SaveSysInfo, SYSINFO_VERSION 2)
static void ReadSysInfoCounted(FILE *SysInfoFile, int nType) {
    atom = (struct AtomStr *)calloc(atomnum + 1, sizeof(struct AtomStr));
    for (int i = 0; i <= atomnum; i++) {
        struct PropertyStr *p = calloc(1, sizeof(struct PropertyStr));
        atom[i].property = p;
        atom[i].dynamic = calloc(1, sizeof(struct DynamicStr));

        fread(p->name, sizeof(char), 8, SysInfoFile);
        fread(p->extraProperty, sizeof(char), 10 * 16, SysInfoFile);
        fread(p->resName, sizeof(char), 8, SysInfoFile);
        fread(p->typeName, sizeof(char), 3 * 16, SysInfoFile);
        fread(&p->num, sizeof(int), 1, SysInfoFile);
        fread(p->typeofAtom, sizeof(int), 3, SysInfoFile);
        fread(&p->typeofRes, sizeof(int), 1, SysInfoFile);
        fread(&p->charge, sizeof(double), 1, SysInfoFile);
        fread(&p->mass, sizeof(double), 1, SysInfoFile);
        fread(p->color, sizeof(double), 3, SysInfoFile);
        fread(&p->sequence, sizeof(struct AtomSeqStr), 1, SysInfoFile);

        int nConstr = 0;
        fread(&nConstr, sizeof(int), 1, SysInfoFile);
        struct ConstraintStr **tail = &p->constr;
        for (int k = 0; k < nConstr; k++) {
            struct ConstraintStr *c = calloc(1, sizeof(struct ConstraintStr));
            ReadConstraint(c, SysInfoFile);
            *tail = c;
            tail = &c->next;
        }

        fread(&atom[i].dynamic->HB, sizeof(struct HBStr), 1, SysInfoFile);
    }

    fread(&HBPotential.typeNum, sizeof(int), 1, SysInfoFile);
    fread(HBPotential.scale, sizeof(double), HB_CLASSNUM, SysInfoFile);
    fread(HBPotential.barrier, sizeof(double), HB_CLASSNUM, SysInfoFile);
    HBPotential.classOf = (int *)calloc(HBPotential.typeNum + 1, sizeof(int));
    fread(HBPotential.classOf, sizeof(int), HBPotential.typeNum + 1, SysInfoFile);

    _CALLOC2(struct ConstraintStr, potentialPairCollision, nType, nType);
    for (int i = 0; i < nType; i++)
        for (int n = 0; n < nType; n++)
            ReadConstraint(&potentialPairCollision[i][n], SysInfoFile);

    potentialPairHB = (struct ConstraintStr ***)calloc(HBPotential.typeNum + 1, sizeof(struct ConstraintStr **));
    for (int i = 0; i <= HBPotential.typeNum; i++) {
        _CALLOC2(struct ConstraintStr, potentialPairHB[i], nType, nType);
    }
    for (int type = 0; type <= HBPotential.typeNum && HBPotential.typeNum > 0; type++)
        for (int i = 0; i < nType; i++)
            for (int n = 0; n < nType; n++)
                ReadConstraint(&potentialPairHB[type][i][n], SysInfoFile);
}


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
    
    // files from sDMD since the counted format start with a marker and a version;
    // older files start directly with the atom number and store raw structs (ReadSysInfoLegacy)
    int version = 0;
    fread(&atomnum,  sizeof (int), 1, SysInfoFile);
    if (atomnum == SYSINFO_MAGIC) {
        fread(&version, sizeof (int), 1, SysInfoFile);
        if (version != SYSINFO_VERSION) {
            printf("!!ERROR!!: unknown SysInfo.dat format version %i! %s:%i\n", version, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        fread(&atomnum, sizeof (int), 1, SysInfoFile);
    }
    fread(&numofMol, sizeof (int), 1, SysInfoFile);
    fread(&numofRes, sizeof (int), 1, SysInfoFile);

    _CALLOC2(int, connectionMap, (atomnum + 1), (atomnum + 1));
    
    residue = (struct ResidueStr *)calloc(numofRes + 1, sizeof(struct ResidueStr));
    fread(residue, sizeof (struct ResidueStr), (numofRes + 1), SysInfoFile);
    
    molecule = (struct MolStr *)calloc(numofMol + 1, sizeof(struct MolStr));
    fread(molecule, sizeof (struct MolStr), (numofMol + 1), SysInfoFile);
    
    //read force field info
    int num[3];
    fread(&num[0], sizeof(int), 1, SysInfoFile); num[0] ++;
    fread(&num[1], sizeof(int), 1, SysInfoFile); num[1] ++;
    FFInfo.group = (struct FFBeadStr *)calloc(num[0], sizeof(struct FFBeadStr));
    FFInfo.bead  = (struct FFBeadStr *)calloc(num[1], sizeof(struct FFBeadStr));
    FFInfo.LJ    = (struct FFLJStr   *)calloc(num[1], sizeof(struct FFLJStr));
    _CALLOC2(struct FFBondStr,     FFInfo.bond,     num[1], num[1]);
    _CALLOC2(struct FFAngleStr,    FFInfo.angle,    num[1], num[1]);
    _CALLOC2(struct FFDihedralStr, FFInfo.dihedral, num[1], num[1]);
    fread(FFInfo.group,    sizeof(struct FFBeadStr),     num[0],          SysInfoFile);
    fread(FFInfo.bead,     sizeof(struct FFBeadStr),     num[1],          SysInfoFile);
    fread(&num[2], sizeof(int), 1, SysInfoFile);
    if (num[2]) {
        fread(FFInfo.bond[0],     sizeof(struct FFBondStr),     num[1] * num[1], SysInfoFile);
    }
    fread(&num[2], sizeof(int), 1, SysInfoFile);
    if (num[2]) {
        fread(FFInfo.angle[0],    sizeof(struct FFAngleStr),    num[1] * num[1], SysInfoFile);
    }
    fread(&num[2], sizeof(int), 1, SysInfoFile);
    if (num[2]) {
        fread(FFInfo.dihedral[0], sizeof(struct FFDihedralStr), num[1] * num[1], SysInfoFile);
    }
    fread(&num[2], sizeof(int), 1, SysInfoFile);
    if (num[2]) {
        fread(FFInfo.LJ,          sizeof(struct FFLJStr),       num[1],          SysInfoFile);
    }

    for (int i = 0; i < num[1] * num[1]; i++) FFInfo.dihedral[0][i].list = NULL;   // stale pointer in the file

    if (version >= 2) {
        ReadSysInfoCounted(SysInfoFile, num[1]);
    } else {
        ReadSysInfoLegacy(SysInfoFile, num[1]);
    }

    CreateGroup();

    fclose(SysInfoFile);
        
    printf("Done!\n");
    fflush(stdout);

    return;
}


int ReadLog(int id) {
    int bufferSize = 1024 * 50;
    char *buffer, sign[4];
    FILE *inputFile;
    
    buffer = (char *)calloc(bufferSize, sizeof(char));
    for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
        sprintf(buffer, "%s%s", path, fileList[id].list[sectNum + 1]);
        inputFile = fopen(buffer, "r");
        if (inputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", fileList[id].list[sectNum + 1], path, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        struct SectionStr *sect = &sectInfo[sectNum];
        
        while (fgets(buffer, bufferSize, inputFile) != NULL) {
            if (buffer[0] == '#') {
                fgets(buffer, bufferSize, inputFile);   //System
                fgets(buffer, bufferSize, inputFile);   //Protein Num
                fgets(buffer, bufferSize, inputFile);   //Protein Seq
                
                for (int i = 0; i < numofMol; i ++) {
                    fgets(buffer, bufferSize, inputFile);
                }
                
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->oldTime);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->targetTime);

                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->temp);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->outputRate);
                
                if (sectNum == 0) {
                    fscanf(inputFile, "%s%s%lf\n", buffer, sign, &cutoffr); //cut off radius cannot change between sections
                } else {
                    fgets(buffer, bufferSize, inputFile);
                }
                
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->methodType);
                
                fgets(buffer, bufferSize, inputFile);   //ThermostatType
                fgets(buffer, bufferSize, inputFile);   //ThermostatParameter
                
                fgets(buffer, bufferSize, inputFile);   //considering the compatibility
                if (strncmp(buffer, "DMDMethod", 9) == 0) {
                    fgets(buffer, bufferSize, inputFile);
                    fgets(buffer, bufferSize, inputFile);
                }
                
                if (strncmp(buffer, "CGModel", 7) == 0) { //considering the compatibility
                    char tmp[1024] = { 0 };
                    sscanf(buffer, "%s%s%s%s\n", tmp, sign, CG.type[0], CG.type[1]);
                    if (strncmp(CG.type[0], "no", 2) != 0) CG.mark = 1;
                    fgets(buffer, bufferSize, inputFile);
                }

                char tmp[1024] = { 0 };
                sscanf(buffer, "%s%s%s\n",  tmp, sign, sect->wallObj.wallExist);
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->wallObj.wallType);
                
                if (strcmp(sect->wallObj.wallExist, "smooth") == 0) {
                    ReadWall(id, sectNum);
                }
                
                fscanf(inputFile, "%s%s%s\n",  buffer, sign, sect->wallDyn.mark);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->wallDyn.size);
                fscanf(inputFile, "%s%s%lf\n", buffer, sign, &sect->wallDyn.rate);
                
                if (sectNum == 0) { //original box size cannot change between sections
                    if (strcmp(sect->wallDyn.mark, "no") != 0) {
                        fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                               &boxOrigDim[1],
                               &boxOrigDim[2],
                               &boxOrigDim[3]);
                        fgets(buffer, bufferSize, inputFile);
                    } else {
                        fgets(buffer, bufferSize, inputFile);
                        fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                               &boxOrigDim[1],
                               &boxOrigDim[2],
                               &boxOrigDim[3]);
                    }
                } else {
                    fgets(buffer, bufferSize, inputFile);
                    fgets(buffer, bufferSize, inputFile);
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
                    
                    _CALLOC2(double, sect->flow.charge.position, sect->flow.charge.num, 4);
                    _CALLOC2(double, sect->flow.charge.velocity, sect->flow.charge.num, 4);
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
                    fgets(buffer, bufferSize, inputFile);   //the last \n
                } else {
                    fgets(buffer, bufferSize, inputFile);
                }
                
                fscanf(inputFile, "%s%s%i\n", buffer, sign, &sect->obstObj.num);
                if (sect->obstObj.num) {
                    sect->obstObj.mark = 1;
                    _CALLOC2(double, sect->obstObj.position, sect->obstObj.num, 4);
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
                    }
                    fgets(buffer, bufferSize, inputFile);   //the last \n
                    
                    ReadObst(id, sectNum);
                } else {
                    fgets(buffer, bufferSize, inputFile);
                }
                
                fscanf(inputFile, "%s%s%i\n", buffer, sign, &sect->tunlObj.mark);
                if (sect->tunlObj.mark) {
                    fscanf(inputFile, "%s%s%lf%lf%lf\n", buffer, sign,
                           &sect->tunlObj.startPosition,
                           &sect->tunlObj.endPosition,
                           &sect->tunlObj.diameter);
                    ReadWall(id, sectNum);
                } else {
                    fgets(buffer, bufferSize, inputFile);
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
    
    for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
        if (sectNum < fileList[id].count - 1) {
            sectInfo[sectNum].frameCount = sectInfo[sectNum + 1].oldTime - sectInfo[sectNum].oldTime;
        } else {
            sectInfo[sectNum].frameCount = sectInfo[sectNum].targetTime - sectInfo[sectNum].oldTime;
        }
    }
    
    free(buffer);
    return 1;
}


void ReadWall(int id, int sectNum) {
    char directory[1024];
    FILE *WallInfoFile;
    
    sprintf(directory, "%s%s", path, files[inWall][id].name);
    WallInfoFile = fopen(directory, "rb");
    if (WallInfoFile == NULL) {
        return;
    }
    
    struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
    fread(property, sizeof(struct PropertyStr), 1, WallInfoFile);
    sectInfo[sectNum].wallObj.wall.property = property;
    
    struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));
    fread(dynamic, sizeof(struct DynamicStr), 1, WallInfoFile);
    sectInfo[sectNum].wallObj.wall.dynamic = dynamic;
    
    fclose(WallInfoFile);
    return;
}


void ReadObst(int id, int sectNum) {
    char directory[1024];
    FILE *ObstInfoFile;
    
    sprintf(directory, "%s%s", path, files[inObst][id].name);
    ObstInfoFile = fopen(directory, "rb");
    if (ObstInfoFile == NULL) {
        return;
    }
    
    for (int i = 0; i < sectInfo[sectNum].obstObj.num; i ++) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        fread(property, sizeof(struct PropertyStr), 1, ObstInfoFile);
        sectInfo[sectNum].obstObj.obst[i].property = property;
        
        struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));
        fread(dynamic, sizeof(struct DynamicStr), 1, ObstInfoFile);
        sectInfo[sectNum].obstObj.obst[i].dynamic = dynamic;
    }
    
    fclose(ObstInfoFile);
    return;
}


void ReadCGPot(int id) {
    char directory[1024];
    FILE *CGInfoFile;
    
    printf("Collecting CG potential information...");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[inCG][id].name);
    CGInfoFile = fopen(directory, "rb");
    if (CGInfoFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inCG][id].name, path, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    _CALLOC2(struct ConstraintStr, potentialPairCG, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
    for (int i = 0; i <= FFInfo.bead[0].val[0]; i ++) {
        for (int n = 0; n <= FFInfo.bead[0].val[0]; n ++) {
            fread(&potentialPairCG[i][n], sizeof(struct ConstraintStr), 1, CGInfoFile);
            struct StepPotenStr *thisStep = NULL;
            struct StepPotenStr *nextStep = potentialPairCG[i][n].step;
            
            int flag2 = 0;
            while (nextStep) {
                struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
                fread(newStep, sizeof(struct StepPotenStr), 1, CGInfoFile);
                
                if (flag2 == 0) {
                    potentialPairCG[i][n].step = newStep;
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
    
	printf("Done!\n");
	fflush(stdout);

    fclose(CGInfoFile);
    return;
}

void FreeVariables(void) {
    free(residue);    residue = NULL;
    free(molecule);      molecule   = NULL;
    
    for (int i = 0; i <= atomnum; i ++) {
        free(atom[i].dynamic);
        free(atom[i].property);
    }
    free(atom);     atom     = NULL;
    free(celllist); celllist = NULL;

    for (int i = 0; i < NumofFileType; i++) {
        free(files[i]);
    }
    free(files);        files        = NULL;
    free(fileList);     fileList     = NULL;
    free(analysisList); analysisList = NULL;
    
    free(FFInfo.group);
    free(FFInfo.bead);
    
    if (FFInfo.LJ)       {  free(FFInfo.LJ);}
    if (FFInfo.bond)     {_FREE2(FFInfo.bond);}
    if (FFInfo.angle)    {_FREE2(FFInfo.angle);}
    if (FFInfo.dihedral) {_FREE2(FFInfo.dihedral);}
    if (potentialPairCG) {_FREE2(potentialPairCG);}
    
    _FREE2(connectionMap);
    _FREE2(potentialPairCollision);
    for (int i = 0; i <= HBPotential.typeNum; i ++) {
        _FREE2(potentialPairHB[i]);
    }
    free(potentialPairHB);

    free(group.list_Mol);
    free(group.list_Atom);
    free(group.molList);
    free(group.atomList);
    
    return;
}
