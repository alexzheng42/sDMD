//
//  Energy.c
//  Analysis
//
//  Created by Size Zheng on 11/2/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"

int HBModel(struct AtomStr *atom1, struct AtomStr *atom2);
int AtomTypeChange(int originalType, int direct);
double CalKeEnergy(void);
double CalPoEnergy(void);
double CalWlEnergy(struct SectionStr *sect);
double CalTemp(double TKE);
double FindPair (struct AtomStr *atom1, struct AtomStr *atom2, char *type, void *obstInfo);
struct ConstraintStr *RightPair(int type1, int type2, int flag);


void EnergyInfo(int id) {
    double TKE = 0, TPE = 0, TWE = 0, Temp;
    long step = 0;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE *TrjInputFile, *CntInputFile;
    FILE *KEOutputFile, *PEOutputFile, *WEOutputFile, *TolEOutputFile, *TempOutputFile;
    
    printf("Analyzing energy information... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outTemp][id].name);
    TempOutputFile = fopen(directory, "w");
    fprintf(TempOutputFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(TempOutputFile, "%s\n", "@    yaxis  label \"T (K)\"");
    
    sprintf(directory, "%s%s", path, files[outKE][id].name);
    KEOutputFile = fopen(directory, "w");
    fprintf(KEOutputFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(KEOutputFile, "%s\n", "@    yaxis  label \"E (kcal/mol)\"");
    
    sprintf(directory, "%s%s", path, files[outPE][id].name);
    PEOutputFile = fopen(directory, "w");
    fprintf(PEOutputFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(PEOutputFile, "%s\n", "@    yaxis  label \"E (kcal/mol)\"");
    
    sprintf(directory, "%s%s", path, files[outWE][id].name);
    WEOutputFile = fopen(directory, "w");
    fprintf(WEOutputFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(WEOutputFile, "%s\n", "@    yaxis  label \"E (kcal/mol)\"");
    
    sprintf(directory, "%s%s", path, files[outTolE][id].name);
    TolEOutputFile = fopen(directory, "w");
    fprintf(TolEOutputFile, "%s\n", "@    xaxis  label \"t\"");
    fprintf(TolEOutputFile, "%s\n", "@    yaxis  label \"E (kcal/mol)\"");
    
    for (int sectNum = 0; sectNum < fileList.count; sectNum ++) {
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inTrj][id].name, fileList.list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        TrjInputFile = fopen(directory, "r");
        if (TrjInputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inTrj][id].name, path, __FILE__, __LINE__);
            printf("           You may need to execute -rPBC first to generate the required file!\n");
            exit(EXIT_FAILURE);
        }
        
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inCnt][id].name, fileList.list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        CntInputFile = fopen(directory, "r");
        if (CntInputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inCnt][id].name, path, __FILE__, __LINE__);
            printf("           You may need to specify the connection information file by using -cnt flag!\n");
            exit(EXIT_FAILURE);
        }
        
        for (step = 0; step < sectInfo[sectNum].frameCount; step ++) {
            if (ReadGro(TrjInputFile) || ReadConnectionMap(CntInputFile)) break;
            
            PrintProcess(step);
            
            sect = &sectInfo[sectNum];
            TKE = CalKeEnergy();
            TPE = CalPoEnergy();
            if (strcmp(sect->wallObj.wallExist, "no") ||
                sect->tunlObj.mark) {
                TWE = CalWlEnergy(sect);
            }
            Temp = CalTemp(TKE);
            
            fprintf(KEOutputFile,   "%8.2lf%10.4lf\n", step * sect->outputRate + sect->oldTime, TKE);
            fprintf(PEOutputFile,   "%8.2lf%10.4lf\n", step * sect->outputRate + sect->oldTime, TPE + TWE);
            fprintf(WEOutputFile,   "%8.2lf%10.4lf\n", step * sect->outputRate + sect->oldTime, TWE);
            fprintf(TolEOutputFile, "%8.2lf%10.4lf\n", step * sect->outputRate + sect->oldTime, TKE + TPE + TWE);
            fprintf(TempOutputFile, "%8.2lf%10.4lf\n", step * sect->outputRate + sect->oldTime, Temp);
        }
        fclose(TrjInputFile);
        fclose(CntInputFile);
    }
    
    fclose(KEOutputFile);
    fclose(PEOutputFile);
    fclose(WEOutputFile);
    fclose(TolEOutputFile);
    fclose(TempOutputFile);
    
    printf("\bDone!\n");
    fflush(stdout);
    return;
}

double CalKeEnergy(void) { //per frame
    double temp;
    double TKE = 0;
    double speed[4] = {0};
    
    for (int n = 1; n <= atomnum; n++) {
        transfer_vector(speed, atom[n].dynamic->velocity);
        temp = dotprod(speed, speed);
        TKE += 0.5 * atom[n].property->mass * temp;
    }
    
    return TKE;
}

double CalPoEnergy(void) { //per frame
    signed int x[27]={0,1,0,0,0,1,1,1,-1,0,0,0,-1,-1,-1,0,0,1,-1,1,-1,-1,-1,1,1,1,-1};
    signed int y[27]={0,0,1,0,1,0,1,1,0,-1,0,-1,0,-1,-1,-1,1,0,0,-1,1,-1,1,-1,1,-1,1};
    signed int z[27]={0,0,0,1,1,1,0,1,0,0,-1,-1,-1,0,-1,1,-1,-1,1,0,0,1,-1,-1,-1,1,1};
    int pair;
    int cell_neighbor[4];
    int cellindex;
    int *atom_CellAxis;
    double TPE = 0;
    double *boxsize;
    double positionshift[4] = {0};
    

    for (int i = 1; i <= atomnum; i ++) {
        atom[i].property->type = AtomTypeChange(atom[i].property->type, 0);
        if (CheckHBConnection(i)) {
            atom[i].property->type = AtomTypeChange(atom[i].property->type, 1);
        }
    }
    
    LinkList();
    pointToStruct(boxsize, boxOrigDim);
    
    for (int i = 1; i <= atomnum; i ++) {
        pointToStruct(atom_CellAxis, atom[i].dynamic->cellIndex);
        
        for (int n = 0; n <= 26; n ++) {
            cell_neighbor[3] = atom_CellAxis[3] + z[n];
            cell_neighbor[2] = atom_CellAxis[2] + y[n];
            cell_neighbor[1] = atom_CellAxis[1] + x[n]; //scan the neighborhood 27 subcells, include the target subcell itself
            
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
            while (pair != 0 && pair > i) {
                dotplus(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);
                TPE += FindPair(&atom[i], &atom[pair], "normal", NULL);
                
                pair = celllist[pair];
            }
        }
    }
    
    /*
    for (int i = 1; i <= atomnum; i ++) {
        for (int n = i + 1; n <= atomnum; n ++) {
            TPE += FindPair(&atom[i], &atom[n], "normal");
        }
    }*/
    
    return TPE;
}

double CalWlEnergy(struct SectionStr *sect) { //per frame
    double TWE = 0;
    
    for (int i = 1; i <= atomnum; i ++) {
        atom[i].property->type = AtomTypeChange(atom[i].property->type, 0);
        if (CheckHBConnection(i)) {
            atom[i].property->type = AtomTypeChange(atom[i].property->type, 1);
        }
    }
    
    if (strcmp(sect->wallObj.wallExist, "smooth") == 0) {
        for (int i = 1; i <= atomnum; i ++) {
            TWE += FindPair(&atom[i], &sect->wallObj.wall, "wall", &sect->wallObj);
        }
    }
    
    if (sect->tunlObj.mark) {
        for (int i = 1; i <= atomnum; i ++) {
            if (atom[i].dynamic->coordinate[1] >= sect->tunlObj.startPosition &&
                atom[i].dynamic->coordinate[1] <= sect->tunlObj.endPosition) {
                TWE += FindPair(&atom[i], &sect->tunlObj.tunnel, "tunnel", &sect->tunlObj);
            }
        }
    }
    
    return TWE;
}

double CalTemp(double TKE) { //per frame
    return 2 * TKE / (3 * atomnum) * 500; //remove "500" will be the reduced T
}

double FindPair (struct AtomStr *atom1, struct AtomStr *atom2, char *type, void *obstInfo)
{
    int type1 = 0, type2 = 0;
    int atom_i = atom1->property->num;
    int atom_j = atom2->property->num;
    double tolerance, floatzero = 10E-8;
    double position_i[4], position_j[4];
    double speed_i[4], speed_j[4];
    double r_ij[4], v_ij[4], direction, distance2;
    double PE = 0, accumPotential = 0;
    struct ConstraintStr *thisConstr = NULL;
    
    transfer_vector(position_i, atom1->dynamic->coordinate);
    transfer_vector(position_j, atom2->dynamic->coordinate);
    
    if (strcmp(type, "wall") == 0 && strcmp(((struct WallStr *)obstInfo)->wallExist, "smooth") == 0) {
        if (strncmp(((struct WallStr *)obstInfo)->wallType, "parallel", 1) == 0) { //parallel plates of wall
            position_i[1] = position_i[3] = 0;
        } else if (strncmp(((struct WallStr *)obstInfo)->wallType, "cylinder", 1) == 0) {
            position_i[1] = 0;
        }
    } else if (strcmp(type, "tunnel") == 0) {
        position_i[1] = 0;
    }
    
    dotminus(position_i, position_j, r_ij);
    distance2 = dotprod(r_ij, r_ij);
    if (strcmp(type, "wall") == 0 && strcmp(((struct WallStr *)obstInfo)->wallExist, "smooth") == 0) {
        distance2 = boxCurtDim[2] * boxCurtDim[2] * 0.25 + distance2 - boxCurtDim[2] * sqrt(distance2);
    } else if (strcmp(type, "tunnel") == 0) {
        double diameter = ((struct TunnelStr *)obstInfo)->diameter;
        distance2 = diameter * diameter * 0.25 + distance2 - diameter * sqrt(distance2);
    } else if (strcmp(type, "normal") == 0) {
    } else {
        printf("!!ERROR!!: FindPair type is invalid! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    if (distance2 <= cutoffr * cutoffr && !(connectionMap[atom_i][atom_j] & bondConnect)) {
        transfer_vector(speed_i, atom1->dynamic->velocity);
        transfer_vector(speed_j, atom2->dynamic->velocity);
        
        if (strcmp(type, "wall") == 0 && strcmp(((struct WallStr *)obstInfo)->wallExist, "smooth") == 0) {
            if (strncmp(((struct WallStr *)obstInfo)->wallType, "parallel", 1) == 0) { //parallel plates of wall
                speed_i[1] = speed_i[3] = 0;
            } else if (strncmp(((struct WallStr *)obstInfo)->wallType, "cylinder", 1) == 0) {
                speed_i[1] = 0;
            }
        } else if (strcmp(type, "tunnel") == 0) {
            speed_i[1] = 0;
        }
        
        dotminus(speed_i, speed_j, v_ij);
        direction = dotprod(r_ij, v_ij);
        
        if ((strcmp(type, "wall") == 0 && strcmp(((struct WallStr *)obstInfo)->wallExist, "smooth") == 0) ||
            strcmp(type, "tunnel") == 0)
            direction *= -1;
        
        if (direction < 0) {
            direction = 1; // -> <-
        } else {
            direction = -1; // <- ->
        }
        tolerance = direction * floatzero;
        
        if (connectionMap[atom_i][atom_j] & HBConnect) {
            type1 = HBModel(atom1, atom2);
            if (type1 == 1) {
                PE -= HBPotential.BB;
            } else {
                PE -= HBPotential.BS;
            }
            
            thisConstr = &potentialPairHB[type1][0][0];
        } else if (connectionMap[atom_i][atom_j] & constraintConnect) {
            thisConstr = atom1->property->constr;
            while (thisConstr->connection != atom_j)
                thisConstr = thisConstr->next;
        } else if (connectionMap[atom_i][atom_j] & neighborConnect) {
            int typeNum = 0;
            int HBPartner = 0;
            
            if ((HBPartner = CheckHBConnection(atom_i)) > 0 && atom[HBPartner].dynamic->HB.neighbor == atom_j) {
                typeNum = HBModel(atom1, &atom[HBPartner]);
                type1 = AtomTypeChange(atom2->property->type, 0);
                type2 = AtomTypeChange(atom1->property->type, 0);
            } else {
                HBPartner = CheckHBConnection(atom_j);
                typeNum = HBModel(atom2, &atom[HBPartner]);
                type1 = AtomTypeChange(atom1->property->type, 0);
                type2 = AtomTypeChange(atom2->property->type, 0);
            }
            
            thisConstr = &potentialPairHB[typeNum][type1][type2];
        } else {
            type1 = atom1->property->type;
            type2 = atom2->property->type;
            thisConstr = RightPair(type1, type2, 0);
        }
        
        struct StepPotenStr *thisStep = thisConstr->step;
        while (thisStep != NULL && thisStep->d != 0) {
            if (distance2 < thisStep->d + tolerance) {
                accumPotential = thisStep->accumulated;
                return PE + accumPotential;
            }
            thisStep = thisStep->next;
        }
    }
    
    return 0;
}

int HBModel(struct AtomStr *atom1, struct AtomStr *atom2) {
    struct AtomStr *switcher;
    
    if (atom1->dynamic->HB.role == 'A') { //atom1 -> donator, atom2 -> acceptor
        switcher = atom1;
        atom1 = atom2;
        atom2 = switcher;
    }
    
    if (AtomTypeChange(atom1->property->type, 0) == 2 /*HB*/ &&
        AtomTypeChange(atom2->property->type, 0) == 26 /*OZB*/) {
        return 1;
    } else if (AtomTypeChange(atom1->property->type, 0) == 2) {
        if (strcmp(atom2->property->nameOfAA, "ASP") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLU") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLN") == 0 ||
            strcmp(atom2->property->nameOfAA, "ASN") == 0) {
            return 2;
        } else if (strcmp(atom2->property->nameOfAA, "HIS") == 0 ||
                   strcmp(atom2->property->nameOfAA, "TRP") == 0) {
            return 4;
        } else if (AtomTypeChange(atom2->property->type, 0) == 31) {
            return 10;
        } else {
            return 3;
        }
    } else if (AtomTypeChange(atom2->property->type, 0) == 26) {
        return 5;
    } else if (AtomTypeChange(atom1->property->type, 0) == 31 &&
               AtomTypeChange(atom2->property->type, 0) == 31) {
        return 9;
    } else {
        if (strcmp(atom2->property->nameOfAA, "ASP") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLU") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLN") == 0 ||
            strcmp(atom2->property->nameOfAA, "ASN") == 0) {
            return 6;
        } else if (strcmp(atom2->property->nameOfAA, "HIS") == 0 ||
                   strcmp(atom2->property->nameOfAA, "TRP") == 0) {
            return 8;
        } else {
            return 7;
        }
    }
    
    return -1;
}


int NeighborModel(char *atom1, char *atom2) {
    int num = 0;
    
    if (strncmp(atom1, "C", 1) == 0) {
        num = 3;
    } else if (strcmp(atom1, "OW") == 0) {
        num = 5;
    } else {
        if (strncmp(atom1, "N", 1) == 0) {
            if (strncmp(atom2, "O", 1) == 0) {
                num = 2;
            } else {
                num = 6;
            }
        } else {
            if (strncmp(atom2, "O", 1) == 0) {
                num = 5;
            } else {
                num = 2;
            }
        }
    }
    
    return num;
}

int AtomTypeChange(int originalType, int direct) {
    if (direct == 1) {
        switch (originalType) {
            case 1:             //H
                return 3;       //H_HB
                break;
            case 2:             //HB
                return 4;       //HB_HB
                break;
            case 26:
                return 27;
                break;
                
            default:
                break;
        }
    } else if (direct == 0) {
        switch (originalType) {
            case 3:             //H_HB
                return 1;       //H
                break;
            case 4:             //HB_HB
                return 2;       //HB
                break;
            case 27:            //OZB_HB
                return 26;      //OZB
                break;
                
            default:
                break;
        }
    } else {
        printf("!!ERROR!!: atom type change direction is invalid!\n");
        exit(EXIT_FAILURE);
    }
    
    return originalType;
}

struct ConstraintStr *RightPair(int type1, int type2, int flag) {
    
    if (potentialPairCollision[type1][type2].dmin > 0) {
        return &potentialPairCollision[type1][type2];
    } else if (flag) {
        return NULL;
    }
    
    struct ConstraintStr *targetPair = NULL;
    if ((targetPair = RightPair(AtomTypeChange(type1, 0), type2, 1)) != NULL) {}
    else if ((targetPair = RightPair(type1, AtomTypeChange(type2, 0), 1)) != NULL) {}
    else if ((targetPair = RightPair(AtomTypeChange(type1, 0), AtomTypeChange(type2, 0), 1)) != NULL) {}
    
    return targetPair;
}

