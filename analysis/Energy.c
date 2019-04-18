//
//  Energy.c
//  Analysis
//
//  Created by Size Zheng on 11/2/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"

struct EnergyReadStr energy;

static int HBModel(struct AtomStr *atom1, struct AtomStr *atom2);
static int AtomTypeChange(int originalType, int direct);
static double CalKeEnergy(void);
static double CalPoEnergy(void);
static double CalWlEnergy(struct SectionStr *sect);
static double CalTemp(double TKE);
static double FindPair (struct AtomStr *atom1, struct AtomStr *atom2, char *type, void *obstInfo);
static double HBBarrier(struct AtomStr *atom1, struct AtomStr *atom2);
static struct ConstraintStr *RightPair(int type1, int type2, int flag);


static int beginP, endP, offSet, thisAtomNum;


void EnergyInfo(int id) {
    double TKE = 0, TPE = 0, TWE = 0, Temp = 0;
    long step = 0, totalFrame;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE *TrjInputFile, *CntInputFile;
    FILE *EneOutputFile;
    
    beginP = 1;
    endP = numofprotein;
    thisAtomNum = atomnum;
    offSet = 0;
    
    if (nPP) {
        beginP = endP = nPP;
        thisAtomNum = protein[nPP].endAtomNum - protein[nPP].startAtomNum + 1;
        offSet = protein[nPP].startAtomNum - 1;
    }
    
    printf("Analyzing energy information... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outEne][id].name);
    EneOutputFile = fopen(directory, "w");
    fprintf(EneOutputFile, "# %13s %17s %17s %17s %17s %17s\n",
            "Step",
            "Temp (K)",
            "KinE (kCal/mol)",
            "PotE (kCal/mol)",
            "WalE (kCal/mol)",
            "TotE (kCal/mol)");
    
    for (int sectNum = 0; sectNum < fileList.count; sectNum ++) {
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inTrj][id].name, fileList.list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        TrjInputFile = fopen(directory, "r");
        if (TrjInputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", buffer, path, __FILE__, __LINE__);
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
        
        totalFrame = (sectInfo[sectNum].oldTime + sectInfo[sectNum].frameCount) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
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
            
            fprintf(EneOutputFile, "%15.2lf %17.4lf %17.4lf %17.4lf %17.4lf %17.4lf\n",
                    step * sect->outputRate,
                    Temp,
                    TKE,
                    TPE + TWE,
                    TWE,
                    TKE + TPE + TWE);
        }
        fclose(TrjInputFile);
        fclose(CntInputFile);
    }
    fclose(EneOutputFile);
    
    printf("\bDone!\n");
    fflush(stdout);
    return;
}

double CalKeEnergy(void) { //per frame
    double temp;
    double TKE = 0;
    double speed[4] = {0};
    
    for (int n = 1 + offSet; n <= thisAtomNum + offSet; n++) {
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
        atom[i].property->typeofAtom = AtomTypeChange(atom[i].property->typeofAtom, 0);
        if (CheckHBConnection(i)) {
            atom[i].property->typeofAtom = AtomTypeChange(atom[i].property->typeofAtom, 1);
        }
    }
    
    LinkList();
    pointToStruct(boxsize, boxOrigDim);
    
    for (int i = 1 + offSet; i <= thisAtomNum + offSet; i ++) {
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
            if (!nPP) {
                while (pair != 0) {
                    if (pair > i) {
                        dotplus(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);
                        TPE += FindPair(&atom[i], &atom[pair], "normal", NULL);
                        dotminus(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);
                    }
                    
                    pair = celllist[pair];
                }
            } else {
                while (pair != 0) {
                    if (pair != i /*&& atom[pair].property->sequence.proteinNum == atom[i].property->sequence.proteinNum*/) {
                        dotplus(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);
                        TPE += FindPair(&atom[i], &atom[pair], "normal", NULL) * 0.5;
                        dotminus(atom[pair].dynamic->coordinate, positionshift, atom[pair].dynamic->coordinate);
                    }
                    
                    pair = celllist[pair];
                }
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
        atom[i].property->typeofAtom = AtomTypeChange(atom[i].property->typeofAtom, 0);
        if (CheckHBConnection(i)) {
            atom[i].property->typeofAtom = AtomTypeChange(atom[i].property->typeofAtom, 1);
        }
    }
    
    if (strcmp(sect->wallObj.wallExist, "smooth") == 0) {
        for (int i = 1 + offSet; i <= thisAtomNum + offSet; i ++) {
            TWE += FindPair(&atom[i], &sect->wallObj.wall, "wall", &sect->wallObj);
        }
    }
    
    if (sect->tunlObj.mark) {
        for (int i = 1 + offSet; i <= thisAtomNum + offSet; i ++) {
            if (atom[i].dynamic->coordinate[1] >= sect->tunlObj.startPosition &&
                atom[i].dynamic->coordinate[1] <= sect->tunlObj.endPosition) {
                TWE += FindPair(&atom[i], &sect->tunlObj.tunnel, "tunnel", &sect->tunlObj);
            }
        }
    }
    
    return TWE;
}

double CalTemp(double TKE) { //per frame
    return 2 * TKE / (3 * thisAtomNum) / BOLTZMANN; //remove "BOLTZMANN" will be the reduced T
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
            PE -= HBBarrier(atom1, atom2);
            type1 = HBModel(atom1, atom2);
            
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
                type1 = AtomTypeChange(atom2->property->typeofAtom, 0);
                type2 = AtomTypeChange(atom1->property->typeofAtom, 0);
            } else {
                HBPartner = CheckHBConnection(atom_j);
                typeNum = HBModel(atom2, &atom[HBPartner]);
                type1 = AtomTypeChange(atom1->property->typeofAtom, 0);
                type2 = AtomTypeChange(atom2->property->typeofAtom, 0);
            }
            
            thisConstr = &potentialPairHB[typeNum][type1][type2];
        } else {
            type1 = atom1->property->typeofAtom;
            type2 = atom2->property->typeofAtom;
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
    
    if (atom1->dynamic->HB.role == 'A' ||
        atom1->dynamic->HB.role == 'S') { //atom1 -> donator, atom2 -> acceptor
        switcher = atom1;
        atom1 = atom2;
        atom2 = switcher;
    }
    
    int atom1_OrigType = AtomTypeChange(atom1->property->typeofAtom, 0);
    int atom2_OrigType = AtomTypeChange(atom2->property->typeofAtom, 0);
    
    if (atom1_OrigType == 2 /*HB*/ &&
        atom2_OrigType == 26 /*OZB*/) {
        return 1;
    } else if (atom1_OrigType == 2) { //backbone H, HB
        switch (atom2->property->typeofAA) {
            case 4:  //ASP
            case 7:  //GLU
                return 2;
                break;
                
            case 16: //SER
            case 17: //THR
            case 19: //TYR
                return 3;
                
            case 9:  //HIS
                return 4;
                
            case 5: //CYS
                return 10;
                
            case 3: //ASN
            case 6: //GLN
                //no HB between ASN/GLN and backbone H
                return -1;
                
            default:
                break;
        }
    } else if (atom2_OrigType == 26) { //backbone O, OZB
        return 5;
    } else if (atom1_OrigType == 31 &&
               atom2_OrigType == 31) { //SG
        return 9;
    } else if (atom1_OrigType == 1) { //sidechain H
        switch (atom2->property->typeofAA) {
            case 4:  //ASP
            case 7:  //GLU
                return 6;
                break;
                
            case 16: //SER
            case 17: //THR
            case 19: //TYR
                return 7;
                
            case 9:  //HIS
                return 8;
                
            case 3:  //ASN
            case 6:  //GLN
                return 11;
                
            default:
                break;
        }
    }
    
    printf("!!ERROR!!: cannot find the matched HB pair for atom %2i(%2i%2s) and atom %2i(%2i%2s)!\n",
           atom1->property->num, atom1->property->sequence.aminoacidNum, atom1->property->name,
           atom2->property->num, atom2->property->sequence.aminoacidNum, atom2->property->name);
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    
    return -1;
}


double HBBarrier(struct AtomStr *atom1, struct AtomStr *atom2) {
    struct AtomStr *switcher;
    
    if (atom1->dynamic->HB.role == 'A' ||
        atom1->dynamic->HB.role == 'S') { //atom1 -> donator, atom2 -> acceptor
        switcher = atom1;
        atom1 = atom2;
        atom2 = switcher;
    }
    
    int atom1_OrigType = AtomTypeChange(atom1->property->typeofAtom, 0);
    int atom2_OrigType = AtomTypeChange(atom2->property->typeofAtom, 0);
    
    if (atom1_OrigType == 2 /*HB*/ &&
        atom2_OrigType == 26 /*OZB*/) {
        return HBPotential.BB_v;
    } else if (atom1_OrigType != 2 &&
               atom2_OrigType != 26) {
        return HBPotential.SS_v;
    } else {
        return HBPotential.BS_v;
    }
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


int ReadEnergyFile(FILE *EnergyInputFile, struct EnergyReadStr *thisE) {
    char buffer[1024];
    
repeat:
    if (fgets(buffer, sizeof(buffer), EnergyInputFile) == NULL)
        return -1;
    
    if (buffer[0] == '@' || buffer[0] == '#') {
        goto repeat;
    }
    
    sscanf(buffer, "%lf%lf%lf%lf%lf%lf",
           &thisE->step,
           &thisE->T,
           &thisE->KE,
           &thisE->PE,
           &thisE->WE,
           &thisE->TE);
    
    return 0;
}
