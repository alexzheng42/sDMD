//
//  ToolFunctions.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"


void TimeForward(double time, struct ThreadStr *thisThread) {
    for (int i = 1; i <= atomnum; ++i) {
        thisThread->raw[i]->dynamic->event.time -= time;
    }
}


void UpdateData(double time, char *type, struct ThreadStr *thisThread) {
    if (strncmp(type, "atom", 1) == 0) {
        for (int i = 1; i <= atomnum; i++) {
            thisThread->raw[i]->dynamic->coordinate[1] += thisThread->raw[i]->dynamic->velocity[1] * time;
            thisThread->raw[i]->dynamic->coordinate[2] += thisThread->raw[i]->dynamic->velocity[2] * time;
            thisThread->raw[i]->dynamic->coordinate[3] += thisThread->raw[i]->dynamic->velocity[3] * time;
        }
    } else if (strncmp(type, "partial", 2) == 0) {
        list *atomList = &thisThread->atomList;
        listElem *elem = listFirst(atomList);
        for (int n = 1; n <= atomList->num_members; n ++) {
            ((struct AtomStr *)elem->obj)->dynamic->coordinate[1] += ((struct AtomStr *)elem->obj)->dynamic->velocity[1] * time;
            ((struct AtomStr *)elem->obj)->dynamic->coordinate[2] += ((struct AtomStr *)elem->obj)->dynamic->velocity[2] * time;
            ((struct AtomStr *)elem->obj)->dynamic->coordinate[3] += ((struct AtomStr *)elem->obj)->dynamic->velocity[3] * time;
            
            elem = listNext(atomList, elem);
        }
    } else {
        printf("!ERROR!: UpdateData function type has something wrong!\n");
    }
}


int FindPair (struct AtomStr *atom1, struct AtomStr *atom2, char *interactionType, double direction, double distance2, double *lowerLimit, double *upperLimit, double *lowerPotential, double *upperPotential, double* accumPotential, struct AtomStr *HBPartner, int typeChange)
{
    int type1 = 0, type2 = 0;
    int coreorShell = -1; //0 -> core, 1 -> shell, -1 -> nothing
    double tolerance;
    struct ConstraintStr *thisConstr = NULL;
    
    if (direction < 0) {
        direction = 1; // -> <-
    } else {
        direction = -1; // <- ->
    }
    tolerance = direction * ZERO;
    
    if (strncmp(interactionType, "collision", 3) == 0) {
        
        if (typeChange == 0) { //not change
            type1 = atom1->property->type;
        } else if (typeChange == 1) { //change to _HB
            type1 = AtomTypeChange(atom1->property->type, 1);
        } else if (typeChange == -1) { //change back to normal
            type1 = AtomTypeChange(atom1->property->type, 0);
        } else {
            printf("!!ERROR!!: typeChange value is invalid. %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        type2 = atom2->property->type;
        thisConstr = RightPair(type1, type2, 0);
        
    } else if (strncmp(interactionType, "constraint", 3) == 0) {
        
        thisConstr = atom1->property->constr;
        while (thisConstr->connection != atom2->property->num) {
            thisConstr = thisConstr->next;
            
            if (unlikely(thisConstr == NULL)) {
                printf("!!ERROR!!: cannot find the atom of constraint connection atom #%i for atom #%i! %s:%i\n", atom1->property->num, atom2->property->num, __FILE__, __LINE__);
            }
        }
        
    } else if (strncmp(interactionType, "HB", 1)==0) { //don't need to assign potential here
        
        type1 = HBModel(atom1, atom2);
        thisConstr = &potentialPairHB[type1][0][0];
        
    } else if (strncmp(interactionType, "neighbor", 1)==0) {
        
        int typeNum = HBModel(atom2, HBPartner);
        type1 = AtomTypeChange(atom1->property->type, 0);
        type2 = AtomTypeChange(atom2->property->type, 0);
        
        thisConstr = &potentialPairHB[typeNum][type1][type2];
        
    } else {
        printf("!ERROR!: InteractionType assignment does not match!\n");
        exit(EXIT_FAILURE);
    }
    
    /*
     since some of the constrains only have dmin and dmax
     its coreorShell assignment would be special
     need to consider the moving direction
     -> <-: 0
     <- ->: 1
     */
    int flag = (thisConstr->step == NULL);
    struct StepPotenStr *thisStep = thisConstr->step;
    
    if (unlikely(thisConstr->dmin == 0)) {
        printf("!!ERROR!!: the potential pair is invalid! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (distance2 < thisConstr->dmin - ZERO) {
        
        *lowerLimit = INFINIT;
        *upperLimit = INFINIT;
        *accumPotential = INFINIT;
        *lowerPotential = 0;
        *upperPotential = 0;
        
        if (unlikely(!(strcmp(interactionType, "neighbor") == 0 && atom2->dynamic->HB.bondConnection == 0) &&
                     !(strcmp(interactionType, "HB") == 0 && atom1->dynamic->HB.bondConnection == 0) &&
                     !(strcmp(wallDyn.mark, "no")) &&
                     !tunlObj.mark &&
                     !(codeNum > 1))) {
            printf("!WARNING!: the distance between atom %2i(%2i%2s) and %2i(%2i%2s) is too small!\n",
                   atom1->property->num, atom1->property->sequence.aminoacidNum, atom1->property->name,
                   atom2->property->num, atom2->property->sequence.aminoacidNum, atom2->property->name);
            printf("           pair type: %s. the current distance is %.4lf, while the min is %.4lf.\n", interactionType, sqrt(distance2), sqrt(thisConstr->dmin));
            printf("           %s:%i\n", __FILE__, __LINE__);
            warningsum ++;
        }
        return coreorShell;
    }
    
    *lowerLimit = thisConstr->dmin;
    *lowerPotential = 0;
    
    while (thisStep != NULL && thisStep->d != 0) {
        if (distance2 < thisStep->d + tolerance) {
            *upperLimit = thisStep->d;
            *upperPotential = -1 * thisStep->e;
            *accumPotential = thisStep->accumulated;
            
            if (thisStep == thisConstr->step) {
                coreorShell = 0;
                return coreorShell;
            } else if (strcmp(interactionType, "HB") == 0 &&
                       thisStep->next == NULL && atom1->dynamic->HB.bondConnection > 0) {
                coreorShell = 1;
                return coreorShell;
            } else if (strcmp(interactionType, "neighbor") == 0 &&
                       thisStep->next == NULL && atom2->dynamic->HB.bondConnection > 0) {
                *upperPotential = 0;
                coreorShell = 1; //for HBNeighbor event, the paired atoms will NOT move out of the range, i.e. the last step->d = dmax
                return coreorShell;
            }
            
            return coreorShell;
        }
        *lowerLimit = thisStep->d;
        *lowerPotential = thisStep->e;
        thisStep = thisStep->next;
    }
    
    if (unlikely(strcmp(interactionType, "HB") == 0 && atom1->dynamic->HB.bondConnection > 0)) {
        printf("!WARNING!: HB length between atom %i and %i is too large!\n",
               atom1->property->sequence.atomNum,
               atom2->property->sequence.atomNum);
        printf("           the current is %.4lf, while the max is %.4lf\n", sqrt(distance2), sqrt(*upperLimit));
        printf("           %s:%i\n", __FILE__, __LINE__);
        warningsum ++;
    }
    
    if (thisConstr->dmax > 0) {
        if (distance2 < thisConstr->dmax + ZERO) {
            *upperLimit = thisConstr->dmax;
            *upperPotential = 0;
            *accumPotential = 0;
            
            coreorShell = 1;
            
            //only have dmin and dmax && -> <-
            //here direction = 1 => -> <-
            if (flag == 1 && direction == 1) {
                coreorShell = 0;
            }
            
            return coreorShell;
        } else {
            *lowerLimit=INFINIT;
            *upperLimit=INFINIT;
            *accumPotential = INFINIT;
            *lowerPotential=0;
            *upperPotential=0;
            
            if (unlikely(!(strcmp(interactionType, "neighbor") == 0 && atom2->dynamic->HB.bondConnection == 0) &&
                         !(strcmp(interactionType, "HB") == 0 && atom1->dynamic->HB.bondConnection == 0))) {
                printf("!WARNING!: the distance between atom %2i(%2i%2s) and %2i(%2i%2s) is too large!\n",
                       atom1->property->num, atom1->property->sequence.aminoacidNum, atom1->property->name,
                       atom2->property->num, atom2->property->sequence.aminoacidNum, atom2->property->name);
                printf("           pair type: %s. the current distance is %.4lf, while the min is %.4lf.\n", interactionType, sqrt(distance2), sqrt(thisConstr->dmax));
                printf("           %s:%i\n", __FILE__, __LINE__);
                warningsum ++;
            }
        }
    } else {
        
        //if still not return, then the distance is out of potential interaction
        *upperLimit = INFINIT;
        *upperPotential = INFINIT;
        *accumPotential = INFINIT;
        
        coreorShell = 1;
        return coreorShell;
    }
    
    return coreorShell;
}


/*
 after HB forming, some of the atoms' types will change
 after HB breaking, the types will change back
 direct: 0 -> back to normal
         1 -> change to _HB
 */
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


/*
 after HB forming, the type changed atoms may have no
 cooresponding potential pair for collision event,
 we need to use their original types to find the pair
 */
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


void PrintCellList(struct ThreadStr *thisThread) {
    int num;
    
    num = atomnum + cellnum[1] * cellnum[2] * cellnum[3];
    
    for (int i = 1; i <= num; ++i) {
        printf("%3i = %3i", i, celllist[i]);
        if (i > atomnum && celllist[i] != 0) {
            printf(", %2i:%2i:%2i\n",
                   thisThread->raw[celllist[i]]->dynamic->cellIndex[1],
                   thisThread->raw[celllist[i]]->dynamic->cellIndex[2],
                   thisThread->raw[celllist[i]]->dynamic->cellIndex[3]);
        } else {
            printf("\n");
        }
    }
}


void DisplayTime(char * timer) {
    time_t timetodisplay;
    struct tm* tm_info;
    
    time(&timetodisplay);
    tm_info = localtime(&timetodisplay);
    
    strftime(timer, 30, "%Y_%m_%d_%H_%M_%S", tm_info);
}


void SDEnergyMin(long int stepNum, struct ThreadStr *thisThread) {
    struct interStr {
        int connect[512];
        double coordinate[4];
        double hiLimit[512];
        double loLimit[512];
        double length[512];
    };
    
    long step = 0;
    long totalCN = INFINIT * INFINIT;
    int l, flag;
    double factor = 0.1;
    double diff;
    double energy;
    double direction[4];
    double length2;
    char directory[100];
    struct interStr *interaction;
    struct FileStr outputFile;
    
    interaction = (struct interStr *)calloc(atomnum + 1, sizeof(struct interStr));
    sprintf(directory, "%s/MiniEnergy%s.gro", datadir, REMDInfo.REMD_ExtraName);
    outputFile.file = fopen(directory, "w");
    
    for (int i = 1; i <= atomnum; ++i) {
        
        TRANSFER_VECTOR(interaction[i].coordinate, thisThread->raw[i]->dynamic->coordinate);
        
        //transfer bond information
        int m = 0;
        struct ConstraintStr *thisBond = thisThread->raw[i]->property->bond;
        while (thisBond != NULL) {
            m++;
            interaction[i].connect[m] = thisBond->connection;
            interaction[i].hiLimit[m] = sqrt(thisBond->dmax);
            interaction[i].loLimit[m] = sqrt(thisBond->dmin);
            interaction[i].length[m]  = 0.5 * (interaction[i].hiLimit[m] + interaction[i].loLimit[m]);
            thisBond = thisBond->next;
        }
        
        //transfer constrains information
        struct ConstraintStr *thisConstr = thisThread->raw[i]->property->constr;
        while (thisConstr != NULL) {
            m++;
            interaction[i].connect[m] = thisConstr->connection;
            interaction[i].hiLimit[m] = sqrt(thisConstr->dmax);
            interaction[i].loLimit[m] = sqrt(thisConstr->dmin);
            interaction[i].length[m]  = 0.5 * (interaction[i].hiLimit[m] + interaction[i].loLimit[m]);
            thisConstr = thisConstr->next;
        }
        
        //transfer interaction potential information
        for (int n = 1; n <= atomnum; n ++) {
            direction[1] = thisThread->raw[n]->dynamic->coordinate[1] - thisThread->raw[i]->dynamic->coordinate[1];
            direction[2] = thisThread->raw[n]->dynamic->coordinate[2] - thisThread->raw[i]->dynamic->coordinate[2];
            direction[3] = thisThread->raw[n]->dynamic->coordinate[3] - thisThread->raw[i]->dynamic->coordinate[3];
            
            length2 = DOT_PROD(direction, direction);
            if (length2 <= 144 && //radius of 12 A
                !(connectionMap[i][n] & BOND_CONNECT) &&
                !(connectionMap[i][n] & CONSTRAINT_CONNECT)) {
                FIND_DUPLICATED(n, interaction[i].connect, l, flag);
                if (i != n && flag == 0) {
                    m++;
                    interaction[i].connect[m] = n; //all the surrounding atoms
                    interaction[i].hiLimit[m] = INFINIT; //interaction potential only has low limit
                    
                    int type1 = thisThread->raw[i]->property->type;
                    int type2 = thisThread->raw[n]->property->type;
                    
                    if (type1 > type2) {
                        int i = type2;
                        type2 = type1;
                        type1 = i;
                    }
                    
                    interaction[i].loLimit[m] = sqrt(potentialPairCollision[type1][type2].dmin);
                    interaction[i].length[m]  = interaction[i].loLimit[m] + 0.001;
                }
            }
        }
        
        if (m > 510) {
            printf("!!ERROR!! too many interacting atoms surround atom #%i, which has %i. code needs to be optimized! %s:%i\n", i, m, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
    
    printf("\nMinimizing the potential energy:\n");
    while (totalCN > 0) {
        for (int i = 1; i <= atomnum; i ++) {
            int num = 1;
            double change[4] = {0};
            while (interaction[i].connect[num] != 0) {
                direction[1] = interaction[interaction[i].connect[num]].coordinate[1] - interaction[i].coordinate[1];
                direction[2] = interaction[interaction[i].connect[num]].coordinate[2] - interaction[i].coordinate[2];
                direction[3] = interaction[interaction[i].connect[num]].coordinate[3] - interaction[i].coordinate[3];
                
                length2 = DOT_PROD(direction, direction);
                length2 = sqrt(length2);
                
                if (length2 >= interaction[i].hiLimit[num] || length2 <= interaction[i].loLimit[num]) {
                    diff = interaction[i].length[num] - length2;
                    energy = RandomValue(0, 1) * factor * diff;
                    
                    change[1] -= energy * direction[1] / interaction[i].length[num];
                    change[2] -= energy * direction[2] / interaction[i].length[num];
                    change[3] -= energy * direction[3] / interaction[i].length[num];
                }
                num ++;
            }
            
            interaction[i].coordinate[1] += change[1];
            interaction[i].coordinate[2] += change[2];
            interaction[i].coordinate[3] += change[3];
        }
        
        //recheck the state/energy
        totalCN = 0;
        for (int i = 1; i <= atomnum; i ++) {
            int num = 1;
            while (interaction[i].connect[num] != 0) {
                if (interaction[i].connect[num] < i) {
                    num ++;
                    continue;
                }
                
                direction[1] = interaction[interaction[i].connect[num]].coordinate[1] - interaction[i].coordinate[1];
                direction[2] = interaction[interaction[i].connect[num]].coordinate[2] - interaction[i].coordinate[2];
                direction[3] = interaction[interaction[i].connect[num]].coordinate[3] - interaction[i].coordinate[3];
                
                length2 = DOT_PROD(direction, direction);
                length2 = sqrt(length2);
                
                if (length2 >= interaction[i].hiLimit[num] || length2 <= interaction[i].loLimit[num]) {
                    if (step == stepNum) {
                        printf("atoms between %8i (%s [%4i%3s]) and %8i (%s [%4i%3s]), cur dis = %.4lf, min = %.4lf, max = %.4lf\n",
                               i,
                               atom[i].property->nameOfAA,
                               atom[i].property->sequence.aminoacidNum,
                               atom[i].property->name,
                               interaction[i].connect[num],
                               atom[interaction[i].connect[num]].property->nameOfAA,
                               atom[interaction[i].connect[num]].property->sequence.aminoacidNum,
                               atom[interaction[i].connect[num]].property->name,
                               length2,
                               interaction[i].loLimit[num],
                               interaction[i].hiLimit[num]);
                    }
                    totalCN ++;
                }
                num ++;
            }
            TRANSFER_VECTOR(thisThread->raw[i]->dynamic->coordinate, interaction[i].coordinate);
        }
        
        printf("step = %5li, total conflicts = %li\n", ++step, totalCN);
        if (step % 1000 == 0) SaveGRO(&outputFile);
        if (step % 1000 == 0 && factor > 0.05) factor -= 0.01;
        if (step > stepNum) break;
    }
    
    SaveGRO(&outputFile);
    
    free(interaction);
    fclose(outputFile.file);
    printf("Done!\n\n");
}


void AtomDataCpy(struct AtomStr * dest, struct AtomStr * sour, int flag) {
    memcpy(dest->dynamic, sour->dynamic, sizeof(struct DynamicStr));
    if (flag > 0) {
        dest->property = sour->property;
    }
    
    return;
}


void ResetTarget(int *renewList, struct ThreadStr *thisThread) {
    struct AtomStr *target = NULL;
    
    for (int n = 1; n <= renewList[0]; n ++) {
        target = thisThread->listPtr[renewList[n]];
        
        target->dynamic->event.eventType = Invd_Event;
        target->dynamic->event.time = INFINIT;
        target->dynamic->event.partner = INVALID;
    }
}


void RenewCellList(struct AtomStr *newTargetAtom, struct AtomStr *oldTargetAtom) {
    int newCellIndex;
    int oldCellIndex;
    int targetNum = newTargetAtom->property->num;
    int *sourCellList = celllist;
    int nextAtom;
    
    newCellIndex = newTargetAtom->dynamic->cellIndex[3] * cellnum[1] * cellnum[2]
    + newTargetAtom->dynamic->cellIndex[2] * cellnum[1]
    + newTargetAtom->dynamic->cellIndex[1] + 1;
    
    oldCellIndex = oldTargetAtom->dynamic->cellIndex[3] * cellnum[1] * cellnum[2]
    + oldTargetAtom->dynamic->cellIndex[2] * cellnum[1]
    + oldTargetAtom->dynamic->cellIndex[1] + 1;
    
    nextAtom = oldCellIndex + atomnum;
    
    if (unlikely(sourCellList[nextAtom] == 0)) {
        printf("!!ERROR!!: link list has errors!\n");
    }
    
    //remove the target atom at the old position
    while (sourCellList[nextAtom] != targetNum) {
        nextAtom = sourCellList[nextAtom];
        
        if (unlikely(nextAtom == 0)) {
            printf("!!ERROR!!: link list has errors!\n");
        }
    }
    sourCellList[nextAtom] = sourCellList[targetNum];
    
    //add the target atom into the new position
    nextAtom = newCellIndex + atomnum;
    sourCellList[targetNum] = sourCellList[nextAtom];
    sourCellList[nextAtom] = targetNum;
}


double RandomValue(int low, int high) {
    double value;
    value = low + (double)rand()/((double)RAND_MAX/(high - low));
    
    return value;
}


void Rotation(double * old, double * new, double angle, char fixedAxis) {
    
    if (fixedAxis == 'x') {
        new[1] = old[1];
        new[2] = old[2] * cos(angle) - old[3] * sin(angle);
        new[3] = old[2] * sin(angle) + old[3] * cos(angle);
    } else if (fixedAxis == 'y') {
        new[1] = old[1] * cos(angle) + old[3] * sin(angle);
        new[2] = old[2];
        new[3] = (-1) * old[1] * sin(angle) + old[3] * cos(angle);
    } else if (fixedAxis == 'z') {
        new[1] = old[1] * cos(angle) - old[2] * sin(angle);
        new[2] = old[1] * sin(angle) + old[2] * cos(angle);
        new[3] = old[3];
    }
    
}


double Maxwell_Boltzmann_Distribution(double mean, double variance) {
    double randomNum[3] = {0};
    double choose;
    double output_Value;

repeat:
    do {
        randomNum[0] = 2 * RandomValue(0, 1) - 1;
        randomNum[1] = 2 * RandomValue(0, 1) - 1;
        randomNum[2] = randomNum[0] * randomNum[0] + randomNum[1] * randomNum[1];
    } while (randomNum[2] >= 1);
    randomNum[2] = sqrt(-2 * log(randomNum[2]) / randomNum[2]);
    
    choose = RandomValue(0, 1);
    if (choose >= 0.5) {
        output_Value = randomNum[0] * randomNum[2];
    } else {
        output_Value = randomNum[1] * randomNum[2];
    }
    
    if (output_Value == 0) {
        goto repeat;
    }
    
    return output_Value * sqrt(variance) + mean;
}


void ListRefresh(int targetatom, int * targetlist, int start, int end) {
    int i, n;
    int flag;
    
    if (start == 0 && end == 0) {
        
        FIND_DUPLICATED(targetatom, targetlist, i, flag);
        if (flag == 1) {
            n = i + 1;
            while (targetlist[n] != 0) {
                targetlist[i] = targetlist[n];
                i++;
                n++;
            }
            targetlist[i] = 0;
        }
        
    } else {
        
        FIND_DUPLICATED_IN_THE_RANGE(targetatom, targetlist, start, end, i, flag);
        if (flag == 1) {
            n = i + 1;
            while (targetlist[n] != 0 && n <= end) {
                targetlist[i] = targetlist[n];
                i++;
                n++;
            }
            targetlist[i] = 0;
        }
    }
}


double RandomVelocity(double mass) {
    double speed;
    double variance, mean;
    
    mean = 0;
    variance = targetTemperature / mass;
    
    speed = Maxwell_Boltzmann_Distribution(mean, variance);
    
    return speed;
}


void CalCOMV(int proteinNum, double *netV) {
    for (int n = protein[proteinNum].startAtomNum; n <= protein[proteinNum].endAtomNum; n ++) {
        for (int j = 0; j <= 3; j ++) {
            netV[j] += atom[n].dynamic->velocity[j] * atom[n].property->mass;
        }
    }
    
    for (int j = 0; j <= 3; j ++) {
        netV[j] /= protein[proteinNum].mass;
    }
}


double CalKinetE(struct ThreadStr *thisThread) {
    double temp[4];
    double TKE = 0;
    
    if (flow.mark < 2) {
        for (int i = 1; i <= atomnum; i ++) {
            DOT_MINUS(atom[i].dynamic->velocity, flow.constV.v, temp);
            temp[0] = DOT_PROD(temp, temp);
            
            TKE += 0.5 * atom[i].property->mass * temp[0];
        }
    } else {
        for (int i = 1; i <= numofprotein; i ++) {
            double netV[4] = {0};
            CalCOMV(i, netV);
            
            for (int n = protein[i].startAtomNum; n <= protein[i].endAtomNum; n ++) {
                DOT_MINUS(atom[n].dynamic->velocity, netV, temp);
                temp[0] = DOT_PROD(temp, temp);
                
                TKE += 0.5 * atom[n].property->mass * temp[0];
            }
        }
    }

    if (unlikely(isnan(TKE))) {
        printf("!!ERROR!!: temperature value is not valid!\n");
    }
    
    return TKE;
}


double CalPotenE(struct ThreadStr *thisThread) {
    signed int x[27]={0,1,0,0,0,1,1,1,-1,0,0,0,-1,-1,-1,0,0,1,-1,1,-1,-1,-1,1,1,1,-1};
    signed int y[27]={0,0,1,0,1,0,1,1,0,-1,0,-1,0,-1,-1,-1,1,0,0,-1,1,-1,1,-1,1,-1,1};
    signed int z[27]={0,0,0,1,1,1,0,1,0,0,-1,-1,-1,0,-1,1,-1,-1,1,0,0,1,-1,-1,-1,1,1};
    
    int atom_i, atom_j, swap;
    int pair;
    int cell_neighbor[4];
    int cellindex;
    int *atom_CellAxis;
    char eventType[20];
    double tmp;
    double coor[4];
    double speed_i[4], speed_j[4];
    double r_ij[4], v_ij[4], direction, distance2;
    double TPE = 0, accumPotential = 0;
    double positionshift[4] = {0};
    struct AtomStr **raw = thisThread->raw;
    struct AtomStr *HBPartner;
    
    for (int i = 1; i <= atomnum; i ++) {
        POINT_TO_STRUCT(atom_CellAxis, atom[i].dynamic->cellIndex);
        
        for (int n = 0; n <= 26; n ++) {
            cell_neighbor[3] = atom_CellAxis[3] + z[n];
            cell_neighbor[2] = atom_CellAxis[2] + y[n];
            cell_neighbor[1] = atom_CellAxis[1] + x[n]; //scan the neighborhood 27 subcells, include the target subcell itself
            
            for (int j = 1; j <= 3; j ++) {
                positionshift[j] = 0; //PBC position shift
                if (cell_neighbor[j] < 0) {
                    cell_neighbor[j] = cellnum[j] - 1;
                    positionshift[j] = -1 * boxDimension[j];
                } else if (cell_neighbor[j] >= cellnum[j]) {
                    cell_neighbor[j] = 0;
                    positionshift[j] = boxDimension[j];
                }
            }
            
            cellindex = cell_neighbor[3] * cellnum[1] * cellnum[2] +
                        cell_neighbor[2] * cellnum[1] +
                        cell_neighbor[1] + 1;
            
            pair = celllist[atomnum + cellindex];
            while (pair != 0) {
                if (pair > i) {
                    DOT_PLUS(atom[pair].dynamic->coordinate, positionshift, coor);
                    DOT_MINUS(atom[i].dynamic->coordinate, coor, r_ij);
                    
                    distance2 = DOT_PROD(r_ij, r_ij);
                    
                    if (distance2 <= cutoffr * cutoffr && !(connectionMap[i][pair] & BOND_CONNECT)) {
                        
                        atom_i = i;
                        atom_j = pair;
                        HBPartner = NULL;
                        
                        if (connectionMap[atom_i][atom_j] & HB_CONNECT) {
                            int type = HBModel(raw[atom_i], raw[atom_j]);
                            if (type == 1) {
                                TPE -= HBPotential.BB;
                            } else {
                                TPE -= HBPotential.BS;
                            }
                            sprintf(eventType, "HB");
                        } else if (connectionMap[atom_i][atom_j] & CONSTRAINT_CONNECT) {
                            sprintf(eventType, "constraint");
                        } else if (connectionMap[atom_i][atom_j] & NEIGHBOR_CONNECT) {
                            if (raw[raw[atom_i]->dynamic->HB.bondConnection]->dynamic->HB.neighbor == atom_j) {
                                swap = atom_i;
                                atom_i = atom_j;
                                atom_j = swap;
                            }
                            HBPartner = raw[raw[atom_j]->dynamic->HB.bondConnection];
                            sprintf(eventType, "neighbor");
                        } else {
                            sprintf(eventType, "collision");
                        }
                        
                        TRANSFER_VECTOR(speed_i, raw[atom_i]->dynamic->velocity);
                        TRANSFER_VECTOR(speed_j, raw[atom_j]->dynamic->velocity);
                        
                        DOT_MINUS(speed_i, speed_j, v_ij);
                        direction = DOT_PROD(r_ij, v_ij);
                        
                        FindPair(raw[atom_i], raw[atom_j], eventType, direction, distance2, &tmp, &tmp, &tmp, &tmp, &accumPotential, HBPartner, 0);
                        
                        if (accumPotential < INFINIT / 2) {
                            TPE += accumPotential;
                        }
                    }
                }
                
                pair = celllist[pair];
            }
        }
    }
    
    return TPE;
}


double CalSysTem(struct ThreadStr *thisThread) {
    return 2 * CalKinetE(thisThread) / (3 * atomnum); //reduced
}


double GetViscosity(double thisT) //water, range from 265 to 423 K (-8 to 150 C)
{
    double temp, convertedT;
    double thisViscosity;
    
    convertedT = thisT - 273;
    temp = (20 - convertedT) / (convertedT + 96)*(1.2378 - 1.303 / 10E3 * (20 - convertedT) + 3.06 / 10E6 * (20 - convertedT)*(20 - convertedT) + 2.55 / 10E8 * (20 - convertedT)*(20 - convertedT)*(20 - convertedT));
    temp = pow(10, temp);
    
    thisViscosity = temp * 1002 / 1.66053892;
    
    return thisViscosity;
}


void CreateGELCoordinate(int numOnAxis)
{
    int totalNum = 1;
    char buffer[50];
    FILE *outputFile;
    
    sprintf(buffer, "%s/GEL.gro", datadir);
    outputFile = fopen(buffer, "w");
    
    fprintf(outputFile, "Model\n");
    fprintf(outputFile, "%5i\n", numOnAxis * numOnAxis * numOnAxis);
    
    for (int i = 0; i < numOnAxis; ++i) {
        for (int n = 0; n < numOnAxis; ++n) {
            for (int m = 0; m < numOnAxis; ++m) {
                fprintf(outputFile, "%5iGEL%5s%5i%8.3f%8.3f%8.3f\n", totalNum, "GEL", totalNum, (m + 0.5) * 0.8, (n + 0.5) * 0.8, (i + 0.5) * 0.8);
                totalNum ++;
            }
        }
    }
    
    fprintf(outputFile, "%10.5f%10.5f%10.5f\n", numOnAxis * 0.8, numOnAxis * 0.8, numOnAxis * 0.8);
    fclose(outputFile);
}


void PBCShift (struct AtomStr *atom_i, struct AtomStr *atom_j, double *shift) {
    int cellnum_i[4], cellnum_j[4];
    
    TRANSFER_VECTOR(cellnum_i, atom_i->dynamic->cellIndex);
    TRANSFER_VECTOR(cellnum_j, atom_j->dynamic->cellIndex);
    
    for (int j = 1; j <= 3; j++) {
        shift[j] = 0; //PBC position shift
        if (cellnum_i[j] != cellnum_j[j]) {
            if (cellnum_i[j] == 0 && cellnum_j[j] == cellnum[j] - 1) {
                shift[j] = -1 * boxDimension[j];
            } else if (cellnum_i[j] == cellnum[j] - 1 && cellnum_j[j] == 0) {
                shift[j] = boxDimension[j];
            }
        }
    }
}

//print the atom list in the surrounding cells
void PrintList(list *atomList) {
    listElem *elem = listFirst(atomList);
    struct AtomStr *obj = NULL;
    
    for (int i = 1; i <= atomList->num_members; i ++) {
        obj = (struct AtomStr*)elem->obj;
        printf("atom #: %4i, time: %.4lf, cell: %i:%i:%i\n", obj->property->num, obj->dynamic->event.time, obj->dynamic->cellIndex[1], obj->dynamic->cellIndex[2], obj->dynamic->cellIndex[3]);
        elem = listNext(atomList, elem);
    }
}


void PrintPreCalList() {
    for (int i = 1; i < atomnum; i ++) {
        if (preCalList[i].eventStatus != 0) {
            printf("atom = %5i, status = %i\n", i, preCalList[i].eventStatus);
        }
    }
}


double AbsvalueVector(double * vector)
{
    double v12, v22 ,v32;
    
    v12=vector[1]*vector[1];
    v22=vector[2]*vector[2];
    v32=vector[3]*vector[3];
    
    return sqrt(v12+v22+v32);
}


void PrintBonds(int atomNum) {
    double delta, length;
    double dmin, dmax;
    struct ConstraintStr *thisBond = atom[atomNum].property->bond;
    
    while (thisBond != NULL) {
        dmin = sqrt(thisBond->dmin);
        dmax = sqrt(thisBond->dmax);
        delta  = (dmax - dmin) / (dmax + dmin);
        length = (dmax + dmin) * 0.5;
        printf("atom %-5i(%3s) -- atom %-5i(%3s) length: %lf, delta: %lf\n",
               atomNum, atom[atomNum].property->name, thisBond->connection,
               atom[thisBond->connection].property->name, length, delta);
        
        thisBond = thisBond->next;
    }
}

void PrintConstr(int atomNum) {
    struct ConstraintStr *thisConstr = atom[atomNum].property->constr;
    
    while (thisConstr != NULL) {
        printf("atom %-5i -- atom %-5i %.4lf ", atomNum, thisConstr->connection, thisConstr->dmin);
        struct StepPotenStr *thisStep = thisConstr->step;
        while (thisStep != NULL) {
            printf("%.4lf %.4lf ", thisStep->d, thisStep->e);
            thisStep = thisStep->next;
        }
        printf("%.4lf\n", thisConstr->dmax);
        
        thisConstr = thisConstr->next;
    }
}


void PrintStep(struct StepPotenStr *thisStep) {
    while (thisStep != NULL) {
        printf("%10.4lf %10.4lf (%.4lf)\n", sqrt(thisStep->d), thisStep->e, thisStep->accumulated);
        thisStep = thisStep->next;
    }
    printf("\n");
}


double CalDistance(double *p1, double *p2) {
    double distance[4];
    
    DOT_MINUS(p1, p2, distance);
    return DOT_PROD(distance, distance);
}

void CalAccumulatedPoten(struct StepPotenStr *startPoint) {
    double accumulated = 0;
    struct StepPotenStr *thisStep = NULL;
    
    thisStep = startPoint;
    while (thisStep != NULL) {
        accumulated += thisStep->e;
        thisStep = thisStep->next;
    }
    
    thisStep = startPoint;
    while (thisStep != NULL) {
        thisStep->accumulated = accumulated;
        accumulated -= thisStep->e;
        thisStep = thisStep->next;
    }
}

int EnterTunnel(double time, struct AtomStr *targetAtom) {
    double r_ij[4], r_2;
    double position_i[4], position_j[4];
    double speed[4];
    double distanceShift2, distance2;
    double radius2;
    struct AtomStr *thisTunnel = &tunlObj.tunnel;
    
    for (int n = 0; n < tunlObj.num; n ++) {
        radius2 = tunlObj.diameter[n] * 0.5;
        radius2 *= radius2;
        distanceShift2 = FindPotWellWidth(thisTunnel, targetAtom);
        
        TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
        TRANSFER_VECTOR(position_j, tunlObj.position[n]);
        TRANSFER_VECTOR(     speed, targetAtom->dynamic->velocity);
        
        position_i[1] = position_j[1];
        position_i[2] += speed[2] * time;
        position_i[3] += speed[3] * time;
        
        DOT_MINUS(position_i, position_j, r_ij);
        r_2 = DOT_PROD(r_ij, r_ij);
        
        distance2 = radius2 + distanceShift2 - 2 * sqrt(radius2 * distanceShift2);
        if (r_2 <= distance2) {
            return 1;
        }
    }
    
    return 0;
}

double FindPotWellWidth(struct AtomStr *atom_i, struct AtomStr *atom_j) {
    int type_i = atom_i->property->type;
    int type_j = atom_j->property->type;
    double width = 0;
    
    struct ConstraintStr *thisConstr = RightPair(type_i, type_j, 0);
    
    if ((width = thisConstr->dmax) > 0) {
        return width;
    } else {
        struct StepPotenStr *thisStep = thisConstr->step;
        while (thisStep->next != NULL) {
            thisStep = thisStep->next;
        }
        return thisStep->d;
    }
}

void PrintCollisionPotentialTable() {
    for (int n1 = 1; n1 <= NATOMTYPE; n1 ++) {
        for (int n2 = 1; n2 <= NATOMTYPE; n2 ++) {
            if (potentialPairCollision[n1][n2].dmin != 0) {
                printf("%2i : %2i = %6.4lf ", n1, n2, sqrt(potentialPairCollision[n1][n2].dmin));
                struct StepPotenStr *thisStep = potentialPairCollision[n1][n2].step;
                while (thisStep != NULL) {
                    printf("%6.4lf %6.4lf %6.4lf ", sqrt(thisStep->d), thisStep->e, thisStep->accumulated);
                    thisStep = thisStep->next;
                }
                printf("\n");
            }
        }
    }
}

void PrintHBPotentialTable(int HBTypeNum) {
    for (int n1 = 0; n1 <= NATOMTYPE; n1 ++) {
        for (int n2 = 0; n2 <= NATOMTYPE; n2 ++) {
            if (potentialPairHB[HBTypeNum][n1][n2].dmin != 0) {
                printf("%2i : %2i = %6.4lf ", n1, n2, sqrt(potentialPairHB[HBTypeNum][n1][n2].dmin));
                struct StepPotenStr *thisStep = potentialPairHB[HBTypeNum][n1][n2].step;
                while (thisStep != NULL) {
                    printf("%6.4lf %6.4lf %6.4lf ", sqrt(thisStep->d), thisStep->e, thisStep->accumulated);
                    thisStep = thisStep->next;
                }
                if (potentialPairHB[HBTypeNum][n1][n2].dmax > 0) {
                    printf("%6.4lf", sqrt(potentialPairHB[HBTypeNum][n1][n2].dmax));
                }
                printf("\n");
            }
        }
    }
}

void FreeConstr(struct ConstraintStr *thisConstr) {
    struct StepPotenStr *thisStep = thisConstr->step;
    struct StepPotenStr *nextStep = thisStep;
    
    while (nextStep) {
        nextStep = thisStep->next;
        free(thisStep);
        thisStep = nextStep;
    }
    
    free(thisConstr);
}

void PrintData(int *list) {
    double *coordinate;
    
    for (int i = 1; i <= list[0]; i ++) {
        int targetAtom = list[i];
        POINT_TO_STRUCT(coordinate, atom[targetAtom].dynamic->coordinate);
        
        printf("frame #%li, atom = #%i, %lf %lf %lf %.20lf\n",
               frame, targetAtom,
               coordinate[1], coordinate[2], coordinate[3],
               atom[targetAtom].dynamic->event.time);
        fflush(stdout);
    }
}
