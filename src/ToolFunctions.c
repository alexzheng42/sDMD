//
//  ToolFunctions.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

#define NEIGHBORLISTNUM 2000

void ThreadProcess(void) {

    if (REMDInfo.flag == 1) { // REMD is requested
        REMD();
    } else {
        SingleThread();
    }
    return;
}

void UpdateAtom(struct AtomStr *target, double targetTime) {
    double deltaTime = targetTime - target->dynamic->lastUpdateTime;
    if (deltaTime > 0) {
        target->dynamic->coordinate[1] += target->dynamic->velocity[1] * deltaTime;
        target->dynamic->coordinate[2] += target->dynamic->velocity[2] * deltaTime;
        target->dynamic->coordinate[3] += target->dynamic->velocity[3] * deltaTime;
        target->dynamic->lastUpdateTime = targetTime;
    }
}

void UpdateAllAtoms(double targetTime) {
    for (int i = 1; i <= atomnum; i++) {
        UpdateAtom(&atom[i], targetTime);
    }
}

int FindPair(struct AtomStr *atom1, struct AtomStr *atom2, enum EInterPairType interactionType, double direction, double distance2, double *lowerLimit, double *upperLimit, double *lowerPotential, double *upperPotential, double *accumPotential, struct AtomStr *HBPartner, int typeChange) {
    int type1 = 0, type2 = 0;
    int coreorShell = -1; // 0 -> core, 1 -> shell, -1 -> nothing
    double tolerance;
    struct ConstraintStr *thisConstr = NULL;

    if (direction < 0) {
        direction = 1; // -> <-
    } else {
        direction = -1; // <- ->
    }
    tolerance = direction * ZERO;

    switch (interactionType) {
        case InterColli:
            if (typeChange == 0) { // not change
                type1 = atom1->property->typeofAtom[0];
            } else if (typeChange == 1) { // change to _HB
                type1 = 2;
            } else if (typeChange == -1) { // change back to normal
                type1 = 1;
            } else {
                printf("!!ERROR!!: typeChange value is invalid. %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            type2 = atom2->property->typeofAtom[0];
            thisConstr = RightPair(atom1->property->typeofAtom, type1, atom2->property->typeofAtom, type2, 0);

            if (thisConstr == NULL) {
                *lowerPotential = 0;
                *upperPotential = 0;
                *accumPotential = 0;
                return -2;
            }
            break;

        case InterConstr:
            thisConstr = atom1->property->constr;
            while (thisConstr->connection != atom2->property->num) {
                thisConstr = thisConstr->next;
                if (unlikely(thisConstr == NULL)) {
                    printf("!!ERROR!!: cannot find the atom of constraint connection atom #%i for atom #%i! %s:%i\n", atom1->property->num, atom2->property->num, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
            break;

        case InterHB:
            type1 = HBModel(atom1, atom2);
            if (unlikely(type1 < 0)) {
                *lowerPotential = 0;
                *upperPotential = 0;
                *accumPotential = 0;
                return -2;
            }
            thisConstr = &potentialPairHB[type1][0][0];
            break;

        case InterNeighbor: {
            int typeNum = HBModel(atom2, HBPartner);
            if (unlikely(typeNum < 0)) {
                *lowerPotential = 0;
                *upperPotential = 0;
                *accumPotential = 0;
                return -2;
            }
            type1 = atom1->property->typeofAtom[1];
            type2 = atom2->property->typeofAtom[1];
            thisConstr = &potentialPairHB[typeNum][type1][type2];
            break;
        }

        case InterCG:
            type1 = atom1->property->typeofRes;
            type2 = atom2->property->typeofRes;
            thisConstr = &potentialPairCG[type1][type2];
            break;

        default:
            printf("!ERROR!: InteractionType assignment does not match!\n");
            exit(EXIT_FAILURE);
    }

    int flag = (thisConstr->step == NULL);
    struct StepPotenStr *thisStep = thisConstr->step;

    if (unlikely(thisConstr->dmin == 0)) {
        printf("!!ERROR!!: the potential pair between atom %i and %i is invalid! %s:%i\n", atom1->property->num, atom2->property->num, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    if (distance2 < thisConstr->dmin - ZERO) {
        *lowerLimit = INFINIT;
        *upperLimit = INFINIT;
        *accumPotential = INFINIT;
        *lowerPotential = 0;
        *upperPotential = 0;

        if (unlikely(!(interactionType == InterNeighbor && !HBIsPartner(atom2, HBPartner->property->num)) &&
                     !(interactionType == InterHB && !HBIsPartner(atom1, atom2->property->num)) &&
                     !FindElemInList(atom1->dynamic->HBNeighbor.neighborPartner, atom2->property->num, -1, -1, 0) &&
                     !(strcmp(wallDyn.mark, "no")) &&
                     !tunlObj.mark)) {
            printf("!WARNING![%li]: the distance between atom %2i(%2i%2s) and %2i(%2i%2s) is too small!\n",
                   ++warningsum,
                   atom1->property->num, atom1->property->sequence.resNum, atom1->property->name,
                   atom2->property->num, atom2->property->sequence.resNum, atom2->property->name);
            printf("         pair type: %d. the current distance is %.4lf, while the min is %.4lf.\n", interactionType, sqrt(distance2), sqrt(thisConstr->dmin));
            printf("         %s:%i\n", __FILE__, __LINE__);
        }
        return coreorShell;
    }

    *lowerLimit = thisConstr->dmin;
    *lowerPotential = 0;

    // Events land exactly on a step boundary; the direction-dependent tolerance
    // decides which side the pair is on. Potentials have only a few steps, so a
    // linear scan is faster than a lookup table (Optimization_Report.md Section 10).
    double scanLower = *lowerLimit, scanLowerPotential = *lowerPotential;
    struct StepPotenStr *firstStep = thisStep;
rescan:
    while (thisStep != NULL && thisStep->d != 0) {
        if (distance2 < thisStep->d + tolerance) {
            break;
        }
        *lowerLimit = thisStep->d;
        *lowerPotential = thisStep->e;
        thisStep = thisStep->next;
    }

    // A bonded HB sitting exactly on its outer (break) edge and moving outward belongs to the
    // last shell: an HB has no outer side. Without this, a re-prediction at that instant (e.g.
    // another event of the same atom at the same time) finds no step and cannot schedule the break.
    if ((thisStep == NULL || thisStep->d == 0) && tolerance < 0 && interactionType == InterHB &&
        distance2 < *lowerLimit + ZERO && HBIsPartner(atom1, atom2->property->num)) {
        tolerance = ZERO;
        thisStep = firstStep;
        *lowerLimit = scanLower;
        *lowerPotential = scanLowerPotential;
        goto rescan;
    }

    if (thisStep != NULL && thisStep->d != 0) {
        *upperLimit = thisStep->d;
        *upperPotential = -1 * thisStep->e;
        *accumPotential = thisStep->accumulated;

        if (thisStep == thisConstr->step) {
            coreorShell = 0;
            return coreorShell;
        } else if (interactionType == InterHB &&
                   thisStep->next == NULL && HBIsPartner(atom1, atom2->property->num)) {
            coreorShell = 1;
            return coreorShell;
        } else if (interactionType == InterNeighbor &&
                   thisStep->next == NULL && HBIsPartner(atom2, HBPartner->property->num)) {
            *upperPotential = 0;
            coreorShell = 1;
            return coreorShell;
        }

        return coreorShell;
    }

    if (unlikely(interactionType == InterHB && HBIsPartner(atom1, atom2->property->num))) {
        printf("!WARNING![%li]: HB length between atom %i and %i is too large!\n",
               ++warningsum,
               atom1->property->sequence.atomNum,
               atom2->property->sequence.atomNum);
        printf("         the current is %.4lf, while the max is %.4lf\n", sqrt(distance2), sqrt(*upperLimit));
        printf("         %s:%i\n", __FILE__, __LINE__);
    }

    if (thisConstr->dmax > 0) {
        if (distance2 < thisConstr->dmax + ZERO) {
            *upperLimit = thisConstr->dmax;
            *upperPotential = 0;
            *accumPotential = 0;
            coreorShell = 1;
            if (flag == 1 && direction == 1) {
                coreorShell = 0;
            }
            return coreorShell;
        } else {
            *lowerLimit = INFINIT;
            *upperLimit = INFINIT;
            *accumPotential = INFINIT;
            *lowerPotential = 0;
            *upperPotential = 0;

            if (unlikely(!(interactionType == InterNeighbor && !HBIsPartner(atom2, HBPartner->property->num)) &&
                         !(interactionType == InterHB && !HBIsPartner(atom1, atom2->property->num)))) {
                printf("!WARNING![%li]: the distance between atom %2i(%2i%2s) and %2i(%2i%2s) is too large!\n",
                       ++warningsum,
                       atom1->property->num, atom1->property->sequence.resNum, atom1->property->name,
                       atom2->property->num, atom2->property->sequence.resNum, atom2->property->name);
                printf("         pair type: %d. the current distance is %.4lf, while the min is %.4lf.\n", interactionType, sqrt(distance2), sqrt(thisConstr->dmax));
                printf("         %s:%i\n", __FILE__, __LINE__);
            }
        }
    } else {
        *upperLimit = INFINIT;
        *upperPotential = INFINIT;
        *accumPotential = INFINIT;
        coreorShell = 1;
        if (flag == 1 && direction == 1) {
            coreorShell = 0;
        }
        return coreorShell;
    }

    return coreorShell;
}

/*
 after HB forming, the type changed atoms may have no
 cooresponding potential pair for collision event,
 we need to use their original types to find the pair
 */
struct ConstraintStr *RightPair(int *atom1_typeList, int type1, int *atom2_typeList, int type2, int flag) {

    if (potentialPairCollision[atom1_typeList[type1]][atom2_typeList[type2]].dmin > 0) {
        return &potentialPairCollision[atom1_typeList[type1]][atom2_typeList[type2]];
    } else if (flag) {
        return NULL;
    }

    struct ConstraintStr *targetPair = NULL;
    if ((targetPair = RightPair(atom1_typeList, 1, atom2_typeList, type2, 1)) != NULL ||
        (targetPair = RightPair(atom1_typeList, type1, atom2_typeList, 1, 1)) != NULL ||
        (targetPair = RightPair(atom1_typeList, 1, atom2_typeList, 1, 1)) != NULL) {
    }

    return targetPair;
}

/*
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
*/

void DisplayTime(char *timer) {
    time_t timetodisplay;
    struct tm *tm_info;

    time(&timetodisplay);
    tm_info = localtime(&timetodisplay);

    strftime(timer, 30, "%Y_%m_%d_%H_%M_%S", tm_info);
}

void SDEnergyMin(long int stepNum, struct AtomListStr *thisList) {
    struct interStr {
        int connect[NEIGHBORLISTNUM];
        double coordinate[4];
        double hiLimit[NEIGHBORLISTNUM];
        double loLimit[NEIGHBORLISTNUM];
        double length[NEIGHBORLISTNUM];
    };

    char directory[PATH_LEN];
    int cell_neighbor[4], cellindex, atom_CellAxis[2][4];
    int num, pair;
    long step = 0;
    long totalCN = INFINIT * INFINIT;
    double factor = 0.5;
    double diff, energy, length2;
    double direction[4], newCoor[4];
    double positionshift[4] = {0};
    struct interStr *interaction;
    struct FileStr outputFile;

    PBC("full", NULL, &atomList);
    LinkList("full", NULL, NULL, &atomList);

    interaction = (struct interStr *)calloc(atomnum + 1, sizeof(struct interStr));
    sprintf(directory, "%s/MiniEnergy%s.gro", datadir, REMDInfo.REMD_ExtraName);
    outputFile.file = fopen(directory, "w");

    for (int i = 1; i <= atomnum; ++i) {
        TRANSFER_VECTOR(interaction[i].coordinate, thisList->ptr[i]->dynamic->coordinate);

        // prefixed atoms will not be optimized
        if (FindStrInList(thisList->ptr[i]->property->extraProperty, "_CONS", 0)) {
            continue;
        }

        // transfer bond information
        num = 0;
        struct ConstraintStr *thisBond = thisList->ptr[i]->property->bond;
        while (thisBond != NULL) {
            num++;
            interaction[i].connect[num] = thisBond->connection;
            interaction[i].hiLimit[num] = sqrt(thisBond->dmax);
            interaction[i].loLimit[num] = sqrt(thisBond->dmin);
            interaction[i].length[num] = 0.5 * (interaction[i].hiLimit[num] + interaction[i].loLimit[num]);
            thisBond = thisBond->next;
        }

        // transfer constrains information
        struct ConstraintStr *thisConstr = thisList->ptr[i]->property->constr;
        while (thisConstr != NULL) {
            num++;
            interaction[i].connect[num] = thisConstr->connection;
            interaction[i].hiLimit[num] = sqrt(thisConstr->dmax);
            interaction[i].loLimit[num] = sqrt(thisConstr->dmin);
            interaction[i].length[num] = 0.5 * (interaction[i].hiLimit[num] + interaction[i].loLimit[num]);
            thisConstr = thisConstr->next;
        }

        // transfer interaction potential information
        TRANSFER_VECTOR(atom_CellAxis[0], thisList->ptr[i]->dynamic->cellIndex)
        for (int n = 0; n < atomList.stencilNum; n++) {
            cell_neighbor[3] = atom_CellAxis[0][3] + atomList.z[n];
            cell_neighbor[2] = atom_CellAxis[0][2] + atomList.y[n];
            cell_neighbor[1] = atom_CellAxis[0][1] + atomList.x[n]; // scan the neighborhood 27 subcells, include the target subcell itself

            for (int j = 1; j <= 3; j++) {
                positionshift[j] = 0; // PBC position shift
                if (cell_neighbor[j] < 0) {
                    cell_neighbor[j] += cellnum[j];
                    positionshift[j] = -1 * boxDimension[j];
                } else if (cell_neighbor[j] >= cellnum[j]) {
                    cell_neighbor[j] -= cellnum[j];
                    positionshift[j] = boxDimension[j];
                }
            }

            cellindex = cell_neighbor[3] * cellnum[1] * cellnum[2] +
                        cell_neighbor[2] * cellnum[1] +
                        cell_neighbor[1] + 1;

            pair = celllist[atomnum + cellindex];
            while (pair != 0) {
                if (pair != i) {
                    DOT_PLUS(thisList->ptr[pair]->dynamic->coordinate, positionshift, newCoor);
                    DOT_MINUS(interaction[i].coordinate, newCoor, direction);
                    length2 = DOT_PROD(direction, direction);

                    if (length2 <= cutoffr * cutoffr &&
                        !(connectionMap[i][pair] & BOND_CONNECT) &&
                        !(connectionMap[i][pair] & CONSTRAINT_CONNECT)) {
                        if (!FindElemInList(interaction[i].connect, pair, -1, -1, 0)) {
                            num++;
                            interaction[i].connect[num] = pair;    // all the surrounding atoms
                            interaction[i].hiLimit[num] = INFINIT; // interaction potential only has low limit

                            int type1 = thisList->ptr[i]->property->typeofAtom[thisList->ptr[i]->property->typeofAtom[0]];
                            int type2 = thisList->ptr[pair]->property->typeofAtom[thisList->ptr[pair]->property->typeofAtom[0]];

                            if (type1 > type2) {
                                int s = type2;
                                type2 = type1;
                                type1 = s;
                            }

                            interaction[i].loLimit[num] = sqrt(potentialPairCollision[type1][type2].dmin);
                            interaction[i].length[num] = interaction[i].loLimit[num] + 0.05;
                        }
                    }
                }
                pair = celllist[pair];
            }
        }

        if (num > (NEIGHBORLISTNUM - 2)) {
            printf("!!ERROR!! too many interacting atoms surround atom #%i, which has %i. code needs to be optimized! %s:%i\n", i, num, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }

    printf("\nMinimizing the potential energy:\n");
    while (totalCN > 0) {
        for (int i = 1; i <= atomnum; i++) {
            double change[4] = {0};

            // prefixed atoms will not be optimized
            if (FindStrInList(thisList->ptr[i]->property->extraProperty, "_CONS", 0)) {
                continue;
            }

            num = 1;
            while (interaction[i].connect[num] != 0) {
                int n = interaction[i].connect[num];
                TRANSFER_VECTOR(atom_CellAxis[0], thisList->ptr[i]->dynamic->cellIndex);
                TRANSFER_VECTOR(atom_CellAxis[1], thisList->ptr[n]->dynamic->cellIndex);

                CellPBCShift(atom_CellAxis[0], atom_CellAxis[1], positionshift);

                DOT_PLUS(interaction[n].coordinate, positionshift, newCoor);
                DOT_MINUS(interaction[i].coordinate, newCoor, direction);
                length2 = DOT_PROD(direction, direction);
                length2 = sqrt(length2);

                if (length2 >= interaction[i].hiLimit[num] || length2 <= interaction[i].loLimit[num]) {

#ifdef DEBUG_PRINTF
                    printf("atom %-5i(%3s) -- atom %-5i(%3s) orig_length: %lf, min: %lf, max: %lf, curr_length: %lf\n",
                           i, atom[i].property->name,
                           interaction[i].connect[num], atom[interaction[i].connect[num]].property->name,
                           interaction[i].length[num],
                           interaction[i].loLimit[num],
                           interaction[i].hiLimit[num],
                           length2);
#endif

                    diff = length2 - interaction[i].length[num];
                    energy = RandomValue(0, 1) * factor * diff;

                    tmpDouble = (interaction[i].length[num] < 1.0 ? 1.0 : interaction[i].length[num]);
                    change[1] -= energy * direction[1] / tmpDouble;
                    change[2] -= energy * direction[2] / tmpDouble;
                    change[3] -= energy * direction[3] / tmpDouble;
                }
                num++;
            }

            interaction[i].coordinate[1] += change[1];
            interaction[i].coordinate[2] += change[2];
            interaction[i].coordinate[3] += change[3];

            for (int dim = 1; dim <= 3; dim++) {
                if (isnan(interaction[i].coordinate[dim]) ||
                    fabs(interaction[i].coordinate[dim]) > INFINIT) {
                    printf("!!ERROR!!: energy minimization failed! data overflow! Please check the system structure.\n");
                    printf("           %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
        }

        // recheck the state/energy
        totalCN = 0;
        for (int i = 1; i <= atomnum; i++) {
            num = 1;
            while (interaction[i].connect[num] != 0) {
                int n = interaction[i].connect[num];
                TRANSFER_VECTOR(atom_CellAxis[0], thisList->ptr[i]->dynamic->cellIndex);
                TRANSFER_VECTOR(atom_CellAxis[1], thisList->ptr[n]->dynamic->cellIndex);

                CellPBCShift(atom_CellAxis[0], atom_CellAxis[1], positionshift);

                DOT_PLUS(interaction[n].coordinate, positionshift, newCoor);
                DOT_MINUS(interaction[i].coordinate, newCoor, direction);
                length2 = DOT_PROD(direction, direction);
                length2 = sqrt(length2);

                if (length2 >= interaction[i].hiLimit[num] || length2 <= interaction[i].loLimit[num]) {
                    if (step == stepNum) {
                        printf("atoms between %8i (%s [%4i%3s]) and %8i (%s [%4i%3s]), cur dis = %.4lf, min = %.4lf, max = %.4lf\n",
                               i,
                               atom[i].property->resName,
                               atom[i].property->sequence.resNum,
                               atom[i].property->name,
                               interaction[i].connect[num],
                               atom[interaction[i].connect[num]].property->resName,
                               atom[interaction[i].connect[num]].property->sequence.resNum,
                               atom[interaction[i].connect[num]].property->name,
                               length2,
                               interaction[i].loLimit[num],
                               interaction[i].hiLimit[num]);
                    }
                    totalCN++;
                }
                num++;
            }
            TRANSFER_VECTOR(thisList->ptr[i]->dynamic->coordinate, interaction[i].coordinate);
        }

        printf("step = %5li, total conflicts = %li\n", ++step, totalCN);
        if (step % 1000 == 0)
            SaveGRO(&outputFile);
        if (step % 1000 == 0 && factor > 0.05)
            factor -= 0.01;
        if (step > stepNum)
            break;
    }

    SaveGRO(&outputFile);

    free(interaction);
    free(celllist);
    celllist = NULL; // so that the next LinkList call can be valid.
    fclose(outputFile.file);
    printf("Done!\n\n");

    return;
}

void AtomDataCpy(struct AtomStr *dest, struct AtomStr *sour, int flag) {
    memcpy(dest->dynamic, sour->dynamic, sizeof(struct DynamicStr));
    if (flag > 0 || FindStrInList(sour->property->extraProperty, "_CONS", 0)) {
        dest->property = sour->property;
    }

    return;
}

void ResetTarget(int *renewList, struct AtomListStr *thisList) {
    struct AtomStr *target = NULL;

    for (int n = 1; renewList[n]; n++) {
        target = thisList->ptr[renewList[n]];

        target->dynamic->event.eventType = Invd_Event;
        target->dynamic->event.time = INFINIT;
        target->dynamic->event.partner = INVALID;
        ClearBackupEvents(renewList[n]);
    }
}

double RandomValue(int low, int high) {
    double value;
    value = low + (double)rand() / ((double)RAND_MAX / (high - low));

    return value;
}

void Rotation(double *old, double *new, double angle, char fixedAxis) {

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

void ListRefresh(int target, int *targetlist, int start, int end) {
    int i, n;

    if (start == 0 && end == 0) {
        if ((i = FindElemInList(targetlist, target, -1, -1, 0))) {
            n = i + 1;
            while (targetlist[n] != 0) {
                targetlist[i] = targetlist[n];
                i++;
                n++;
            }
            targetlist[i] = 0;
        }

    } else {
        if ((i = FindElemInList(targetlist, target, -1, -1, 0))) {
            n = i + 1;
            while (targetlist[n] != 0 && n <= end) {
                targetlist[i] = targetlist[n];
                i++;
                n++;
            }
            targetlist[i] = 0;
        }
    }
    
    return;
}

double RandomVelocity(double mass) {
    double speed;
    double variance, mean;

    mean = 0;
    variance = targetTemperature / mass;

    speed = Maxwell_Boltzmann_Distribution(mean, variance);

    return speed;
}

void CalCOMV(int molNum, double *netV) {
    for (int n = molecule[molNum].startAtomNum; n <= molecule[molNum].endAtomNum; n++) {
        for (int j = 0; j <= 3; j++) {
            netV[j] += atom[n].dynamic->velocity[j] * atom[n].property->mass;
        }
    }

    for (int j = 0; j <= 3; j++) {
        netV[j] /= molecule[molNum].mass;
    }
}

double CalKinetE(void) {
    double temp[4];
    double TKE = 0;

    if (flow.mark < 2) {
        for (int i = 1; i <= atomnum; i++) {
            DOT_MINUS(atom[i].dynamic->velocity, flow.constV.v, temp);
            temp[0] = DOT_PROD(temp, temp);

            TKE += 0.5 * atom[i].property->mass * temp[0];
        }
    } else {
        for (int i = 1; i <= numofMol; i++) {
            double netV[4] = {0};
            CalCOMV(i, netV);

            for (int n = molecule[i].startAtomNum; n <= molecule[i].endAtomNum; n++) {
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

double CalPotenE(void) {

    int atom_i, atom_j, swap;
    int pair;
    int cell_neighbor[4];
    int cellindex;
    int *atom_CellAxis;
    enum EInterPairType interactionType;
    double tmp;
    double coor[4];
    double speed_i[4], speed_j[4];
    double r_ij[4], v_ij[4], direction, distance2;
    double TPE = 0, accumPotential = 0;
    double positionshift[4] = {0};
    struct AtomStr *HBPartner;

    for (int i = 1; i <= atomnum; i++) {
        POINT_TO_STRUCT(atom_CellAxis, atom[i].dynamic->cellIndex);

        for (int n = 0; n < atomList.stencilNum; n++) {
            cell_neighbor[3] = atom_CellAxis[3] + atomList.z[n];
            cell_neighbor[2] = atom_CellAxis[2] + atomList.y[n];
            cell_neighbor[1] = atom_CellAxis[1] + atomList.x[n]; // scan the neighborhood 27 subcells, include the target subcell itself

            for (int j = 1; j <= 3; j++) {
                positionshift[j] = 0; // PBC position shift
                if (cell_neighbor[j] < 0) {
                    cell_neighbor[j] += cellnum[j];
                    positionshift[j] = -1 * boxDimension[j];
                } else if (cell_neighbor[j] >= cellnum[j]) {
                    cell_neighbor[j] -= cellnum[j];
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
                            interactionType = InterHB;
                            TPE -= HBBarrier(&atom[atom_i], &atom[atom_j]);
                        } else if (connectionMap[atom_i][atom_j] & CONSTRAINT_CONNECT) {
                            interactionType = InterConstr;
                        } else if (connectionMap[atom_i][atom_j] & NEIGHBOR_CONNECT) {
                            int hostPartner;
                            if ((hostPartner = HBNeighborHostPartner(&atom[atom_i], atom_j)) > 0) {
                                swap = atom_i;
                                atom_i = atom_j;
                                atom_j = swap;
                            } else if ((hostPartner = HBNeighborHostPartner(&atom[atom_j], atom_i)) == 0) {
                                hostPartner = atom[atom_j].dynamic->HB.bondConnection[1];
                            }
                            HBPartner = &atom[hostPartner];
                            interactionType = InterNeighbor;
                        } else {
                            interactionType = InterColli;
                        }

                        TRANSFER_VECTOR(speed_i, atom[atom_i].dynamic->velocity);
                        TRANSFER_VECTOR(speed_j, atom[atom_j].dynamic->velocity);

                        DOT_MINUS(speed_i, speed_j, v_ij);
                        direction = DOT_PROD(r_ij, v_ij);

                        FindPair(&atom[atom_i], &atom[atom_j], interactionType, direction, distance2, &tmp, &tmp, &tmp, &tmp, &accumPotential, HBPartner, 0);

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

double CalSysTem(void) {
    return 2 * CalKinetE() / (3 * atomnum); // reduced
}

double GetViscosity(double thisT) // water, range from 265 to 423 K (-8 to 150 C)
{
    double temp, convertedT;
    double thisViscosity;

    convertedT = thisT - 273;
    temp = (20 - convertedT) / (convertedT + 96) * (1.2378 - 1.303 / 10E3 * (20 - convertedT) + 3.06 / 10E6 * (20 - convertedT) * (20 - convertedT) + 2.55 / 10E8 * (20 - convertedT) * (20 - convertedT) * (20 - convertedT));
    temp = pow(10.0, temp);

    thisViscosity = temp * 1002 / 1.66053892;

    return thisViscosity;
}

void CreateGELCoordinate(int numOnAxis) {
    int totalNum = 1;
    char buffer[PATH_LEN];
    FILE *outputFile;

    sprintf(buffer, "%s/GEL.gro", datadir);
    outputFile = fopen(buffer, "w");

    fprintf(outputFile, "Model\n");
    fprintf(outputFile, "%5i\n", numOnAxis * numOnAxis * numOnAxis);

    for (int i = 0; i < numOnAxis; ++i) {
        for (int n = 0; n < numOnAxis; ++n) {
            for (int m = 0; m < numOnAxis; ++m) {
                fprintf(outputFile, "%5iGEL%5s%5i%8.3f%8.3f%8.3f\n", totalNum, "GEL", totalNum, (m + 0.5) * 0.8, (n + 0.5) * 0.8, (i + 0.5) * 0.8);
                totalNum++;
            }
        }
    }

    fprintf(outputFile, "%10.5f%10.5f%10.5f\n", numOnAxis * 0.8, numOnAxis * 0.8, numOnAxis * 0.8);
    fclose(outputFile);
}

// PBC shift to add to atom j's position when seen from atom i, from their cell
// indices (minimum image). Exact for any pair within the scan stencil, because
// cellnum >= 2 * reach + 1; for the 27-cell stencil it is the classic
// "cell 0 next to cell cellnum - 1" rule.
void CellPBCShift(int *cell_i, int *cell_j, double *shift) {
    for (int j = 1; j <= 3; j++) {
        int d = cell_j[j] - cell_i[j];
        if (2 * d > cellnum[j]) {
            shift[j] = -1 * boxDimension[j];
        } else if (-2 * d > cellnum[j]) {
            shift[j] = boxDimension[j];
        } else {
            shift[j] = 0;
        }
    }
}

void PBCShift(struct AtomStr *atom_i, struct AtomStr *atom_j, double *shift) {
    CellPBCShift(atom_i->dynamic->cellIndex, atom_j->dynamic->cellIndex, shift);
}

// print the atom list in the surrounding cells
/*
void PrintList(list *atomList) {
    listElem *elem = listFirst(atomList);
    struct AtomStr *obj = NULL;

    for (int i = 1; i <= atomList->num_members; i ++) {
        obj = (struct AtomStr*)elem->obj;
        printf("atom #: %4i, time: %.4lf, cell: %i:%i:%i\n", obj->property->num, obj->dynamic->event.time, obj->dynamic->cellIndex[1], obj->dynamic->cellIndex[2], obj->dynamic->cellIndex[3]);
        elem = listNext(atomList, elem);
    }
}
 */

double AbsvalueVector(double *vector) {
    double v12, v22, v32;

    v12 = vector[1] * vector[1];
    v22 = vector[2] * vector[2];
    v32 = vector[3] * vector[3];

    return sqrt(v12 + v22 + v32);
}

void PrintBonds(int atomNum) {
    double delta, length;
    double dmin, dmax;
    struct ConstraintStr *thisBond = atom[atomNum].property->bond;

    while (thisBond != NULL) {
        dmin = sqrt(thisBond->dmin);
        dmax = sqrt(thisBond->dmax);
        delta = (dmax - dmin) / (dmax + dmin);
        length = (dmax + dmin) * 0.5;
        printf("atom %-5i(%3s) -- atom %-5i(%3s) length: %lf, delta: %lf, min: %lf, max: %lf\n",
               atomNum, atom[atomNum].property->name, thisBond->connection,
               atom[thisBond->connection].property->name, length, delta, dmin, dmax);

        thisBond = thisBond->next;
    }
}

void PrintConstr(int atomNum) {
    struct ConstraintStr *thisConstr = atom[atomNum].property->constr;

    while (thisConstr != NULL) {
        printf("atom %-5i(%-6s) -- atom %-5i(%-6s)\n",
               atomNum, atom[atomNum].property->name,
               thisConstr->connection, atom[thisConstr->connection].property->name);
        printf("           %8s %8s %8s\n", "dis", "acu", "dif");
        printf("     dmin: %8.4lf\n", sqrt(thisConstr->dmin));

        struct StepPotenStr *thisStep = thisConstr->step;
        while (thisStep != NULL) {
            printf("           %8.4lf %8.4lf %8.4lf\n", sqrt(thisStep->d), thisStep->accumulated, thisStep->e);
            thisStep = thisStep->next;
        }
        printf("     dmax: %8.4lf\n", sqrt(thisConstr->dmax));

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
    /*
     the potential table is read from the right to the left.
     the most right potential (start potential) = delta_E_most_right
     the other potentials = E_current + delta_E
     */

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

    return;
}

int EnterTunnel(double time, struct AtomStr *target) {
    double r_ij[4], r_2;
    double position_i[4], position_j[4];
    double speed[4];
    double distanceShift2, distance2;
    double radius2;
    struct AtomStr *thisTunnel = &tunlObj.tunnel;

    for (int n = 0; n < tunlObj.num; n++) {
        radius2 = tunlObj.diameter[n] * 0.5;
        radius2 *= radius2;
        distanceShift2 = FindPotWellWidth(thisTunnel, target);

        TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
        TRANSFER_VECTOR(position_j, tunlObj.position[n]);
        TRANSFER_VECTOR(speed, target->dynamic->velocity);

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
    double width = 0;

    struct ConstraintStr *thisConstr = RightPair(atom_i->property->typeofAtom,
                                                 atom_i->property->typeofAtom[0],
                                                 atom_j->property->typeofAtom,
                                                 atom_j->property->typeofAtom[0],
                                                 0);

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

void PrintCollisionPotentialTable(void) {
    for (int n1 = 1; n1 <= FFInfo.bead[0].val[0]; n1++) {
        for (int n2 = 1; n2 <= FFInfo.bead[0].val[0]; n2++) {
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
    for (int n1 = 0; n1 <= FFInfo.bead[0].val[0]; n1++) {
        for (int n2 = 0; n2 <= FFInfo.bead[0].val[0]; n2++) {
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

void PrintListData(int *thisList) {
    double *coordinate;

    for (int i = 1; i <= thisList[0]; i++) {
        int target = thisList[i];
        POINT_TO_STRUCT(coordinate, atom[target].dynamic->coordinate);

        printf("frame #%li, atom = #%i, %lf %lf %lf %.20lf\n",
               frame, target,
               coordinate[1], coordinate[2], coordinate[3],
               atom[target].dynamic->event.time);
        fflush(stdout);
    }
}

void PrintHBData(void) {
    for (int i = 1; i <= atomnum; i++) {
        for (int n = 1; n <= atomnum; n++) {
            if (HBIsPartner(&atom[i], n) ||
                HBIsPartner(&atom[n], i) ||
                connectionMap[i][n] > 2 ||
                connectionMap[n][i] > 2) {

                printf("#%4i ", i);
                printf("HB: %4i %4i, neighbors:", atom[i].dynamic->HB.bondConnection[1], atom[i].dynamic->HB.bondConnection[2]);
                for (int m = 1; m <= HB_NEIGHBOR_MAX; m++) {
                    printf(" %4i,", atom[i].dynamic->HBNeighbor.neighborPartner[m]);
                }
                printf(" record:");
                for (int m = 1; m <= atomnum; m++) {
                    if (connectionMap[i][m] > 2) {
                        printf(" %5i-%5i[%i],", i, m, connectionMap[i][m]);
                    }
                }
                printf("\n");
            }
        }
    }

    return;
}

int FindElemInList(int *list, int value, int start, int end, int add) {
    if (end >= 0 && start >= 0 && end >= start) {
        for (int i = start; i <= end; i++) {
            if (list[i] == value) {
                return i;
            } else if (list[i] == 0) {
                if (add)
                    list[i] = value;
                return 0;
            }
        }
    } else {
        if (add) {
            printf("!!ERROR!!: unbounded search does not support \"add\"!\n");
            printf("%s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        int i = 0;
        while (list[++i]) {
            if (list[i] == value) {
                return i;
            }
        }
    }
    return 0;
}


double CalExponential(double value, int times) {
    double res = 1.0;

    for (int i = 0; i < times; i++) {
        res *= value;
    }

    return res;
}


// new -> old per-atom layout (first HB partner only), for the files sAnalysis reads
void DynamicToLegacy(const struct DynamicStr *dynamic, struct DynamicStrLegacy *old) {
    memset(old, 0, sizeof(struct DynamicStrLegacy));
    memcpy(old->cellIndex, dynamic->cellIndex, sizeof(old->cellIndex));
    memcpy(old->coordinate, dynamic->coordinate, sizeof(old->coordinate));
    memcpy(old->velocity, dynamic->velocity, sizeof(old->velocity));
    old->lastUpdateTime = dynamic->lastUpdateTime;
    old->HB.bondConnection = dynamic->HB.bondConnection[0] ? dynamic->HB.bondConnection[1] : 0;
    memcpy(old->HB.neighbor, dynamic->HB.neighbor, sizeof(old->HB.neighbor));
    old->HB.role = dynamic->HB.role;
    old->HB.interactionType = dynamic->HB.interactionType;
    memcpy(old->HBNeighbor.neighborPartner, dynamic->HBNeighbor.neighborPartner, sizeof(old->HBNeighbor.neighborPartner));
    old->HBNeighbor.partnerNum = dynamic->HBNeighbor.partnerNum;
    old->HBNeighbor.neighborStatus = dynamic->HBNeighbor.neighborStatus;
    memcpy(old->HBNeighbor.removedAtom, dynamic->HBNeighbor.removedAtom, sizeof(old->HBNeighbor.removedAtom));
    old->event = dynamic->event;
}
