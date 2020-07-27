//
//  TimePrediction.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

#define HBPAIR_AD(atom1, atom2) \
(atom1->dynamic->HB.role == 'A' && atom2->dynamic->HB.role == 'D') || \
(atom1->dynamic->HB.role == 'D' && atom2->dynamic->HB.role == 'A')

#define HBPAIR_SD(atom1, atom2) \
(atom1->dynamic->HB.role == 'S' && atom2->property->typeofAtom == 2) || \
(atom2->dynamic->HB.role == 'S' && atom1->property->typeofAtom == 2)

#define HBPAIR_SS(atom1, atom2) \
(atom1->dynamic->HB.role == 'S' && atom2->dynamic->HB.role == 'S')


struct ParameterStr {
    int coreorShell;
    int s_time;
    double r_temp[4], v_ij[4], r_ij[4], b_ij, v_2, r_2, d_ij2;
    double positionshift[4];
    double lowerPotential, upperPotential;
    double shortlimit2, longlimit2;
    double calibratedcutoffr;
    double position_i[4], speed_i[4];
    double position_j[4], speed_j[4];
};

void InitializeEvent(struct InteractionEventStr *event);
void CalculateParameters(struct AtomStr *atom_i, struct AtomStr *atom_j, struct ParameterStr *parameters);
void CalculateVCOM(int startNum, int endNum, struct AtomStr *atomList, double *vCOM);
int CheckDuplication(struct AtomStr *atom_i, struct AtomStr *atom_j, enum EEventType type);
int CheckPairMatch(struct AtomStr *atom_i, struct AtomStr *atom_j, struct AtomListStr *thisList);
double CalculateDisc(struct ParameterStr *parameters);
double CalculateTime(struct ParameterStr *parameters);


void CollisionTime(struct AtomStr *collision_i, struct AtomListStr *thisList) {
    int target_i, target_j;
    int cell_neighbor[4], cellIndex;
    int *sourCellList = celllist;
    struct AtomStr *collision_j;
    struct InteractionEventStr event;
    struct ParameterStr parameters;
    
    InitializeEvent(&event);
    parameters.calibratedcutoffr = cutoffr * EXTEND_RATIO; //this radius should be less than the length of subcell
    parameters.calibratedcutoffr *= parameters.calibratedcutoffr;
    
    target_i = collision_i->property->num;
    TRANSFER_VECTOR(parameters.position_i, collision_i->dynamic->coordinate); //easier to use coordinate and velocity data
    TRANSFER_VECTOR(parameters.speed_i,    collision_i->dynamic->velocity);
    
    for (int i = 0; i < 27; ++i) {
        //scan the neighborhood 27 subcells, include the target subcell itself
        cell_neighbor[3] = collision_i->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = collision_i->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = collision_i->dynamic->cellIndex[1] + thisList->x[i];
        
        for (int j = 1; j <= 3; j++) {
            if (unlikely(cell_neighbor[j] < 0)) {
                cell_neighbor[j] = cellnum[j] - 1;
            } else if (unlikely(cell_neighbor[j] >= cellnum[j])) {
                cell_neighbor[j] = 0;
            }
        }
        
        cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2]
        + cell_neighbor[2] * cellnum[1]
        + cell_neighbor[1] + 1;
        
        target_j = sourCellList[atomnum + cellIndex];
        while (target_j) {
            collision_j = thisList->ptr[target_j];
            
            if (target_i != target_j) {
                //---------------------------------
                //check if there is any duplicate case or if this pair has been already registered
                //new code allows the target atom and its partner both have the same event information
                //which will reduce the complexity of the list maintenance
                if (CheckDuplication(collision_i, collision_j, Coli_Event)) {
                    target_j = sourCellList[target_j];
                    continue;
                }
                //---------------------------------
                
                //---------------------------------
                //if the two target beads are connected by bonds, or in HB relationship, skip them
                if (connectionMap[target_i][target_j] & BOND_CONNECT ||
                    collision_i->dynamic->HB.bondConnection == target_j ||
                    FindElemInList(collision_i->dynamic->HBNeighbor.neighborPartner, target_j, -1, -1, 0)) {
                    target_j = sourCellList[target_j];
                    continue;
                }
                //---------------------------------
                
                TRANSFER_VECTOR(parameters.position_j, collision_j->dynamic->coordinate); //easier to use coordinate and velocity data
                TRANSFER_VECTOR(parameters.speed_j,    collision_j->dynamic->velocity);
                CalculateParameters(collision_i, collision_j, &parameters);
                
                if (parameters.r_2 <= parameters.calibratedcutoffr) { //if the distance between two atoms <= scan radius
                    
                    if (connectionMap[target_i][target_j] & CONSTRAINT_CONNECT) {
                        parameters.coreorShell = FindPair(collision_i, collision_j, "constraint",
                                                          parameters.b_ij, parameters.r_2,
                                                          &parameters.shortlimit2,    &parameters.longlimit2,
                                                          &parameters.lowerPotential, &parameters.upperPotential,
                                                          &tmpDouble, NULL, 0);
                    } else {
                        parameters.coreorShell = FindPair(collision_i, collision_j, "collision",
                                                          parameters.b_ij, parameters.r_2,
                                                          &parameters.shortlimit2, &parameters.longlimit2,
                                                          &parameters.lowerPotential, &parameters.upperPotential,
                                                          &tmpDouble, NULL, 0);
                        if (parameters.lowerPotential == 0 && parameters.upperPotential == 0) {
                            target_j = sourCellList[target_j];
                            continue;
                        }
                    }
                    
                    //judge which type of collisions will happen
                    //core hit? well capture? well bounce? well escape?
                    //Supplement to Introductory Chemical Engineering Thermodynamics, 2nd ED.
                    if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) { //approaching and will collide
                        
                        parameters.s_time    = -1;
                        parameters.d_ij2     = parameters.shortlimit2;
                        event.potential = parameters.lowerPotential;
                        
                        if (parameters.coreorShell == 0) { //inner well contact
                            event.subEventType = CoreColli;
                        } else { //well contact
                            event.subEventType = TBDCaptr;
                        }
                        
                    } else { //move further away
                        
                        parameters.s_time    = 1;
                        parameters.d_ij2     = parameters.longlimit2;
                        event.potential = parameters.upperPotential;
                        
                        if (strncmp(Methodtype, "Ding", 1) == 0) {
                            if (connectionMap[target_i][target_j] & CONSTRAINT_CONNECT) {
                                if (parameters.coreorShell == 1) { //most outer well contact
                                    event.subEventType = WellBounc;
                                } else { //well contact (pull)
                                    event.subEventType = TBDEscap;
                                }
                            } else {
                                if (parameters.coreorShell == 1) {
                                    target_j = sourCellList[target_j];
                                    continue;
                                } else {
                                    event.subEventType = TBDEscap;
                                }
                            }
                        }
                    }
                    
                    event.time = CalculateTime(&parameters);
                    
                    if (unlikely(event.time < 0)) {
                        printf("!!ERROR!!: collision time is less than zero bewteen atoms %i and %i\n", target_i, target_j); //for debug and error check
                        printf("           the current distance is %.4lf\n", sqrt(parameters.r_2));
                        printf("           the min is %.4lf and the max is %.4lf\n", sqrt(parameters.shortlimit2), sqrt(parameters.longlimit2));
                        printf("           %s:%i\n", __FILE__, __LINE__);
                    }
                    
                    JobAssign(collision_i, collision_j, &event, Coli_Event);
                }
            }
            target_j = sourCellList[target_j];
        }
    }
    
    return;
}


void BondTime(struct AtomStr *bond_i, struct AtomListStr *thisList) {
    int target_i, target_j;
    int connect[16] = {0}; //4+10 real bonds and pseudo bonds, 15 -> safe num
    double dmin[16] = {0}, dmax[16] = {0};
    struct AtomStr *bond_j;
    struct InteractionEventStr event;
    struct ParameterStr parameters;
    
    InitializeEvent(&event);
    target_i = bond_i->property->num;
    TRANSFER_VECTOR(parameters.position_i, bond_i->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i,    bond_i->dynamic->velocity);
    
    //---------------
    //transfer bond information
    int m = 0;
    struct ConstraintStr *thisBond = bond_i->property->bond;
    while (thisBond != NULL) {
        m++;
        connect[m] = thisBond->connection;
        dmin[m]    = thisBond->dmin;
        dmax[m]    = thisBond->dmax;
        thisBond   = thisBond->next;
    }
    
    //separate list between loop
    m++;
    connect[m] = 0;
    dmin[m]    = 0;
    dmax[m]    = 0;
    //---------------
    
    int n = 1;
    while (connect[n] != 0) {
        
        target_j = connect[n];
        bond_j = thisList->ptr[target_j];
        
        if (CheckDuplication(bond_i, bond_j, Bond_Event)) {
            n++;
            continue;
        }
        
        /*
         here although we are using the parameters x2, they are not
         squared. instead, they are the real length
         */
        parameters.longlimit2  = dmax[n]; //bonds can fluctuate a little bit
        parameters.shortlimit2 = dmin[n];
        
        TRANSFER_VECTOR(parameters.position_j, bond_j->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j,    bond_j->dynamic->velocity);
        CalculateParameters(bond_i, bond_j, &parameters);
        
        if (unlikely(parameters.r_2 - parameters.longlimit2 > ZERO)) {
            printf("!!ERROR!!: long-error between atoms %i(%i%s) and %i(%i%s)\n",
                   target_i, bond_i->property->sequence.aminoacidNum, bond_i->property->name,
                   target_j, bond_j->property->sequence.aminoacidNum, bond_j->property->name);
            printf("           the current distance is %.4lf, while the max is %.4lf\n", sqrt(parameters.r_2), sqrt(parameters.longlimit2));
            printf("           %s:%i\n", __FILE__, __LINE__);
            
            if (parameters.b_ij >= 0) { //try to correct the error
                event.subEventType = WellBounc;
                event.time = 1e-16 * target_i * target_j;
                JobAssign(bond_i, bond_j, &event, Bond_Event);
                n ++;
                continue;
            }
            
        } else if (unlikely(parameters.shortlimit2 - parameters.r_2 > ZERO)) {
            printf("!!ERROR!!: short-error between atoms %i(%i%s) and %i(%i%s)\n",
                   target_i, bond_i->property->sequence.aminoacidNum, bond_i->property->name,
                   target_j, bond_j->property->sequence.aminoacidNum, bond_j->property->name);
            printf("           the current distance is %.4lf, while the min is %.4lf\n", sqrt(parameters.r_2), sqrt(parameters.shortlimit2));
            printf("           %s:%i\n", __FILE__, __LINE__);
            
            if (parameters.b_ij < 0) { //try to correct the error
                event.subEventType = CoreColli;
                event.time = 1e-16 * target_i * target_j;
                JobAssign(bond_i, bond_j, &event, Bond_Event);
                n ++;
                continue;
            }
        }
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) {
            
            parameters.s_time  = -1;
            parameters.d_ij2   = parameters.shortlimit2;
            event.subEventType = CoreColli;
            
        } else {
            
            parameters.s_time  = 1;
            parameters.d_ij2   = parameters.longlimit2;
            event.subEventType = WellBounc;
            
        }
        
        event.time = CalculateTime(&parameters);

        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: bond time is less than zero between atoms %i and %i %s:%i\n", target_i, target_j, __FILE__, __LINE__);
            event.time = INFINIT;
        }
        
        JobAssign(bond_i, bond_j, &event, Bond_Event);
        n++;
    } //while (connect[n]!=0)
}


void HBTime(struct AtomStr *HB_i, struct AtomListStr *thisList) {
    int target_i, target_j;
    struct AtomStr *HB_j;
    struct InteractionEventStr event;
    struct ParameterStr parameters;
    
    InitializeEvent(&event);
    target_i = HB_i->property->num;
    TRANSFER_VECTOR(parameters.position_i, HB_i->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i,    HB_i->dynamic->velocity);
    
    if (HB_i->dynamic->HB.interactionType == HBSkip) {
        /*
         this means the previous predicted event is HB event,
         but due to the inappropriate positions of the neighbor atoms,
         the HB formation could not be executed.
         After that, in the normal condition, the target atom will keep moving
         as the original speed and direction until the next event.
         
         In this case, there would be a chance that the target atoms, say A and B,
         are very close after the previous HB event but their distance is still larger
         than potentialPairHB[][].dmax. Thus, if we do HB time prediction again,
         the next timestep for these atoms will be INFINITLY small,
         which will get the whole process stuck...
         
         Therefore, here we should skip the time prediction for one time.
         */
        HB_i->dynamic->HB.interactionType = NoEvent;
        return;
    }
    
    HB_i->dynamic->HB.interactionType = NoEvent;
    if (HB_i->dynamic->HB.bondConnection > 0) { //check if there has been already a connection between the two target beads
        
        target_j = HB_i->dynamic->HB.bondConnection;
        HB_j = thisList->ptr[target_j];

        if (CheckDuplication(HB_i, HB_j, HB_Event)) {
            return;
        }
        
        TRANSFER_VECTOR(parameters.position_j, HB_j->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j, HB_j->dynamic->velocity);
        CalculateParameters(HB_i, HB_j, &parameters);
        
        parameters.coreorShell = FindPair(HB_i, HB_j, "HB",
                                          parameters.b_ij, parameters.r_2,
                                          &parameters.shortlimit2, &parameters.longlimit2,
                                          &parameters.lowerPotential, &parameters.upperPotential,
                                          &tmpDouble, NULL, 0);
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) {
            
            parameters.s_time    = -1;
            parameters.d_ij2     = parameters.shortlimit2;
                 event.potential = parameters.lowerPotential;
            
            if (parameters.coreorShell == 0) {
                event.subEventType = CoreColli;
            } else {
                event.subEventType = TBDCaptr;
            }
            
        } else {
            parameters.s_time    = 1;
            parameters.d_ij2     = parameters.longlimit2;
                 event.potential = parameters.upperPotential;
            
            if (strncmp(Methodtype, "Ding", 1) == 0) {
                if (parameters.coreorShell == 1) { //most outer shell, may break
                    event.subEventType = TBDHBBrik;
                } else {
                    event.subEventType = TBDEscap;
                }
            }
        }
        
        event.time = CalculateTime(&parameters);
        
        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: HB time is less than zero between atoms %i and %i! %s:%i\n", target_i, target_j, __FILE__, __LINE__);
            event.time = INFINIT;
        }
        
        JobAssign(HB_i, HB_j, &event, HB_Event);
        
    } else if (target_i != 0 &&
               (HB_i->dynamic->HB.role == 'A' ||
                HB_i->dynamic->HB.role == 'D' ||
                HB_i->dynamic->HB.role == 'S') &&
               HB_i->dynamic->HB.bondConnection == 0) {
        
        parameters.calibratedcutoffr = cutoffr * EXTEND_RATIO;
        parameters.calibratedcutoffr *= parameters.calibratedcutoffr;
        
        int cell_neighbor[4], cellIndex;
        int *sourCellList = celllist;
        
        for (int i = 0; i < 27; ++i) {
            //scan the neighborhood 27 subcells, include the target subcell itself
            cell_neighbor[3] = HB_i->dynamic->cellIndex[3] + thisList->z[i];
            cell_neighbor[2] = HB_i->dynamic->cellIndex[2] + thisList->y[i];
            cell_neighbor[1] = HB_i->dynamic->cellIndex[1] + thisList->x[i];
            
            for (int j = 1; j <= 3; j++) {
                if (unlikely(cell_neighbor[j] < 0)) {
                    cell_neighbor[j] = cellnum[j] - 1;
                } else if (unlikely(cell_neighbor[j] >= cellnum[j])) {
                    cell_neighbor[j] = 0;
                }
            }
            
            cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2]
            + cell_neighbor[2] * cellnum[1]
            + cell_neighbor[1] + 1;
        
            target_j = sourCellList[atomnum + cellIndex];
            while (target_j) {
                HB_j = thisList->ptr[target_j];
                
                if (target_i != target_j &&
                    (/*normal A and D*/
                     HBPAIR_AD(HB_i, HB_j) ||
                     /*S and D, here D has to be in backbone*/
                     HBPAIR_SD(HB_i, HB_j) ||
                     /*S and S*/
                     HBPAIR_SS(HB_i, HB_j)
                     ) &&
                    HB_j->dynamic->HB.bondConnection == 0 &&
                    ((HB_i->property->sequence.proteinNum == HB_j->property->sequence.proteinNum &&
                      ABSVALUE((HB_i->property->sequence.aminoacidNum - HB_j->property->sequence.aminoacidNum)) >= 4) ||
                     (HB_i->property->sequence.proteinNum != HB_j->property->sequence.proteinNum)))
                {
                    
                    if (CheckDuplication(HB_i, HB_j, HB_Event)) {
                        target_j = sourCellList[target_j];
                        continue;
                    }
                    
                    TRANSFER_VECTOR(parameters.position_j, HB_j->dynamic->coordinate);
                    TRANSFER_VECTOR(parameters.speed_j,    HB_j->dynamic->velocity);
                    CalculateParameters(HB_i, HB_j, &parameters);
                    
                    if (parameters.r_2 <= parameters.calibratedcutoffr) { //if the distance between two atoms <= scan radius
                        
                        if (CheckPairMatch(HB_i, HB_j, thisList) < 0) { //the pair is not in the reaction library list
                            target_j = sourCellList[target_j];
                            continue;
                        }
                        
                        parameters.coreorShell = FindPair(HB_i, HB_j, "HB",
                                                          parameters.b_ij, parameters.r_2,
                                                          &parameters.shortlimit2, &parameters.longlimit2,
                                                          &parameters.lowerPotential, &parameters.upperPotential,
                                                          &tmpDouble, NULL, 0);
                        
                        if (parameters.b_ij < 0 && parameters.coreorShell == 1 &&
                            CalculateDisc(&parameters) > 0) { //approaching
                            
                            parameters.s_time       = -1;
                            parameters.d_ij2        = parameters.shortlimit2;
                            event.potential         = parameters.lowerPotential;
                            event.subEventType      = TBDHBForm; //may form HB
                            
                        } else {
                            target_j = sourCellList[target_j];
                            continue;
                        }
                        
                        event.time = CalculateTime(&parameters);
                        
                        if (unlikely(event.time < 0)) {
                            printf("!!ERROR!!: HB time is less than zero bewteen atoms %i and %i! %s:%i\n", target_i, target_j, __FILE__, __LINE__);
                            event.time = INFINIT;
                        }
                        
                        JobAssign(HB_i, HB_j, &event, HB_Event);
                    } //if (absvalue_vector(r_HB_ij[0])<=calibratedcutoffr)
                }
                target_j = sourCellList[target_j];
            }
        }
    }
    
    return;
}


void HBNeighborTime(struct AtomStr *neighbor_i, struct AtomListStr *thisList) {
    int flag = 0;
    int target_i, target_j;
    struct AtomStr *neighbor_j;
    struct AtomStr *virtualAtom[3];
    struct InteractionEventStr event;
    struct ParameterStr parameters;
    
    InitializeEvent(&event);
    target_i = neighbor_i->property->num;
    TRANSFER_VECTOR(parameters.position_i, neighbor_i->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i,    neighbor_i->dynamic->velocity);
    
    for (int n = 1; n <= neighbor_i->dynamic->HBNeighbor.neighborStatus; n++) {
        
        target_j = neighbor_i->dynamic->HBNeighbor.neighborPartner[n];
        neighbor_j   = thisList->ptr[target_j];

        if (flag == 0) {
            neighbor_i->dynamic->HBNeighbor.partnerNum = n; //decide which host will be used to do calculation
        }
        
        if (CheckDuplication(neighbor_i, neighbor_j, HBNe_Event)) {
            neighbor_i->dynamic->HBNeighbor.partnerNum = n;
            flag = 1;
            continue;
        }
        
        TRANSFER_VECTOR(parameters.position_j, neighbor_j->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j,    neighbor_j->dynamic->velocity);
        CalculateParameters(neighbor_i, neighbor_j, &parameters);
        
        //==============================
        //check which atom is HOST and which atom is NEIGHBOR
        if ((neighbor_i->dynamic->HB.role == 'A' ||
             neighbor_i->dynamic->HB.role == 'D' ||
             neighbor_i->dynamic->HB.role == 'S') &&
            (neighbor_i->dynamic->HB.bondConnection && thisList->ptr[neighbor_i->dynamic->HB.bondConnection]->dynamic->HB.neighbor == target_j)) { //both could be either 'A' or 'D'
            virtualAtom[0] = neighbor_i; //host
            virtualAtom[1] = neighbor_j; //neighbor
            virtualAtom[2] = thisList->ptr[neighbor_i->dynamic->HB.bondConnection];
        } else {
            virtualAtom[0] = neighbor_j;
            virtualAtom[1] = neighbor_i;
            virtualAtom[2] = thisList->ptr[neighbor_j->dynamic->HB.bondConnection];
        }
        //==============================
        
        parameters.coreorShell = FindPair(virtualAtom[1], virtualAtom[0], "neighbor",
                                          parameters.b_ij, parameters.r_2,
                                          &parameters.shortlimit2, &parameters.longlimit2,
                                          &parameters.lowerPotential, &parameters.upperPotential,
                                          &tmpDouble, virtualAtom[2], 0);
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) {
            
            parameters.s_time    = -1;
            parameters.d_ij2     = parameters.shortlimit2;
                 event.potential = parameters.lowerPotential;
            
            if (parameters.coreorShell == 0) { //inner well contact
                event.subEventType = CoreColli;
            } else { //outer well contact
                event.subEventType = TBDCaptr;
            }
            
        } else {
            
            parameters.s_time    = 1;
            parameters.d_ij2     = parameters.longlimit2;
                 event.potential = parameters.upperPotential;
            
            if (strncmp(Methodtype, "Ding", 1) == 0) {
                if (parameters.coreorShell == 1) { //most outer well contact
                    event.subEventType = WellBounc;
                } else { //outer well contact (pull)
                    event.subEventType = TBDEscap;
                }
            }

        }
        
        event.time = CalculateTime(&parameters);
        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: HB neighbor time is less than zero between atoms %i and %i! %s:%i\n", target_i, target_j, __FILE__, __LINE__);
            event.time = INFINIT;
        }
        
        if (JobAssign(neighbor_i, neighbor_j, &event, HBNe_Event)) {
            neighbor_i->dynamic->HBNeighbor.partnerNum = n;
            flag = 1;
        }
    }
}


void ThermostatTime(struct AtomStr *target) {
    double temp;
    struct InteractionEventStr event;
    
    if (strcmp(thermostatType, "no") == 0) {
        return;
    }
    
    InitializeEvent(&event);
    temp = Maxwell_Boltzmann_Distribution(0, 1);
    event.time = thermoF * ABSVALUE(temp);
    
    JobAssign(target, NULL, &event, Ther_Event);
}


void PBCandCrossCellTime (struct AtomStr *target) {
    double position[4], speed[4];
    double boundary[2][4]; //0: lower; 1: higher
    double shortertime = INFINIT;
    struct InteractionEventStr event;
    
    InitializeEvent(&event);
    TRANSFER_VECTOR(position, target->dynamic->coordinate);
    TRANSFER_VECTOR(speed, target->dynamic->velocity);
    
    for (int n = 1; n <= 3; n++) {
        boundary[0][n] = target->dynamic->cellIndex[n] * cellsize[n];
        boundary[1][n] = (target->dynamic->cellIndex[n] + 1) * cellsize[n];
    }
    
    for (int n = 1; n <= 3; n++) {
        if (speed[n] > 0) {
            shortertime = (boundary[1][n] - position[n] + ZERO) / speed[n];
            
            if (event.time > shortertime && shortertime > 0) {
                event.time = shortertime;
            }
        } else if (speed[n] < 0) {
            shortertime = (boundary[0][n] - position[n] - ZERO) / speed[n];
            
            if (event.time > shortertime && shortertime > 0) {
                event.time = shortertime;
            }
        } else {
            printf("PCC speed is equal to zero?! %s:%i\n", __FILE__, __LINE__);
        }
    }
    
    JobAssign(target, NULL, &event, PCC_Event);
}


void WallTime(struct AtomStr *target) {
    int atomNum = target->property->num;
    double boxSize[4];
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
ReCal:
    InitializeEvent(&event);
    TRANSFER_VECTOR(boxSize, boxDimension);
    TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i, target->dynamic->velocity);
    
    if (strncmp(wallExist, "smooth", 2) == 0) {
        
        double distance2 = 0;
        struct AtomStr *thisWall = &wall[0];
        TRANSFER_VECTOR(parameters.position_j, thisWall->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j, thisWall->dynamic->velocity);
        
        if (strncmp(wallType, "parallel", 1) == 0) { //parallel plates of wall
            parameters.position_i[1] = parameters.position_i[3] = 0;
            parameters.speed_i[1]    = parameters.speed_i[3]    = 0;
        } else if (strncmp(wallType, "cylinder", 1) == 0) {
            parameters.position_i[1] = 0;
            parameters.speed_i[1]    = 0;
        }
        
        DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
        DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
        parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
        parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
        parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
        distance2 = boxSize[2] * boxSize[2] * 0.25 + parameters.r_2 - boxSize[2] * sqrt(parameters.r_2);
        
        parameters.coreorShell = FindPair(target, thisWall, "collision",
                                          -1 * parameters.b_ij, distance2,
                                          &parameters.shortlimit2, &parameters.longlimit2,
                                          &parameters.lowerPotential, &parameters.upperPotential,
                                          &tmpDouble, NULL, 0);
        
        if (parameters.lowerPotential == 0 &&
            parameters.upperPotential == 0) {
            
            if (strcmp(wallDyn.mark, "no")) {
                wallDyn.touch ++;
                
                if (strncmp(wallType, "parallel", 1) == 0) {
                    boxDimension[2] += wallDyn.step;
                } else if (strncmp(wallType, "cylinder", 1) == 0) {
                    boxDimension[2] += wallDyn.step;
                    boxDimension[3] += wallDyn.step;
                } else if (strncmp(wallType, "sphere", 1) == 0) {
                    boxDimension[1] += wallDyn.step;
                    boxDimension[2] += wallDyn.step;
                    boxDimension[3] += wallDyn.step;
                }
                
                if (boxDimension[2] > wallDyn.origBoxDim[2]) {
                    printf("!!ERROR!!: box size cannot be expanded beyond the original! %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                
                goto ReCal;
            }

            printf("!!ERROR!!: atom #%i is outside of the wall! %s:%i\n", atomNum, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) { //approaching to the center
            if (parameters.coreorShell == 1) {
                return;
            } else {
                event.subEventType = TBDCaptr;
            }
            
            parameters.s_time    = -1;
            parameters.d_ij2     = boxSize[2] * boxSize[2] * 0.25 + parameters.longlimit2 - boxSize[2] * sqrt(parameters.longlimit2);
            event.potential = parameters.upperPotential;
            
        } else { //approaching to the wall edge
            if (parameters.coreorShell == 0) {
                event.subEventType = WellBounc;
            } else {
                event.subEventType = TBDEscap;
            }
            
            parameters.s_time    = 1;
            parameters.d_ij2     = boxSize[2] * boxSize[2] * 0.25 + parameters.shortlimit2 - boxSize[2] * sqrt(parameters.shortlimit2);
            event.potential = parameters.lowerPotential;
        }
        
        event.time = CalculateTime(&parameters);
        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: wall time is less than zero, atom #%i! %s:%i\n", atomNum, __FILE__, __LINE__);
            event.time = INFINIT;
        }
        
        JobAssign(target, NULL, &event, Wall_Event);
    }
    
    return;
}


void ObstTime(struct AtomStr *target) {
    int atomNum = target->property->num;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    for (int n = 0; n < obstObj.num; n ++) {
        InitializeEvent(&event);

        struct AtomStr *thisWall = &obstObj.obst[n];
        for (int dim = 1; dim <= 3; dim ++) {
            if (obstObj.position[n][dim] < 0) {
                continue;
            }
            
            double distance2 = 0;
            TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
            TRANSFER_VECTOR(parameters.speed_i, target->dynamic->velocity);
            
            TRANSFER_VECTOR(parameters.position_j, thisWall->dynamic->coordinate);
            TRANSFER_VECTOR(parameters.speed_j, thisWall->dynamic->velocity);
            
            for (int i = 1; i <= 3; i ++) {
                if (i != dim) {
                    parameters.position_i[i] = 0;
                    parameters.position_j[i] = 0;
                    
                    parameters.speed_i[i] = 0;
                    parameters.speed_j[i] = 0;
                }
            }
            
            DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
            DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
            parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
            parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
            parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
            distance2 = parameters.r_2;
            
            parameters.coreorShell = FindPair(target, thisWall, "collision",
                                              parameters.b_ij, distance2,
                                              &parameters.shortlimit2, &parameters.longlimit2,
                                              &parameters.lowerPotential, &parameters.upperPotential,
                                              &tmpDouble, NULL, 0);
            
            if (parameters.lowerPotential == 0 &&
                parameters.upperPotential == 0) {
                //already in the holes, skip
                continue;
            }
            
            if (parameters.b_ij < 0) { //approaching
                if (parameters.coreorShell == 0) {
                    event.subEventType = CoreColli;
                } else {
                    event.subEventType = TBDCaptr;
                }
                
                parameters.s_time    = -1;
                parameters.d_ij2     = parameters.shortlimit2;
                event.potential = parameters.lowerPotential;
                
            } else { //leaving
                if (parameters.coreorShell == 1) {
                    continue;
                } else {
                    event.subEventType = TBDEscap;
                }
                
                parameters.s_time    = 1;
                parameters.d_ij2     = parameters.longlimit2;
                event.potential = parameters.upperPotential;
            }
            
            event.time = CalculateTime(&parameters);
            if (unlikely(event.time < 0)) {
                printf("!!ERROR!!: obstacle time is less than zero, atom #%i! %s:%i\n", atomNum, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            //at the next collision, the target atom would have been in the holes,
            //then the current prediction should be invalid
            if (!(obstObj.hole[n] && EnterTunnel(event.time, target)) &&
                JobAssign(target, NULL, &event, Obst_Event)) {
                target->dynamic->event.partner = -1 * dim - 10 * n;
            }
        }
    }
    
    return;
}


void TunnelTime(struct AtomStr *target) {
    int atomNum = target->property->num;
    double radius2;
    double distance2 = 0;
    double distanceShift2;
    double extraRange;
    double futurePos;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    struct AtomStr *thisWall = &tunlObj.tunnel;
    struct ConstraintStr *thisConstr = RightPair(target->property->typeofAtom, thisWall->property->typeofAtom, 0);
    
    distanceShift2 = thisConstr->dmin;
    extraRange = sqrt(FindPotWellWidth(thisWall, target));
    
    TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i,    target->dynamic->velocity);
    parameters.position_i[0] = parameters.position_i[1];
    parameters.speed_i[0]    = parameters.speed_i[1];

    for (int n = 0; n < tunlObj.num; n ++) {
        radius2 = tunlObj.diameter[n] / 2;
        radius2 *= radius2;
        
        InitializeEvent(&event);
        TRANSFER_VECTOR(parameters.position_j, tunlObj.position[n]);
        TRANSFER_VECTOR(parameters.speed_j,    thisWall->dynamic->velocity);
        
        parameters.position_i[1] = parameters.position_j[1]; //end position
        parameters.speed_i[1]    = 0;
        
        DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
        DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
        parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
        parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
        parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
        
        if (parameters.r_2 > radius2 + distanceShift2 - 2 * sqrt(radius2 * distanceShift2) + ZERO) {
            continue; //currently outside of the tunnel
        }
        distance2 = radius2 + parameters.r_2 - 2 * sqrt(parameters.r_2 * radius2);
        
        parameters.coreorShell = FindPair(target, thisWall, "collision",
                                          -1 * parameters.b_ij, distance2,
                                          &parameters.shortlimit2, &parameters.longlimit2,
                                          &parameters.lowerPotential, &parameters.upperPotential,
                                          &tmpDouble, NULL, 0);
        
        if (unlikely(parameters.lowerPotential == 0 &&
                     parameters.upperPotential == 0)) {
            printf("!!ERROR!!: atom #%i is outside of the wall! %s:%i\n", atomNum, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) { //approaching to the center
            if (parameters.coreorShell == 1) {
                continue;
            } else {
                event.subEventType = TBDCaptr;
            }
            
            parameters.s_time    = -1;
            parameters.d_ij2     = radius2 + parameters.longlimit2 - 2 * sqrt(parameters.longlimit2 * radius2);
            event.potential = parameters.upperPotential;
            
        } else { //approaching to the wall edge
            if (parameters.coreorShell == 0) {
                event.subEventType = WellBounc;
            } else {
                event.subEventType = TBDEscap;
            }
            
            parameters.s_time    = 1;
            parameters.d_ij2     = radius2 + parameters.shortlimit2 - 2 * sqrt(parameters.shortlimit2 * radius2);
            event.potential = parameters.lowerPotential;
        }
        
        event.time = CalculateTime(&parameters);
        
        //at the next event, the target atom may have left the tunnel,
        //then the current prediction should be invalid
        futurePos = parameters.position_i[0] + parameters.speed_i[0] * event.time;
        if (futurePos < tunlObj.position[n][0] - extraRange ||
            futurePos > tunlObj.position[n][1] + extraRange) {
            continue;
        }
        
        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: wall time is less than zero, atom #%i! %s:%i\n", atomNum, __FILE__, __LINE__);
            event.time = INFINIT;
        }
        
        if (JobAssign(target, NULL, &event, Tunl_Event)) {
            target->dynamic->event.partner = -10000 * n;
        }
    }
    
    return;
}


void ChargeTime(struct AtomStr *target) {
    int flag = 0;
    int gapNum;
    int totalCNum = (flow.charge.PBCMark) ? flow.charge.num : flow.charge.num / 2;
    double charge;
    double distance2, gap2;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    if (!target->property->charge) {
        return;
    }
    charge = target->property->charge;
    
    for (int chargeNum = 0; chargeNum < totalCNum; chargeNum ++) {
        if (target->dynamic->coordinate[1] > flow.charge.position[chargeNum][1] - ZERO) {
            continue;
        }
        
        //if the atom have interacted with the charge in the system, then it will not interact
        //with the pseudo-charge outside of the pbc box any more
        if (flow.charge.position[chargeNum][1] < boxDimension[1]) {
            flag = 1;
        } else if (!(flag == 0 && flow.charge.position[chargeNum][1] > boxDimension[1])) {
            break;
        }
        
        InitializeEvent(&event);
        
        TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_i, target->dynamic->velocity);
        
        TRANSFER_VECTOR(parameters.position_j, flow.charge.position[chargeNum]);
        TRANSFER_VECTOR(parameters.speed_j, flow.charge.velocity[chargeNum]);
        
        DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
        DOT_MINUS(parameters.speed_i, parameters.speed_j, parameters.v_ij);
        parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
        parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
        parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
        
        gapNum = sqrt(parameters.r_2) / flow.charge.gap[chargeNum];
        gap2 = flow.charge.gap[chargeNum] * flow.charge.gap[chargeNum];
        distance2 = gapNum * gapNum * gap2;
        
        if (parameters.b_ij < 0 &&
            parameters.r_2 < distance2 + ZERO) {
            gapNum --;
        } else if (parameters.b_ij > 0 &&
                   parameters.r_2 > (gapNum + 1) * (gapNum + 1) * gap2 - ZERO) {
            gapNum ++;
        }
        
        parameters.upperPotential = flow.charge.potGrad[chargeNum] * charge;
        parameters.longlimit2  = (gapNum + 1) * (gapNum + 1) * gap2;
        
        parameters.lowerPotential = -1 * flow.charge.potGrad[chargeNum] * charge;
        parameters.shortlimit2 = gapNum * gapNum * gap2;
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) {
            parameters.s_time    = -1;
            parameters.d_ij2     = parameters.shortlimit2;
            
            event.potential      = parameters.lowerPotential;
            event.subEventType   = TBDCaptr;
        } else {
            parameters.s_time    = 1;
            parameters.d_ij2     = parameters.longlimit2;
            
            event.potential      = parameters.upperPotential;
            event.subEventType   = TBDEscap;
        }
        
        event.time = CalculateTime(&parameters);
        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: charge interaction time is less than zero between atom %i and charge %i! %s:%i\n", target->property->num, chargeNum, __FILE__, __LINE__); //for debug and error check
            event.time = INFINIT;
        }
        
        if (JobAssign(target, NULL, &event, Chrg_Event))
            target->dynamic->event.partner = -100 * chargeNum;
    }
    
    return;
}


void SphObstTime(struct AtomStr *target) {
    int atomNum = target->property->num;
    double distance2 = 0;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    for (int n = 0; n < SphObstObj.num; n ++) {
        InitializeEvent(&event);
        
        TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_i,    target->dynamic->velocity);
        
        TRANSFER_VECTOR(parameters.position_j, SphObstObj.position[n]);
        parameters.speed_j[1] = 0;
        parameters.speed_j[2] = 0;
        parameters.speed_j[3] = 0;
        
        DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
        DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
        parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
        parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
        parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
        distance2 = parameters.r_2;
        
        parameters.shortlimit2 = SphObstObj.radius[n] * SphObstObj.radius[n];
        parameters.lowerPotential = 0;
        
        parameters.longlimit2 = INFINIT;
        parameters.upperPotential = INFINIT;
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) { //approaching
            event.subEventType   = CoreColli;
            parameters.s_time    = -1;
            parameters.d_ij2     = parameters.shortlimit2;
            event.potential      = parameters.lowerPotential;
        } else {
            continue;
        }
        
        event.time = CalculateTime(&parameters);
        if (unlikely(event.time < 0)) {
            printf("!!ERROR!!: spherical obstacle time is less than zero, atom #%i! %s:%i\n", atomNum, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        if (JobAssign(target, NULL, &event, SphO_Event)) {
            target->dynamic->event.partner = -1000 * n;
        }
    }
    
    return;
}


void InitializeEvent(struct InteractionEventStr *event) {
    event->time         = INFINIT;
    event->potential    = 0;
    event->eventType    = Invd_Event;
    event->subEventType = Invalid_Inter;
}

void CalculateParameters(struct AtomStr *atom_i, struct AtomStr *atom_j, struct ParameterStr *parameters) {
    PBCShift(atom_i, atom_j, parameters->positionshift);
    DOT_PLUS(parameters->position_j, parameters->positionshift, parameters->r_temp);
    DOT_MINUS(parameters->position_i, parameters->r_temp,        parameters->r_ij);
    DOT_MINUS(parameters->speed_i,    parameters->speed_j,       parameters->v_ij);
    
    parameters->r_2  = DOT_PROD(parameters->r_ij, parameters->r_ij);
    parameters->b_ij = DOT_PROD(parameters->r_ij, parameters->v_ij); //b_ij is before collision here
    parameters->v_2  = DOT_PROD(parameters->v_ij, parameters->v_ij);
}

double CalculateDisc(struct ParameterStr *parameters) {
    return parameters->b_ij * parameters->b_ij - parameters->v_2 * (parameters->r_2 - parameters->shortlimit2);
}

double CalculateTime(struct ParameterStr *parameters) {
    return (-1 * parameters->b_ij + parameters->s_time *
            sqrt(parameters->b_ij * parameters->b_ij - parameters->v_2 * (parameters->r_2 - parameters->d_ij2))) / parameters->v_2;
}

int CheckDuplication(struct AtomStr *atom_i, struct AtomStr *atom_j, enum EEventType type) {
    if (   atom_j->dynamic->event.partner        == atom_i->property->num
        && atom_j->dynamic->event.eventType      == type
        && atom_j->dynamic->event.partnerCounter == atom_i->dynamic->event.counter
        && atom_i->dynamic->event.time           >  atom_j->dynamic->event.time) {
        
        atom_i->dynamic->event.partner        = atom_j->property->num;
        atom_i->dynamic->event.potential      = atom_j->dynamic->event.potential;
        atom_i->dynamic->event.subEventType   = atom_j->dynamic->event.subEventType;
        atom_i->dynamic->event.time           = atom_j->dynamic->event.time;
        atom_i->dynamic->event.eventType      = type;
        atom_i->dynamic->event.partnerCounter = atom_j->dynamic->event.counter;
        return TRUE;
    }
    
    return FALSE;
}

int CheckPairMatch(struct AtomStr *atom_i, struct AtomStr *atom_j, struct AtomListStr *thisList) {
	int type = HBModel(atom_i, atom_j);
	if (type < 0) return -1;

	//atom_i and atom_j must be their original types
	int atom_i_type = atom_i->property->typeofAtom;
	int atom_j_type = atom_j->property->typeofAtom;
	if (potentialPairHB[type][atom_i_type][atom_j_type].step == NULL) {
		return -1;
	}

	int neighbor_i_type = AtomTypeChange(thisList->ptr[atom_i->dynamic->HB.neighbor]->property->typeofAtom, 0);
	int neighbor_j_type = AtomTypeChange(thisList->ptr[atom_j->dynamic->HB.neighbor]->property->typeofAtom, 0);
	if (potentialPairHB[type][atom_i_type][neighbor_j_type].step == NULL ||
		potentialPairHB[type][atom_j_type][neighbor_i_type].step == NULL) {
		return -1;
	}
	
	return 0;
}

int JobAssign(struct AtomStr *atom_i, struct AtomStr *atom_j, struct InteractionEventStr *event, enum EEventType type) {
    
    if (event->time > 0 && atom_i->dynamic->event.time > event->time) {
        atom_i->dynamic->event.potential    = event->potential;
        atom_i->dynamic->event.subEventType = event->subEventType;
        atom_i->dynamic->event.time         = event->time;
        atom_i->dynamic->event.eventType    = type;
        
        if (atom_j) {
            atom_i->dynamic->event.partner        = atom_j->property->num;
            atom_i->dynamic->event.partnerCounter = atom_j->dynamic->event.counter;
        } else {
            atom_i->dynamic->event.partner = INVALID;
        }
        
        return TRUE;
    }
    return FALSE;
}

void CalculateVCOM(int startNum, int endNum, struct AtomStr *atomList, double *vCOM) {
    double mass = 0;
    
    vCOM[1] = vCOM[2] = vCOM[3] = 0;

    for (int n = startNum; n <= endNum; n ++) {
        mass += atomList[n].property->mass;
        for (int j = 0; j <= 3; j ++) {
            vCOM[j] += atomList[n].dynamic->velocity[j] * atomList[n].property->mass;
        }
    }
    
    for (int j = 1; j <= 3; j ++) {
        vCOM[j] /= mass;
    }
    
    return;
}
