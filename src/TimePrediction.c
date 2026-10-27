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
(atom1->dynamic->HB.role == 'S' && atom2->property->typeofAtom[1] == 2) || \
(atom2->dynamic->HB.role == 'S' && atom1->property->typeofAtom[1] == 2)

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

typedef void (*InteractionCallback)(struct AtomStr *atom_i, struct AtomStr *atom_j, struct ParameterStr *parameters, struct AtomListStr *thisList);
// cheap pair test run before the neighbor is updated and measured; NULL accepts every pair
typedef int (*PairFilter)(struct AtomStr *atom_i, struct AtomStr *atom_j);

static void ScanNeighborhood(struct AtomStr *atom_i, struct AtomListStr *thisList, PairFilter filter, InteractionCallback callback) {
    int target_i, target_j;
    int cell_neighbor[4], cellIndex;
    int *sourCellList = celllist;
    struct AtomStr *atom_j;
    struct ParameterStr parameters;
    
    parameters.calibratedcutoffr = cutoffr * EXTEND_RATIO;
    parameters.calibratedcutoffr *= parameters.calibratedcutoffr;
    
    target_i = atom_i->property->num;
    TRANSFER_VECTOR(parameters.position_i, atom_i->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i,    atom_i->dynamic->velocity);
    
    for (int i = 0; i < thisList->stencilNum; ++i) {
        cell_neighbor[3] = atom_i->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = atom_i->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = atom_i->dynamic->cellIndex[1] + thisList->x[i];
        
        for (int j = 1; j <= 3; j++) {
            if (unlikely(cell_neighbor[j] < 0)) cell_neighbor[j] += cellnum[j];
            else if (unlikely(cell_neighbor[j] >= cellnum[j])) cell_neighbor[j] -= cellnum[j];
        }
        
        cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2] + cell_neighbor[2] * cellnum[1] + cell_neighbor[1] + 1;
        
        target_j = sourCellList[atomnum + cellIndex];
        while (target_j) {
            atom_j = thisList->ptr[target_j];
            if (target_i != target_j && (filter == NULL || filter(atom_i, atom_j))) {
                if (!FindStrInList(atom_j->property->extraProperty, "_HBRMV", 0) &&
                    !FindStrInList(atom_j->property->extraProperty, "_COLSKIP", 0)) {
                    
                    UpdateAtom(atom_j, currenttime); // Update neighbor atom
                    TRANSFER_VECTOR(parameters.position_j, atom_j->dynamic->coordinate);
                    TRANSFER_VECTOR(parameters.speed_j,    atom_j->dynamic->velocity);
                    CalculateParameters(atom_i, atom_j, &parameters);
                    
                    if (parameters.r_2 <= parameters.calibratedcutoffr) {
                        callback(atom_i, atom_j, &parameters, thisList);
                    }
                }
            }
            target_j = sourCellList[target_j];
        }
    }
}

static void CollisionCallback(struct AtomStr *atom_i, struct AtomStr *atom_j, struct ParameterStr *parameters, struct AtomListStr *thisList) {
    struct InteractionEventStr event;
    InitializeEvent(&event);
    int target_i = atom_i->property->num;
    int target_j = atom_j->property->num;

    if (CheckDuplication(atom_i, atom_j, Coli_Event)) return;
    
    if (connectionMap[target_i][target_j] & BOND_CONNECT ||
        connectionMap[target_i][target_j] & HB_CONNECT   ||
        FindElemInList(atom_i->dynamic->HBNeighbor.neighborPartner, target_j, -1, -1, 0)) return;

    if (connectionMap[target_i][target_j] & CONSTRAINT_CONNECT) {
        parameters->coreorShell = FindPair(atom_i, atom_j, InterConstr, parameters->b_ij, parameters->r_2, &parameters->shortlimit2, &parameters->longlimit2, &parameters->lowerPotential, &parameters->upperPotential, &tmpDouble, NULL, 0);
    } else {
        parameters->coreorShell = FindPair(atom_i, atom_j, InterColli, parameters->b_ij, parameters->r_2, &parameters->shortlimit2, &parameters->longlimit2, &parameters->lowerPotential, &parameters->upperPotential, &tmpDouble, NULL, 0);
        if (parameters->lowerPotential == 0 && parameters->upperPotential == 0) return;
    }

    if (parameters->b_ij < 0 && CalculateDisc(parameters) > 0) { // approaching and will reach the inner boundary
        parameters->s_time = -1;
        parameters->d_ij2  = parameters->shortlimit2;
        event.potential    = parameters->lowerPotential;
        event.subEventType = (parameters->coreorShell == 0) ? CoreColli : TBDCaptr;
    } else { // moving apart, or approaching but turning back before the inner boundary
        if (connectionMap[target_i][target_j] & CONSTRAINT_CONNECT) {
            event.subEventType = (parameters->coreorShell == 1) ? WellBounc : TBDEscap;
        } else {
            if (parameters->coreorShell == 1) return;
            event.subEventType = TBDEscap;
        }
        parameters->s_time = 1;
        parameters->d_ij2  = parameters->longlimit2;
        event.potential    = parameters->upperPotential;
    }

    event.time = CalculateTime(parameters);
    if (unlikely(event.time < 0)) event.time = INFINIT;
    JobAssign(atom_i, atom_j, &event, Coli_Event);
}

void CollisionTime(struct AtomStr *collision_i, struct AtomListStr *thisList) {
    UpdateAtom(collision_i, currenttime);
    if (FindStrInList(collision_i->property->extraProperty, "_HBRMV", 0) ||
        FindStrInList(collision_i->property->extraProperty, "_COLSKIP", 0)) return;
    ScanNeighborhood(collision_i, thisList, NULL, CollisionCallback);
}


void BondTime(struct AtomStr *bond_i, struct AtomListStr *thisList) {
    int target_i, target_j;
    int connect[16] = {0}; //4+10 real bonds and pseudo bonds, 15 -> safe num
    double dmin[16] = {0}, dmax[16] = {0};
    struct AtomStr *bond_j;
    struct InteractionEventStr event;
    struct ParameterStr parameters;
    
    UpdateAtom(bond_i, currenttime);
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
        
        UpdateAtom(bond_j, currenttime);
        TRANSFER_VECTOR(parameters.position_j, bond_j->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j,    bond_j->dynamic->velocity);
        CalculateParameters(bond_i, bond_j, &parameters);
        
        if (unlikely(parameters.r_2 - parameters.longlimit2 > ZERO)) {
            if (parameters.b_ij >= 0) { //try to correct the error
                printf("!!ERROR!!: long-error between atoms %i(%i%s) and %i(%i%s)\n",
                       target_i, bond_i->property->sequence.resNum, bond_i->property->name,
                       target_j, bond_j->property->sequence.resNum, bond_j->property->name);
                printf("           the current distance is %.4lf, while the max is %.4lf\n", sqrt(parameters.r_2), sqrt(parameters.longlimit2));
                printf("           %s:%i\n", __FILE__, __LINE__);
                event.subEventType = WellBounc;
                event.time = 1e-16 * target_i * target_j;
                JobAssign(bond_i, bond_j, &event, Bond_Event);
                n ++;
                continue;
            }

        } else if (unlikely(parameters.shortlimit2 - parameters.r_2 > ZERO)) {
            if (parameters.b_ij < 0) { //try to correct the error
                printf("!!ERROR!!: short-error between atoms %i(%i%s) and %i(%i%s)\n",
                       target_i, bond_i->property->sequence.resNum, bond_i->property->name,
                       target_j, bond_j->property->sequence.resNum, bond_j->property->name);
                printf("           the current distance is %.4lf, while the min is %.4lf\n", sqrt(parameters.r_2), sqrt(parameters.shortlimit2));
                printf("           %s:%i\n", __FILE__, __LINE__);
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

static int HBPairFilter(struct AtomStr *HB_i, struct AtomStr *HB_j) {
    if (!(HBPAIR_AD(HB_i, HB_j) || HBPAIR_SD(HB_i, HB_j) || HBPAIR_SS(HB_i, HB_j))) return FALSE;
    if (!HBHasFreeSlot(HB_j) || HBIsPartner(HB_i, HB_j->property->num)) return FALSE;
    if (HB_i->property->sequence.molNum == HB_j->property->sequence.molNum &&
        ABSVALUE((HB_i->property->sequence.resNum - HB_j->property->sequence.resNum)) < 4) return FALSE;
    return TRUE;
}

static void HBCallback(struct AtomStr *HB_i, struct AtomStr *HB_j, struct ParameterStr *parameters, struct AtomListStr *thisList) {
    struct InteractionEventStr event;
    InitializeEvent(&event);
    int target_i = HB_i->property->num;
    int target_j = HB_j->property->num;

#ifdef HB_TRACE
    if (currenttime > 340.0) {
        printf("[HB_TRACE] t=%.4f Checking H#%d(res%d) O#%d(res%d) r=%.3f\n",
               currenttime, target_i, HB_i->property->sequence.resNum,
               target_j, HB_j->property->sequence.resNum, sqrt(parameters->r_2));
        fflush(stdout);
    }
#endif

    if (!(HBPAIR_AD(HB_i, HB_j) || HBPAIR_SD(HB_i, HB_j) || HBPAIR_SS(HB_i, HB_j))) {
#ifdef HB_TRACE
        if (currenttime > 340.0) printf("[HB_TRACE] t=%.4f  -> FAIL: HBPAIR role check (roles %c %c)\n", currenttime, HB_i->dynamic->HB.role, HB_j->dynamic->HB.role);
#endif
        return;
    }
    if (!HBHasFreeSlot(HB_j) || HBIsPartner(HB_i, target_j)) {
#ifdef HB_TRACE
        if (currenttime > 340.0) printf("[HB_TRACE] t=%.4f  -> FAIL: O#%d already bonded to H#%d\n", currenttime, target_j, HB_j->dynamic->HB.bondConnection[1]);
#endif
        return;
    }
    if (!((HB_i->property->sequence.molNum == HB_j->property->sequence.molNum && ABSVALUE((HB_i->property->sequence.resNum - HB_j->property->sequence.resNum)) >= 4) || (HB_i->property->sequence.molNum != HB_j->property->sequence.molNum))) {
#ifdef HB_TRACE
        if (currenttime > 340.0) printf("[HB_TRACE]   -> FAIL: sequence resNum diff=%d\n", ABSVALUE((HB_i->property->sequence.resNum - HB_j->property->sequence.resNum)));
#endif
        return;
    }
    if (CheckDuplication(HB_i, HB_j, HB_Event)) {
#ifdef HB_TRACE
        if (currenttime > 340.0) printf("[HB_TRACE]   -> FAIL: duplicate event\n");
#endif
        return;
    }
    if (CheckPairMatch(HB_i, HB_j, thisList) < 0) {
#ifdef HB_TRACE
        if (currenttime > 340.0) printf("[HB_TRACE]   -> FAIL: CheckPairMatch (HB.neighbor[0]=%d)\n", HB_i->dynamic->HB.neighbor[0]);
#endif
        return;
    }

    parameters->coreorShell = FindPair(HB_i, HB_j, InterHB, parameters->b_ij, parameters->r_2, &parameters->shortlimit2, &parameters->longlimit2, &parameters->lowerPotential, &parameters->upperPotential, &tmpDouble, NULL, 0);

    if (parameters->b_ij < 0 && parameters->coreorShell == 1 && CalculateDisc(parameters) > 0) {
        parameters->s_time = -1;
        parameters->d_ij2 = parameters->shortlimit2;
        event.potential = parameters->lowerPotential;
        event.subEventType = TBDHBForm;
        event.time = CalculateTime(parameters);
        if (unlikely(event.time < 0)) event.time = INFINIT;
#ifdef HB_TRACE
        if (currenttime > 340.0) printf("[HB_TRACE]   -> SCHEDULED TBDHBForm at t+%.4f (b=%.4f cs=%d disc=%.4f)\n",
               event.time, parameters->b_ij, parameters->coreorShell, CalculateDisc(parameters));
#endif
        JobAssign(HB_i, HB_j, &event, HB_Event);
    }
#ifdef HB_TRACE
    else if (currenttime > 340.0) {
        printf("[HB_TRACE]   -> NO SCHEDULE: b=%.4f cs=%d disc=%.4f r=%.3f\n",
               parameters->b_ij, parameters->coreorShell, CalculateDisc(parameters), sqrt(parameters->r_2));
    }
#endif
}

void HBTime(struct AtomStr *HB_i, struct AtomListStr *thisList) {
    UpdateAtom(HB_i, currenttime);
    if (FindStrInList(HB_i->property->extraProperty, "_HBRMV", 0)) return;
    // HBSkip (a failed formation) only suppresses the search for new partners;
    // the atom's existing HBs must always be re-predicted
    int skipScan = (HB_i->dynamic->HB.interactionType == HBSkip);
    HB_i->dynamic->HB.interactionType = NoEvent;

    for (int p = 1; p <= HBPartnerNum(HB_i); p++) {
        struct AtomStr *HB_j = thisList->ptr[HB_i->dynamic->HB.bondConnection[p]];
        UpdateAtom(HB_j, currenttime);
        if (FindStrInList(HB_j->property->extraProperty, "_HBRMV", 0)) {
            printf("!!ERROR!!: HB_j #%i has been removed! %s:%i\n", HB_i->dynamic->HB.bondConnection[p], __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if (CheckDuplication(HB_i, HB_j, HB_Event)) continue;

        struct ParameterStr parameters;
        TRANSFER_VECTOR(parameters.position_i, HB_i->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_i,    HB_i->dynamic->velocity);
        TRANSFER_VECTOR(parameters.position_j, HB_j->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j,    HB_j->dynamic->velocity);
        CalculateParameters(HB_i, HB_j, &parameters);

        parameters.coreorShell = FindPair(HB_i, HB_j, InterHB, parameters.b_ij, parameters.r_2, &parameters.shortlimit2, &parameters.longlimit2, &parameters.lowerPotential, &parameters.upperPotential, &tmpDouble, NULL, 0);

        struct InteractionEventStr event;
        InitializeEvent(&event);
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) {
            parameters.s_time = -1;
            if (parameters.coreorShell == 0) {
                parameters.d_ij2 = parameters.shortlimit2;
                event.subEventType = CoreColli;
            } else {
                parameters.d_ij2 = parameters.shortlimit2;
                event.subEventType = TBDCaptr;
            }
            event.potential = parameters.lowerPotential;
        } else {
            // Handles b_ij > 0 (outward) AND b_ij < 0 with disc <= 0 (H can't
            // reach the next inner step; must schedule the outer boundary crossing
            // so the velocity adjustment at that wall is not skipped).
            parameters.s_time = 1;
            if (parameters.coreorShell == 1) {
                parameters.d_ij2 = parameters.longlimit2;
                event.subEventType = TBDHBBrik;
            } else {
                parameters.d_ij2 = parameters.longlimit2;
                event.subEventType = TBDEscap;
            }
            event.potential = parameters.upperPotential;
        }
        event.time = CalculateTime(&parameters);
        // on the outer edge already (FindPair keeps a bonded pair there in the last shell): break now
        if (event.time < 0 && parameters.s_time == 1 && ABSVALUE(parameters.r_2 - parameters.d_ij2) < ZERO)
            event.time = 0;
        if (unlikely(event.time < 0)) event.time = INFINIT;
        JobAssign(HB_i, HB_j, &event, HB_Event);
    }
    if (!skipScan && HBHasFreeSlot(HB_i) && HB_i->property->num != 0 && (HB_i->dynamic->HB.role == 'A' || HB_i->dynamic->HB.role == 'D' || HB_i->dynamic->HB.role == 'S')) {
        ScanNeighborhood(HB_i, thisList, HBPairFilter, HBCallback);
    }
}


void HBNeighborTime(struct AtomStr *neighbor_i, struct AtomListStr *thisList) {
    UpdateAtom(neighbor_i, currenttime);
    int flag = 0;
    int target_i, target_j;
    struct AtomStr *neighbor_j;
    struct AtomStr *virtualAtom[3];
    struct InteractionEventStr event;
    struct ParameterStr parameters;
    
    if (FindStrInList(neighbor_i->property->extraProperty, "_HBRMV", 0)) {
        return;
    }
    
    InitializeEvent(&event);
    target_i = neighbor_i->property->num;
    TRANSFER_VECTOR(parameters.position_i, neighbor_i->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i,    neighbor_i->dynamic->velocity);
    
    for (int n = 1; n <= neighbor_i->dynamic->HBNeighbor.neighborStatus; n++) {
        
        target_j   = neighbor_i->dynamic->HBNeighbor.neighborPartner[n];
        neighbor_j = thisList->ptr[target_j];
        
        if (FindStrInList(neighbor_j->property->extraProperty, "_HBRMV", 0)) {
            printf("!!ERROR!!: the neighbor_j #%i has been removed! Something is wrong!\n", target_j);
            printf("           %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        if (flag == 0) {
            neighbor_i->dynamic->HBNeighbor.partnerNum = n; //decide which host will be used to do calculation
        }
        
        UpdateAtom(neighbor_j, currenttime);

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
        int hostPartner;
        if ((neighbor_i->dynamic->HB.role == 'A' ||
             neighbor_i->dynamic->HB.role == 'D' ||
             neighbor_i->dynamic->HB.role == 'S') &&
            (hostPartner = HBNeighborHostPartner(neighbor_i, target_j)) > 0)
        {
            virtualAtom[0] = neighbor_i; //host
            virtualAtom[1] = neighbor_j; //neighbor
            virtualAtom[2] = thisList->ptr[hostPartner];
        } else {
            virtualAtom[0] = neighbor_j;
            virtualAtom[1] = neighbor_i;
            if ((hostPartner = HBNeighborHostPartner(neighbor_j, target_i)) == 0)
                hostPartner = neighbor_j->dynamic->HB.bondConnection[1];
            virtualAtom[2] = thisList->ptr[hostPartner];
        }
        //==============================
        
        parameters.coreorShell = FindPair(virtualAtom[1], virtualAtom[0], InterNeighbor,
                                          parameters.b_ij, parameters.r_2,
                                          &parameters.shortlimit2, &parameters.longlimit2,
                                          &parameters.lowerPotential, &parameters.upperPotential,
                                          &tmpDouble, virtualAtom[2], 0);
        
        if (parameters.b_ij < 0 && CalculateDisc(&parameters) > 0) {
            parameters.s_time  = -1;
            parameters.d_ij2   = parameters.shortlimit2;
            event.potential    = parameters.lowerPotential;
            event.subEventType = (parameters.coreorShell == 0) ? CoreColli : TBDCaptr;
        } else {
            parameters.s_time  = 1;
            parameters.d_ij2   = parameters.longlimit2;
            event.potential    = parameters.upperPotential;
            event.subEventType = (parameters.coreorShell == 1) ? WellBounc : TBDEscap;
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
    UpdateAtom(target, currenttime);
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
    UpdateAtom(target, currenttime);
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
    
    if (FindStrInList(target->property->extraProperty, "_HBRMV", 0)) {
        return;
    }
    
ReCal:
    InitializeEvent(&event);
    TRANSFER_VECTOR(boxSize, boxDimension);
    TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(parameters.speed_i, target->dynamic->velocity);
    
    if (strcmp(surfInfo.style, "implicit") == 0) {
        
        double distance2 = 0;
        struct AtomStr *thisWall = &wall[0];
        TRANSFER_VECTOR(parameters.position_j, thisWall->dynamic->coordinate);
        TRANSFER_VECTOR(parameters.speed_j, thisWall->dynamic->velocity);
        
        if (strncmp(surfInfo.type, "parallel", 1) == 0) { //parallel plates of wall
            parameters.position_i[1] = parameters.position_i[3] = 0;
            parameters.speed_i[1]    = parameters.speed_i[3]    = 0;
        } else if (strncmp(surfInfo.type, "cylinder", 1) == 0) {
            parameters.position_i[1] = 0;
            parameters.speed_i[1]    = 0;
        }
        
        DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
        DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
        parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
        parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
        parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
        distance2 = boxSize[2] * boxSize[2] * 0.25 + parameters.r_2 - boxSize[2] * sqrt(parameters.r_2);
        
        parameters.coreorShell = FindPair(target, thisWall, InterColli,
                                          -1 * parameters.b_ij, distance2,
                                          &parameters.shortlimit2, &parameters.longlimit2,
                                          &parameters.lowerPotential, &parameters.upperPotential,
                                          &tmpDouble, NULL, 0);
        
        if (parameters.lowerPotential == 0 &&
            parameters.upperPotential == 0) {
            
            if (strcmp(wallDyn.mark, "no") != 0) {
                wallDyn.touch ++;
                
                if (strncmp(surfInfo.type, "parallel", 1) == 0) {
                    boxDimension[2] += wallDyn.step;
                } else if (strncmp(surfInfo.type, "cylinder", 1) == 0) {
                    boxDimension[2] += wallDyn.step;
                    boxDimension[3] += wallDyn.step;
                } else if (strncmp(surfInfo.type, "sphere", 1) == 0) {
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
    UpdateAtom(target, currenttime);
    int atomNum = target->property->num;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    if (FindStrInList(target->property->extraProperty, "_HBRMV", 0)) {
        return;
    }
    
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
            
            parameters.coreorShell = FindPair(target, thisWall, InterColli,
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
    UpdateAtom(target, currenttime);
    int atomNum = target->property->num;
    double radius2;
    double distance2 = 0;
    double distanceShift2;
    double extraRange;
    double futurePos;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    struct AtomStr *thisWall = &tunlObj.tunnel;
    struct ConstraintStr *thisConstr = RightPair(  target->property->typeofAtom,
                                                   target->property->typeofAtom[0],
                                                 thisWall->property->typeofAtom,
                                                 thisWall->property->typeofAtom[0],
                                                 0);
    
    if (FindStrInList(target->property->extraProperty, "_HBRMV", 0)) {
        return;
    }
    
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
        
        parameters.coreorShell = FindPair(target, thisWall, InterColli,
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
    UpdateAtom(target, currenttime);
    int flag = 0;
    int gapNum;
    int totalCNum = (flow.charge.PBCMark) ? flow.charge.num : flow.charge.num / 2;
    double charge;
    double distance2, gap2;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    if (FindStrInList(target->property->extraProperty, "_HBRMV", 0)) {
        return;
    }
    
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
        DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
        parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
        parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
        parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
        
        gapNum = (int)(sqrt(parameters.r_2) / flow.charge.gap[chargeNum]);
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
    UpdateAtom(target, currenttime);
    int atomNum = target->property->num;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    if (FindStrInList(target->property->extraProperty, "_HBRMV", 0)) {
        return;
    }
    
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
    event->flag         = 1;
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

/*
 Backup events (Optimization_Report.md Section 12). Besides its earliest event,
 each atom keeps the next BACKUP_EVENTS best events of its last prediction.
 When the earliest event expires because the partner took part in another
 event (the atom itself is unchanged), the earliest backup whose partner is
 also unchanged is the atom's true next event, so the full re-prediction can
 be skipped. Events with changed partners are re-predicted from the partner's
 side. Thermostat events are not kept: a fresh one is drawn exactly as a full
 re-prediction would.
 */
#define BACKUP_EVENTS 8

struct BackupListStr {
    int num;
    struct InteractionEventStr event[BACKUP_EVENTS];
};

static struct BackupListStr *backupEvents = NULL;
static int backupUsable = 0;

void InitBackupEvents(void) {
    free(backupEvents);
    backupEvents = (struct BackupListStr *)calloc(atomnum + 1, sizeof(struct BackupListStr));
    // events of walls, obstacles, tunnels, charges and CG surfaces encode their
    // partner after JobAssign, so backups are only used without them
    backupUsable = !(strcmp(surfInfo.style, "implicit") == 0 || obstObj.mark || tunlObj.mark ||
                     flow.mark == 3 || SphObstObj.mark || CG.mark);
}

void ShiftBackupEvents(double shift) {
    if (backupEvents == NULL) return;
    for (int i = 1; i <= atomnum; i++)
        for (int k = 0; k < backupEvents[i].num; k++)
            backupEvents[i].event[k].time -= shift;
}

void ClearBackupEvents(int atomNum) {
    if (backupEvents) backupEvents[atomNum].num = 0;
}

// TRUE if an HB-formation event (type/sub/partner/time) must be ordered before
// the existing event: a collision with the same partner at the same time
static int HBFormWinsTie(enum EEventType type, enum EInteractionType sub, int partner, double time,
                         const struct InteractionEventStr *existing) {
    return type == HB_Event && sub == TBDHBForm &&
           existing->eventType == Coli_Event && existing->partner == partner &&
           fabs(existing->time - time) <= 1e-9 * (1 + time);
}

static void KeepBackup(struct AtomStr *atom_i, const struct InteractionEventStr *event) {
    if (backupEvents == NULL || event->eventType == Ther_Event ||
        event->eventType == Invd_Event || event->time >= INFINIT) return;

    struct BackupListStr *list = &backupEvents[atom_i->property->num];
    if (list->num == BACKUP_EVENTS && event->time >= list->event[BACKUP_EVENTS - 1].time) return;
    if (list->num < BACKUP_EVENTS) list->num++;

    int pos = list->num - 1;
    while (pos > 0 && (list->event[pos - 1].time > event->time ||
                       HBFormWinsTie(event->eventType, event->subEventType, event->partner, event->time, &list->event[pos - 1]))) {
        list->event[pos] = list->event[pos - 1];
        pos--;
    }
    list->event[pos] = *event;
}

int TakeBackupEvent(struct AtomStr *target, struct AtomListStr *thisList) {
    if (!backupUsable || backupEvents == NULL) return FALSE;

    struct BackupListStr *list = &backupEvents[target->property->num];
    for (int k = 0; k < list->num; k++) {
        struct InteractionEventStr *backup = &list->event[k];
        if (backup->partner > 0 &&
            thisList->ptr[backup->partner]->dynamic->event.counter != backup->partnerCounter) continue;

        target->dynamic->event.partner        = backup->partner;
        target->dynamic->event.partnerCounter = backup->partnerCounter;
        target->dynamic->event.time           = backup->time;
        target->dynamic->event.potential      = backup->potential;
        target->dynamic->event.eventType      = backup->eventType;
        target->dynamic->event.subEventType   = backup->subEventType;

        list->num -= k + 1;
        memmove(list->event, list->event + k + 1, list->num * sizeof(struct InteractionEventStr));

        UpdateAtom(target, currenttime);
        ThermostatTime(target);
        return TRUE;
    }
    return FALSE;
}

int CheckDuplication(struct AtomStr *atom_i, struct AtomStr *atom_j, enum EEventType type) {
    if (   atom_j->dynamic->event.partner        == atom_i->property->num
        && atom_j->dynamic->event.eventType      == type
        && atom_j->dynamic->event.partnerCounter == atom_i->dynamic->event.counter
        && atom_i->dynamic->event.time           >  atom_j->dynamic->event.time) {
        
        KeepBackup(atom_i, &atom_i->dynamic->event);
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
	int atom_i_type = atom_i->property->typeofAtom[1];
	int atom_j_type = atom_j->property->typeofAtom[1];
	if (potentialPairHB[type][atom_i_type][atom_j_type].step == NULL) {
		return -1;
	}

    int neighbor_i_type, neighbor_j_type, count = 0;
    for (int i = 1; i <= atom_i->dynamic->HB.neighbor[0]; i ++) {
        neighbor_i_type = thisList->ptr[atom_i->dynamic->HB.neighbor[i]]->property->typeofAtom[1];
        for (int j = 1; j <= atom_j->dynamic->HB.neighbor[0]; j ++) {
            neighbor_j_type = thisList->ptr[atom_j->dynamic->HB.neighbor[j]]->property->typeofAtom[1];
            if (potentialPairHB[type][atom_i_type][neighbor_j_type].step != NULL &&
                potentialPairHB[type][atom_j_type][neighbor_i_type].step != NULL) {
                count ++;
            }
        }
    }
    if (count == 0)
        return -1;

	return count;
}

int JobAssign(struct AtomStr *atom_i, struct AtomStr *atom_j, struct InteractionEventStr *event, enum EEventType type) {
    if (event->time <= 0) return FALSE;
    double abs_time = event->time;
    if (abs_time != INFINIT) {
        abs_time += currenttime;
    }
    
    // HB formation and the same pair's nonbonded step at the HB boundary are one
    // physical event: the formation is attempted first, and a failed attempt is
    // processed as the collision (FailedHBToCollision)
    if (atom_i->dynamic->event.time > abs_time ||
        (atom_j && HBFormWinsTie(type, event->subEventType, atom_j->property->num, abs_time, &atom_i->dynamic->event))) {
        KeepBackup(atom_i, &atom_i->dynamic->event);
        atom_i->dynamic->event.potential    = event->potential;
        atom_i->dynamic->event.subEventType = event->subEventType;
        atom_i->dynamic->event.time         = abs_time;
        atom_i->dynamic->event.eventType    = type;
        
        if (atom_j) {
            atom_i->dynamic->event.partner        = atom_j->property->num;
            atom_i->dynamic->event.partnerCounter = atom_j->dynamic->event.counter;
        } else {
            atom_i->dynamic->event.partner = INVALID;
        }
        
        return TRUE;
    }

    struct InteractionEventStr offered = *event;
    offered.time           = abs_time;
    offered.eventType      = type;
    offered.partner        = atom_j ? atom_j->property->num : INVALID;
    offered.partnerCounter = atom_j ? atom_j->dynamic->event.counter : 0;
    KeepBackup(atom_i, &offered);
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
