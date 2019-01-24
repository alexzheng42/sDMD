//
//  Events.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

int HBEvent(struct AtomStr **atomLibrary, struct AtomStr *HB_i, struct AtomStr *HB_j, struct AtomStr *neighbor_i, struct AtomStr *neighbor_j);
void ThermostatEvent(struct AtomStr *targetAtom, int threadID);
void PBCandCrossCellEvent(struct AtomStr *targetAtom, struct ThreadStr *thisThread);
void WallEvent(struct AtomStr *targetAtom);
void ObstEvent(struct AtomStr *targetAtom);
void TunnelEvent(struct AtomStr *targetAtom);
void ChargeEvent(struct AtomStr *targetAtom);
void SphObstEvent(struct AtomStr *targetAtom);
void InteractionEvent(struct AtomStr *targetAtom, struct AtomStr *partner, struct AtomStr *HBNeighbor_i, struct AtomStr *HBNeighbor_j, struct ThreadStr* thisThread, list *atomList);
double HBEnergyChange(struct AtomStr *targetAtom, int skipAtom, char *type, list *atomList, struct ThreadStr *thisThread);


int DoEvent(struct ThreadStr *thisThread) {
    struct AtomStr *targetAtom, *partner;
    struct AtomStr *HBNeighbor_i, *HBNeighbor_j;
    int typenum, interactionType = 0;
    
    targetAtom = thisThread->newTarget;
    partner = thisThread->newPartner;
    typenum = targetAtom->dynamic->event.eventType;
    
    //switch to the event which will happen earliest
    switch (typenum) {
            
        case Coli_Event: //do collision event
        case Bond_Event: //do bond event
        case HBNe_Event: //do HB neighbor event
            InteractionEvent(targetAtom, partner, NULL, NULL, thisThread, &thisThread->atomList);
            break;
            
        case HB_Event: //do HB event
            HBNeighbor_i = thisThread->listPtr[targetAtom->dynamic->HB.neighbor];
            HBNeighbor_j = thisThread->listPtr[partner->dynamic->HB.neighbor];
            InteractionEvent(targetAtom, partner, HBNeighbor_i, HBNeighbor_j, thisThread, &thisThread->atomList);
            interactionType = HBEvent(thisThread->listPtr, targetAtom, partner, HBNeighbor_i, HBNeighbor_j);
            break;
            
        case Ther_Event: //do thermostat ghost collision
            ThermostatEvent(targetAtom, thisThread->tid);
            break;
            
        case PCC_Event: //do PCC event
            PBCandCrossCellEvent(targetAtom, thisThread);
            break;
            
        case Lagv_Event: //do Langevin event
            break;
            
        case Wall_Event: //do wall collision
            WallEvent(targetAtom);
            break;
            
        case Obst_Event: //do obstruction collision
            ObstEvent(targetAtom);
            break;
            
        case Tunl_Event: //do tunnel collision
            TunnelEvent(targetAtom);
            break;
            
        case Chrg_Event: //do charge interaction
            ChargeEvent(targetAtom);
            break;
            
        case SphO_Event: //do column interaction
            SphObstEvent(targetAtom);
            break;
            
        case Invd_Event:
        case Cacl_Event:
        case  TBD_Event:
            printf("!!ERROR!!: do_event switch number has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            
        default:
            break;
            
    }
    
    return interactionType; //if it is a HB event, return HB interaction type (HB forming or breaking)
}


void LinkList(char * type, struct AtomStr *targetAtom, struct ThreadStr *thisThread) {
    int i, n;
    int cellindex;
    double *boxsize;
    double *position;
    
    POINT_TO_STRUCT(boxsize, boxDimension);
    
    if (unlikely(strncmp(type, "full", 1) == 0 && targetAtom == NULL)) {
        
        //-------------------------------------------
        //assign a group of memory to store the link list
        INT_CALLOC(celllist, (atomnum + cellnum[1] * cellnum[2] * cellnum[3] + 1));
        
        //-------------------------------------------
        //calculate which subcells the atoms are in and generate the list
        for (i = 1; i <= atomnum; i++) {
            
            POINT_TO_STRUCT(position, thisThread->raw[i]->dynamic->coordinate);
            
            for (n = 1; n <= 3; n++) {
                
                thisThread->raw[i]->dynamic->cellIndex[n] = position[n] / cellsize[n]; //include 0
                //atom on the edge of box size has been removed by function pbc()
                //integer numbers
                //for a 3*3*3 system, cell index on one side would be:
                //  [0)    [1)   [2)
                //  [0)    [1)   [2)
                //  [0)    [1)   [2)
                
                //just in case the floating-point arithmetic error
                //the PBC-adjusted coordinates may be exactly equal to the upper edge of the PBC box
                if (thisThread->raw[i]->dynamic->cellIndex[n] == cellnum[n]) {
                    thisThread->raw[i]->dynamic->cellIndex[n] --;
                }
            }
            
            cellindex = thisThread->raw[i]->dynamic->cellIndex[3] * cellnum[1] * cellnum[2]
                      + thisThread->raw[i]->dynamic->cellIndex[2] * cellnum[1]
                      + thisThread->raw[i]->dynamic->cellIndex[1] + 1;
            
            celllist[i] = celllist[atomnum + cellindex];
            celllist[atomnum + cellindex] = i;
            
        }
        //-------------------------------------------
        
    } else if (likely(strncmp(type, "part", 1) == 0 && targetAtom != NULL)) {
        
        POINT_TO_STRUCT(position, targetAtom->dynamic->coordinate);
        for (n = 1; n <= 3; n++) {
            targetAtom->dynamic->cellIndex[n] = position[n] / cellsize[n];
            
            if (unlikely(targetAtom->dynamic->cellIndex[n] == cellnum[n])) {
                targetAtom->dynamic->cellIndex[n]--;
            }
        }
        
    } else {
        printf("Link list arguments error!\n");
    }
}


void PBC(char * type, struct AtomStr *targetAtom, struct ThreadStr *thisThread) {
    int dim = 3, step = 1;
    
    if (strcmp(wallExist, "no")) {
        if (strncmp(wallType, "sphere", 1) == 0) return;
        if (strncmp(wallType, "parallel", 1) == 0)
            step = 2;
        else if (strncmp(wallType, "cylinder", 1) == 0)
            dim = 1;
    }
    
    if (unlikely(strncmp(type, "full", 1) == 0 && targetAtom == NULL)) {
        
        for (int i = 1; i <= atomnum; i++) {
            for (int n = 1; n <= dim; n += step) {
                if (thisThread->raw[i]->dynamic->coordinate[n] >= boxDimension[n]) { //coordinate range [0,box size)
                    thisThread->raw[i]->dynamic->coordinate[n] -= boxDimension[n];
                } else if (thisThread->raw[i]->dynamic->coordinate[n] < 0) {
                    thisThread->raw[i]->dynamic->coordinate[n] += boxDimension[n];
                }
            }
        }
        
    } else if (likely(strncmp(type, "part", 1) == 0 && targetAtom != NULL)) {
        
        for (int n = 1; n <= dim; n += step) {            
            if (unlikely(targetAtom->dynamic->coordinate[n] >= boxDimension[n])) { //coordinate range [0,box size)
                targetAtom->dynamic->coordinate[n] -= boxDimension[n];
            } else if (unlikely(targetAtom->dynamic->coordinate[n] < 0)) {
                targetAtom->dynamic->coordinate[n] += boxDimension[n];
            }
        }
        
    } else {
        
        printf("PBC arguments error!\n");
        
    }
}


void PBCandCrossCellEvent(struct AtomStr *targetAtom, struct ThreadStr *thisThread) {
    int n, original[4];
    
    original[1] = targetAtom->dynamic->cellIndex[1];
    original[2] = targetAtom->dynamic->cellIndex[2];
    original[3] = targetAtom->dynamic->cellIndex[3];
    
    for (n = 1; n <= 3; n++) {
        if (unlikely(targetAtom->dynamic->cellIndex[n] == 0 ||
            targetAtom->dynamic->cellIndex[n] == cellnum[n] - 1)) {
            PBC("part", targetAtom, thisThread);
            break;
        }
    }
    
    LinkList("part", targetAtom, thisThread);
    
    if (unlikely(original[1] == targetAtom->dynamic->cellIndex[1] &&
                 original[2] == targetAtom->dynamic->cellIndex[2] &&
                 original[3] == targetAtom->dynamic->cellIndex[3] &&
                 targetAtom->dynamic->event.time >= ZERO * 100)) {
        printf("check PCC calculation! atom: %i time: %lf\n", targetAtom->property->sequence.atomNum, targetAtom->dynamic->event.time);
    }
    
    return;
}


int HBEvent(struct AtomStr **atomLibrary, struct AtomStr *HB_i, struct AtomStr *HB_j, struct AtomStr *neighbor_i, struct AtomStr *neighbor_j) {
    int temp, targetAtom_i, targetAtom_j;
    
    targetAtom_i = HB_i->property->num;
    targetAtom_j = HB_j->property->num;
    
    if (HB_i->dynamic->HB.interactionType == HBForm) {
        
        HB_i->dynamic->HB.bondConnection = targetAtom_j;
        HB_j->dynamic->HB.bondConnection = targetAtom_i;
        
        HB_i->property->type = AtomTypeChange(HB_i->property->type, 1);
        HB_j->property->type = AtomTypeChange(HB_j->property->type, 1);
        
#ifdef VIS
        if (visual) {
            ChangeColor(HB_i->property->type, HB_i->property->color);
            ChangeColor(HB_j->property->type, HB_j->property->color);
        }
#endif
        
        //===================================
        //for HB acceptor and donor
        temp = ++HB_i->dynamic->HBNeighbor.neighborStatus;
        HB_i->dynamic->HBNeighbor.neighborPartner[temp] = HB_j->dynamic->HB.neighbor;
        
        temp = ++HB_j->dynamic->HBNeighbor.neighborStatus;
        HB_j->dynamic->HBNeighbor.neighborPartner[temp] = HB_i->dynamic->HB.neighbor;
        
        //for neighbors
        temp = ++neighbor_i->dynamic->HBNeighbor.neighborStatus;
        neighbor_i->dynamic->HBNeighbor.neighborPartner[temp] = targetAtom_j;
    
        temp = ++neighbor_j->dynamic->HBNeighbor.neighborStatus;
        neighbor_j->dynamic->HBNeighbor.neighborPartner[temp] = targetAtom_i;
        
    } else if (HB_i->dynamic->HB.interactionType == HBBreak) {
        
        HB_i->dynamic->HB.bondConnection = 0;
        HB_j->dynamic->HB.bondConnection = 0;
        
        HB_i->property->type = AtomTypeChange(HB_i->property->type, 0);
        HB_j->property->type = AtomTypeChange(HB_j->property->type, 0);
        
#ifdef VIS
        if (visual) {
            ChangeColor(HB_i->property->type, HB_i->property->color);
            ChangeColor(HB_j->property->type, HB_j->property->color);
        }
#endif
        
        //===================================
        //for HB acceptor and donor
        --HB_i->dynamic->HBNeighbor.neighborStatus;
        ListRefresh(HB_j->dynamic->HB.neighbor, HB_i->dynamic->HBNeighbor.neighborPartner, 1, 4);
        
        --HB_j->dynamic->HBNeighbor.neighborStatus;
        ListRefresh(HB_i->dynamic->HB.neighbor, HB_j->dynamic->HBNeighbor.neighborPartner, 1, 4);
        
        //for neighbor
        --neighbor_i->dynamic->HBNeighbor.neighborStatus;
        ListRefresh(targetAtom_j, neighbor_i->dynamic->HBNeighbor.neighborPartner, 1, 4);
        
        --neighbor_j->dynamic->HBNeighbor.neighborStatus;
        ListRefresh(targetAtom_i, neighbor_j->dynamic->HBNeighbor.neighborPartner, 1, 4);
    }
    
    return HB_i->dynamic->HB.interactionType;
}


void ThermostatEvent(struct AtomStr *targetAtom, int threadID) {
    double speed[4];
    double netV[4] = {0};

    if (chargeAA[COM] && flow.mark == 3) { //!!not parallizable!! will use the raw data
        CalCOMV(targetAtom->property->sequence.proteinNum, netV);
    }
    
    if (strncmp(thermostatType, "Andersen", 1) == 0) {
        for (int n = 1; n <= 3; n++) {
        regenerate:
            speed[n] = RandomVelocity(targetAtom->property->mass);
            speed[n] += flow.constV.v[n] + netV[n] +
                        flow.force.a[targetAtom->property->num][n] *
                        (currenttime - flow.force.timeRec[targetAtom->property->num]);
            
            if (unlikely(speed[n] == 0)) {
                goto regenerate;
            }
        }
        
        TRANSFER_VECTOR(targetAtom->dynamic->velocity, speed);
        flow.force.timeRec[targetAtom->property->num] = currenttime;
    }
}


void WallEvent(struct AtomStr *targetAtom) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, *speed_j = NULL;
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisWall = NULL;
    
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    speed_i = targetAtom->dynamic->velocity;
    
    if (strcmp(wallExist, "smooth") == 0) {
        thisWall = &wall[0];
        TRANSFER_VECTOR(position_j, thisWall->dynamic->coordinate);
        speed_j = thisWall->dynamic->velocity;
        
        if (strncmp(wallType, "parallel", 1) == 0) {
            position_i[1] = position_i[3] = 0;
        } else if (strncmp(wallType, "cylinder", 1) == 0) {
            position_i[1] = 0;
        }
    }
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = targetAtom->property->mass;
    potential = targetAtom->dynamic->event.potential;
    
    switch (targetAtom->dynamic->event.subEventType) {
        case WellBounc: //well bounce
            s_vel = -1;
            potential = 0;
            break;
            
        case TBDCaptr:
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to capture?
                // >0 -> sufficient; <=0 -> no
                potential = 0;
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            break;
            
        case TBDEscap: //TBD outer
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to escape?
                potential = 0;
                s_vel = -1;
            } else {
                s_vel = 1;
            }
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / targetAtom->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, targetAtom->dynamic->velocity);

    return;
}


void DoWallDyn(void) {
    if (currenttime >= wallDyn.curTime) {
        if (strcmp(wallDyn.mark, "pulse") == 0) {
            if ((wallDyn.origBoxDim[2] - boxDimension[2] >= wallDyn.size) ||
                (wallDyn.origBoxDim[2] == boxDimension[2])) {
                wallDyn.sign *= -1;
                wallDyn.touch = 0;
            } else if (fabs(boxDimension[2] + (int)(currenttime / wallDyn.rate + 1) * wallDyn.size - (currenttime - 1) * wallDyn.step - wallDyn.origBoxDim[2]) <= wallDyn.step &&
                       wallDyn.touch >= 5) {  //balance the compressing and expanding time
                wallDyn.sign *= -1;
                wallDyn.touch = 0;
            }
        }
        
        wallDyn.curTime += 1;
        
        if (strcmp(wallType, "parallel") == 0) {
            boxDimension[2] += wallDyn.sign * wallDyn.step;
        } else if (strcmp(wallType, "cylinder") == 0) {
            boxDimension[2] += wallDyn.sign * wallDyn.step;
            boxDimension[3] += wallDyn.sign * wallDyn.step;
        } else if (strcmp(wallType, "sphere") == 0) {
            boxDimension[1] += wallDyn.sign * wallDyn.step;
            boxDimension[2] += wallDyn.sign * wallDyn.step;
            boxDimension[3] += wallDyn.sign * wallDyn.step;
        }
        
        for (int i = 1; i <= atomnum; i ++) {
            if (atom[i].dynamic->event.eventType == Wall_Event) {
                atom[i].dynamic->event.counter ++;
            }
        }
        
        if (((wallDyn.origBoxDim[2] - boxDimension[2] >= wallDyn.size) && strcmp(wallDyn.mark, "compress") == 0) ||
            ((wallDyn.origBoxDim[2] - boxDimension[2] <= wallDyn.step) && strcmp(wallDyn.mark, "expand")   == 0)) {
            sprintf(wallDyn.mark, "no");
        }
    }
}


void ObstEvent(struct AtomStr *targetAtom) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, *speed_j = NULL;
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisWall = &obstObj.obst[targetAtom->dynamic->event.partner / -10];
    
    int dim = (targetAtom->dynamic->event.partner * -1) % 10;
    
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    TRANSFER_VECTOR(position_j,   thisWall->dynamic->coordinate);

    speed_i = targetAtom->dynamic->velocity;
    speed_j =   thisWall->dynamic->velocity;
    
    for (int i = 1; i <= 3; i ++) {
        if (dim != i) {
            position_i[i] = 0;
            position_j[i] = 0;
        }
    }
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(   speed_i,    speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij  = DOT_PROD(r_ij, v_ij);
    reducedMass = targetAtom->property->mass;
    potential = targetAtom->dynamic->event.potential;
    
    switch (targetAtom->dynamic->event.subEventType) {
        case CoreColli:
            s_vel = 1;
            potential = 0;
            break;
            
        case TBDCaptr:
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to capture?
                // >0 -> sufficient; <=0 -> no
                potential = 0;
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            break;
            
        case TBDEscap: //TBD outer
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to escape?
                potential = 0;
                s_vel = -1;
            } else {
                s_vel = 1;
            }
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / targetAtom->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, targetAtom->dynamic->velocity);
    
    return;
}


void TunnelEvent(struct AtomStr *targetAtom) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, *speed_j = NULL;
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisWall = &tunlObj.tunnel;
    
    int num = targetAtom->dynamic->event.partner / -10000;
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, tunlObj.position[num]);
    speed_i = targetAtom->dynamic->velocity;
    speed_j =   thisWall->dynamic->velocity;
    position_i[1] = position_j[1]; //end position

    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = targetAtom->property->mass;
    potential = targetAtom->dynamic->event.potential;
    
    switch (targetAtom->dynamic->event.subEventType) {
        case WellBounc: //well bounce
            s_vel = -1;
            potential = 0;
            break;
            
        case TBDCaptr:
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to capture?
                // >0 -> sufficient; <=0 -> no
                potential = 0;
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            break;
            
        case TBDEscap: //TBD outer
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to escape?
                potential = 0;
                s_vel = -1;
            } else {
                s_vel = 1;
            }
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / targetAtom->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, targetAtom->dynamic->velocity);
    
    return;
}


void ChargeEvent(struct AtomStr *targetAtom) {
    int s_vel = 0;
    int chargeNum = targetAtom->dynamic->event.partner / -100;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, *speed_j = NULL;
    double position_i[4] = {0}, position_j[4] = {0};
    
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, flow.charge.position[chargeNum]);
    speed_i = targetAtom->dynamic->velocity;
    speed_j = flow.charge.velocity[chargeNum];
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = targetAtom->property->mass;
    potential = targetAtom->dynamic->event.potential;
    
    switch (targetAtom->dynamic->event.subEventType) {
        case TBDCaptr:
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to capture?
                // >0 -> sufficient; <=0 -> no
                potential = 0;
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            break;
            
        case TBDEscap: //TBD outer
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to escape?
                potential = 0;
                s_vel = -1;
            } else {
                s_vel = 1;
            }
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / targetAtom->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, targetAtom->dynamic->velocity);
    
    return;
}


void SphObstEvent(struct AtomStr *targetAtom) {
    int s_vel = 0;
    int thisObst = targetAtom->dynamic->event.partner / -1000;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, SphObstObj.position[thisObst]);
    
    speed_i = targetAtom->dynamic->velocity;
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(   speed_i,    speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij  = DOT_PROD(r_ij, v_ij);
    reducedMass = targetAtom->property->mass;
    
    switch (targetAtom->dynamic->event.subEventType) {
        case CoreColli:
            s_vel = 1;
            potential = 0;
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / targetAtom->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, targetAtom->dynamic->velocity);
    
    return;
}


void InteractionEvent (struct AtomStr *targetAtom, struct AtomStr *partner, struct AtomStr *HBNeighbor_i, struct AtomStr *HBNeighbor_j, struct ThreadStr* thisThread, list *atomList)
{
    int s_vel = 0;
    double v_ij[4];
    double b_ij, r_ij[4];
    double term[4], temp;
    double phi, potential, reducedMass;
    double d_ij2;
    double *speed_i, *speed_j;
    double position_i[4], position_j[4];
    double neighborLowerPotential, neighborUpperPotential;
    
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, partner->dynamic->coordinate);
    speed_i = targetAtom->dynamic->velocity;
    speed_j = partner->dynamic->velocity;
    
    for (int n = 1; n <= 3; n++) {
        if (targetAtom->dynamic->cellIndex[n] != partner->dynamic->cellIndex[n]) {
            if (targetAtom->dynamic->cellIndex[n] == 0 && partner->dynamic->cellIndex[n] == cellnum[n] - 1) {
                position_j[n] -= boxDimension[n];
            } else if (targetAtom->dynamic->cellIndex[n] == cellnum[n] - 1 && partner->dynamic->cellIndex[n] == 0) {
                position_j[n] += boxDimension[n];
            }
        }
    } //apply pbc condition
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = targetAtom->property->mass * partner->property->mass / (targetAtom->property->mass + partner->property->mass);
    potential = targetAtom->dynamic->event.potential;
    
    if (unlikely(targetAtom->dynamic->event.eventType    == HB_Event && // HB event only
                (targetAtom->dynamic->event.subEventType == TBDHBForm || targetAtom->dynamic->event.subEventType == TBDHBBrik))) { //may form or break HB
        
        double position_neighbor_i[4];
        double position_neighbor_j[4];
        double r_c12_neighbor_i[4], r_neighbor_i2;
        double r_c12_neighbor_j[4], r_neighbor_j2;
        double accumPoten[6], tmp = 0;
        
        //-------------------------
        //adjust neighbor PBC
        TRANSFER_VECTOR(position_neighbor_i, HBNeighbor_i->dynamic->coordinate);
        TRANSFER_VECTOR(position_neighbor_j, HBNeighbor_j->dynamic->coordinate);
        
        for (int n = 1; n <= 3; n++) {
            
            //------------------
            //for HB_i
            if (targetAtom->dynamic->cellIndex[n] != HBNeighbor_i->dynamic->cellIndex[n]) {
                if (targetAtom->dynamic->cellIndex[n] == 0 && HBNeighbor_i->dynamic->cellIndex[n] == cellnum[n] - 1) {
                    position_neighbor_i[n] -= boxDimension[n];
                } else if (targetAtom->dynamic->cellIndex[n] == cellnum[n] - 1 && HBNeighbor_i->dynamic->cellIndex[n] == 0) {
                    position_neighbor_i[n] += boxDimension[n];
                }
            }
            //------------------
            
            
            //------------------
            //for HB_j
            //all based on the coordinates of HB_i
            if (targetAtom->dynamic->cellIndex[n] != HBNeighbor_j->dynamic->cellIndex[n]) {
                if (targetAtom->dynamic->cellIndex[n] == 0 && HBNeighbor_j->dynamic->cellIndex[n] == cellnum[n] - 1) {
                    position_neighbor_j[n] -= boxDimension[n];
                } else if (targetAtom->dynamic->cellIndex[n] == cellnum[n] - 1 && HBNeighbor_j->dynamic->cellIndex[n] == 0) {
                    position_neighbor_j[n] += boxDimension[n];
                }
            }
            //------------------
        }
        
        //measure the distances from the neighbors and calculate the total energy change
        DOT_MINUS(position_i, position_neighbor_j, r_c12_neighbor_i);
        DOT_MINUS(position_j, position_neighbor_i, r_c12_neighbor_j);
        r_neighbor_i2 = DOT_PROD(r_c12_neighbor_i, r_c12_neighbor_i);   r_neighbor_j2 = DOT_PROD(r_c12_neighbor_j, r_c12_neighbor_j);
        
        accumPoten[0] = HBEnergyChange(targetAtom,    partner->property->num, "_nonHB", atomList, thisThread);
        accumPoten[1] = HBEnergyChange(partner,    targetAtom->property->num, "_nonHB", atomList, thisThread);
        
        accumPoten[2] = HBEnergyChange(targetAtom,    partner->property->num, "_HB", atomList, thisThread);
        accumPoten[3] = HBEnergyChange(partner,    targetAtom->property->num, "_HB", atomList, thisThread);

        FindPair(HBNeighbor_j, targetAtom, "neighbor", tmp, r_neighbor_i2, &temp, &temp, &neighborLowerPotential, &neighborUpperPotential, &accumPoten[4], partner, 0);
        FindPair(HBNeighbor_i, partner,    "neighbor", tmp, r_neighbor_j2, &temp, &temp, &neighborLowerPotential, &neighborUpperPotential, &accumPoten[5], targetAtom, 0);
        
        if (targetAtom->dynamic->event.subEventType == TBDHBForm) { //form
            potential -= (HBModel(targetAtom, partner) == 1) ? HBPotential.BB : HBPotential.BS;
            potential += (accumPoten[2] - accumPoten[0] + accumPoten[3] - accumPoten[1] + accumPoten[4] + accumPoten[5]);
        } else if (targetAtom->dynamic->event.subEventType == TBDHBBrik) { //break
            potential += (HBModel(targetAtom, partner) == 1) ? HBPotential.BB : HBPotential.BS;
            potential += (accumPoten[0] - accumPoten[2] + accumPoten[1] - accumPoten[3] - accumPoten[4] - accumPoten[5]);
        } else {
            printf("!ERROR!: HB event type has something wrong!\n");
            exit(EXIT_FAILURE);
        }
        
        if (potential < -0.5 * INFINIT && targetAtom->dynamic->event.subEventType == TBDHBBrik) { //means this neighbortime must be shorter, i.e. there will be a neighbor event before the current HB event
            targetAtom->dynamic->HB.interactionType = HBSkip;
            return;
        } else if (potential > 0.5 * INFINIT && targetAtom->dynamic->event.subEventType == TBDHBForm) {
            //try to form HB, but there are neighbor atoms that not at the right positions
            targetAtom->dynamic->HB.interactionType = HBSkip;
            return;
        }
        
        if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass > 0) {
            
            if (targetAtom->dynamic->event.subEventType == TBDHBForm) {
                targetAtom->dynamic->event.subEventType = WellCaptr;
                targetAtom->dynamic->HB.interactionType = HBForm;
            } else if (targetAtom->dynamic->event.subEventType == TBDHBBrik) {
                targetAtom->dynamic->event.subEventType = WellEscap;
                targetAtom->dynamic->HB.interactionType = HBBreak;
            } else {
                printf("!!ERROR!!: should not reach here! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
        } else {
            
            //not enough:
            if (targetAtom->dynamic->event.subEventType == TBDHBForm) {
                targetAtom->dynamic->HB.interactionType = HBSkip;
                return;
            } else if (targetAtom->dynamic->event.subEventType == TBDHBBrik) {
                potential = 0;
                targetAtom->dynamic->event.subEventType = WellBounc;
                targetAtom->dynamic->HB.interactionType = NoEvent;
            } else {
                printf("!!ERROR!!: should not reach here! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    switch (targetAtom->dynamic->event.subEventType) {
        case CoreColli: //core collision
            s_vel = 1;
            potential = 0;
            break;
            
        case WellBounc: //well bounce
            if (b_ij < 0) {
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            potential = 0;
            break;
            
        case WellCaptr: //well capture
            s_vel = -1;
            break;
            
        case WellEscap: //well escape
            if (b_ij < 0) {
                s_vel = -1;
            } else {
                s_vel = 1;
            }
            break;
            
        case TBDCaptr: //TBD inner
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to capture?
                // >0 -> sufficient; <=0 -> no
                potential = 0;
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            break;
            
        case TBDEscap: //TBD outer
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to escape?
                potential = 0;
                if (b_ij < 0) {
                    s_vel = 1;
                } else {
                    s_vel = -1;
                }
            } else {
                if (b_ij < 0) {
                    s_vel = -1;
                } else {
                    s_vel = 1;
                }
            }
            break;
            
        default:
            printf("!ERROR!: InteractionEvent input event type error!\n");
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt((b_ij * b_ij - 2 * d_ij2 * potential / reducedMass))) /
    ((targetAtom->property->mass + partner->property->mass) * d_ij2);
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD((partner->property->mass * phi), r_ij, term);
    DOT_PLUS(speed_i, term, targetAtom->dynamic->velocity);
    FACTOR_PROD((targetAtom->property->mass * phi), r_ij, term);
    DOT_MINUS(speed_j, term, partner->dynamic->velocity);
    
    return;
}


double HBEnergyChange(struct AtomStr *targetAtom, int skipAtom, char *type, list *atomList, struct ThreadStr *thisThread) {
    int atom_i, atom_j;
    double cutoffRange = cutoffr * cutoffr;
    double energyChange = 0;
    double position_i[4], position_j[4];
    double r_temp[4], positionshift[4] = {0};
    double r_ij[4], r_2;
    double tmp = 0, accumPoten = 0;
    struct AtomStr *partner;
    listElem *elem = listFirst(atomList);
    
    atom_i = targetAtom->property->num;
    TRANSFER_VECTOR(position_i, targetAtom->dynamic->coordinate);
    
    for (int i = 1; i <= atomList->num_members; ++i) {
        partner = (struct AtomStr *)elem->obj;
        atom_j = partner->property->num;
        
        if (atom_i != atom_j &&
            !(connectionMap[atom_i][atom_j] & BOND_CONNECT) &&
            !(connectionMap[atom_i][atom_j] & CONSTRAINT_CONNECT) &&
            atom_j != skipAtom) {
            
            TRANSFER_VECTOR(position_j, partner->dynamic->coordinate);
            PBCShift(targetAtom, partner, positionshift);

            DOT_PLUS(position_j, positionshift, r_temp);
            DOT_MINUS(position_i, r_temp, r_ij); //here r_ij is before the collision
            r_2 = DOT_PROD(r_ij, r_ij);
        
            if (r_2 <= cutoffRange) {
                if (strcmp(type, "_nonHB") == 0) {
                    FindPair(targetAtom, partner, "collision", tmp, r_2, &tmp, &tmp, &tmp, &tmp, &accumPoten, NULL, -1);
                } else if (strcmp(type, "_HB") == 0) {
                    if (atom_j != thisThread->raw[skipAtom]->dynamic->HB.neighbor) {
                        FindPair(targetAtom, partner, "collision", tmp, r_2, &tmp, &tmp, &tmp, &tmp, &accumPoten, NULL, 1);
                    }
                } else {
                    printf("!!ERROR!!: type is invalid! %s:%i\n", __FILE__, __LINE__);
                }
                if (accumPoten == INFINIT) accumPoten = 0;
                energyChange += accumPoten;
            }
        }
        elem = listNext(atomList, elem);
    }
    
    return energyChange;
}
