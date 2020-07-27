//
//  Events.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

int CheckAlphaHB(struct AtomStr* HB_i, struct AtomStr* HB_j);
void HBEvent(struct AtomStr *HB_i, struct AtomStr *HB_j, struct AtomStr *neighbor_i, struct AtomStr *neighbor_j);
void ThermostatEvent(struct AtomStr *target);
void PBCandCrossCellEvent(struct AtomStr *target, struct AtomListStr *thisList);
void CGSurfaceEvent(struct AtomStr *target, struct AtomListStr *thisThread);
void WallEvent(struct AtomStr *target);
void ObstEvent(struct AtomStr *target);
void TunnelEvent(struct AtomStr *target);
void ChargeEvent(struct AtomStr *target);
void SphObstEvent(struct AtomStr *target);
void InteractionEvent(struct AtomStr *target, struct AtomStr *partner, struct AtomStr *HBNeighbor_i, struct AtomStr *HBNeighbor_j, struct AtomListStr *thisList);
double HBEnergyChange(struct AtomStr *target, int skipAtom, char *type, struct AtomListStr *thisList);


void DoEvent(struct AtomListStr *thisList) {
    struct AtomStr *target, *partner;
    struct AtomStr *HBNeighbor_i, *HBNeighbor_j;
    int typenum;
    
    target  = thisList->target;
    partner = thisList->partner;
    typenum = target->dynamic->event.eventType;
    
    //switch to the event which will happen earliest
    switch (typenum) {
            
        case Coli_Event: //do collision event
            InteractionEvent(target, partner, NULL, NULL, thisList);
            collisioneventsum ++;
            break;
            
        case Bond_Event: //do bond event
            InteractionEvent(target, partner, NULL, NULL, thisList);
            bondeventsum ++;
            break;
            
        case HBNe_Event: //do HB neighbor event
            InteractionEvent(target, partner, NULL, NULL, thisList);
            HBNeighboreventsum ++;
            break;
            
        case HB_Event: //do HB event
            HBNeighbor_i = thisList->ptr[target->dynamic->HB.neighbor];
            HBNeighbor_j = thisList->ptr[partner->dynamic->HB.neighbor];
            InteractionEvent(target, partner, HBNeighbor_i, HBNeighbor_j, thisList);
            HBEvent(target, partner, HBNeighbor_i, HBNeighbor_j);
            HBeventsum ++;
            break;
            
        case Ther_Event: //do thermostat ghost collision
            ThermostatEvent(target);
            thermostateventsum ++;
            break;
            
        case PCC_Event: //do PCC event
            PBCandCrossCellEvent(target, thisList);
            pbcandcrosseventsum ++;
            break;
            
        case Lagv_Event: //do Langevin event
            break;
            
        case Wall_Event: //do wall collision
            WallEvent(target);
            walleventsum ++;
            break;
            
        case CGSu_Event:
            CGSurfaceEvent(target, thisList);
			CGeventsum ++;
            break;
            
        case Obst_Event: //do obstruction collision
            ObstEvent(target);
            walleventsum ++;
            break;
            
        case Tunl_Event: //do tunnel collision
            TunnelEvent(target);
            walleventsum ++;
            break;
            
        case Chrg_Event: //do charge interaction
            ChargeEvent(target);
            break;
            
        case SphO_Event: //do column interaction
            SphObstEvent(target);
            walleventsum ++;
            break;
            
        case Invd_Event:
        case Cacl_Event:
        case  TBD_Event:
            printf("!!ERROR!!: do_event switch number has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            
        default:
            break;
            
    }
    
    return;
}


void LinkList(char * type, struct AtomStr *target, int *oldIdx, struct AtomListStr *thisList) {
    int i, n;
    int cellindex;
    
    if (unlikely(strncmp(type, "full", 1) == 0 && target == NULL)) {
        
        //-------------------------------------------
        //assign a group of memory to store the link list
        INT_CALLOC(celllist, (atomnum + cellnum[1] * cellnum[2] * cellnum[3] + 1));
        
        //-------------------------------------------
        //calculate which subcells the atoms are in and generate the list
        for (i = 1; i <= atomnum; i++) {
            for (n = 1; n <= 3; n++) {
                
                thisList->ptr[i]->dynamic->cellIndex[n] =
                thisList->ptr[i]->dynamic->coordinate[n] / cellsize[n];
                //atom on the edge of box size has been removed by function pbc()
                //integer numbers
                //for a 3*3*3 system, cell index on one side would be:
                //  [0)    [1)   [2)
                //  [0)    [1)   [2)
                //  [0)    [1)   [2)
                
                //just in case the floating-point arithmetic error
                //the PBC-adjusted coordinates may be exactly equal to the upper edge of the PBC box
                if (thisList->ptr[i]->dynamic->cellIndex[n] == cellnum[n]) {
                    thisList->ptr[i]->dynamic->cellIndex[n] --;
                }
            }
            
            cellindex = thisList->ptr[i]->dynamic->cellIndex[3] * cellnum[1] * cellnum[2]
                      + thisList->ptr[i]->dynamic->cellIndex[2] * cellnum[1]
                      + thisList->ptr[i]->dynamic->cellIndex[1] + 1;
            
            celllist[i] = celllist[atomnum + cellindex];
            celllist[atomnum + cellindex] = i;
            
        }
        //-------------------------------------------
        
    } else if (likely(strncmp(type, "part", 1) == 0 && target != NULL)) {
        int newCellIndex;
        int oldCellIndex;
        int nextAtom;
        int *sourCellList = celllist;
        
        oldCellIndex = oldIdx[3] * cellnum[1] * cellnum[2]
                     + oldIdx[2] * cellnum[1]
                     + oldIdx[1] + 1;
        
        for (n = 1; n <= 3; n++) {
            target->dynamic->cellIndex[n] =
            target->dynamic->coordinate[n] / cellsize[n];
            
            if (unlikely(target->dynamic->cellIndex[n] == cellnum[n])) {
                target->dynamic->cellIndex[n]--;
            }
        }
        
        newCellIndex = target->dynamic->cellIndex[3] * cellnum[1] * cellnum[2]
                     + target->dynamic->cellIndex[2] * cellnum[1]
                     + target->dynamic->cellIndex[1] + 1;
        
        nextAtom = oldCellIndex + atomnum;
        if (unlikely(sourCellList[nextAtom] == 0)) {
            printf("!!ERROR!!: link list has errors!\n");
        }
        
        //remove the target atom at the old position
        while (sourCellList[nextAtom] != target->property->num) {
            nextAtom = sourCellList[nextAtom];
            
            if (unlikely(nextAtom == 0)) {
                printf("!!ERROR!!: link list has errors!\n");
            }
        }
        sourCellList[nextAtom] = sourCellList[target->property->num];
        
        //add the target atom into the new position
                                   nextAtom = newCellIndex + atomnum;
        sourCellList[target->property->num] = sourCellList[nextAtom];
        sourCellList[nextAtom]              = target->property->num;
        
    } else {
        printf("Link list arguments error!\n");
    }
}


void PBC(char * type, struct AtomStr *target, struct AtomListStr *thisList) {
    int dim = 3, step = 1;
    
    if (strcmp(wallExist, "no")) {
        if (strncmp(wallType, "sphere", 1) == 0) return;
        if (strncmp(wallType, "parallel", 1) == 0)
            step = 2;
        else if (strncmp(wallType, "cylinder", 1) == 0)
            dim = 1;
    }
    
    if (unlikely(strncmp(type, "full", 1) == 0 && target == NULL)) {
        
        for (int i = 1; i <= atomnum; i++) {
            for (int n = 1; n <= dim; n += step) {
                if (thisList->ptr[i]->dynamic->coordinate[n] >= boxDimension[n]) { //coordinate range [0,box size)
                    thisList->ptr[i]->dynamic->coordinate[n] -= boxDimension[n];
                } else if (thisList->ptr[i]->dynamic->coordinate[n] < 0) {
                    thisList->ptr[i]->dynamic->coordinate[n] += boxDimension[n];
                }
            }
        }
        
    } else if (likely(strncmp(type, "part", 1) == 0 && target != NULL)) {
        
        for (int n = 1; n <= dim; n += step) {            
            if (unlikely(target->dynamic->coordinate[n] >= boxDimension[n])) { //coordinate range [0,box size)
                target->dynamic->coordinate[n] -= boxDimension[n];
            } else if (unlikely(target->dynamic->coordinate[n] < 0)) {
                target->dynamic->coordinate[n] += boxDimension[n];
            }
        }
        
    } else {
        
        printf("PBC arguments error!\n");
        
    }
}


void PBCandCrossCellEvent(struct AtomStr *target, struct AtomListStr *thisList) {
    int n, oldIdx[4];
    
    oldIdx[1] = target->dynamic->cellIndex[1];
    oldIdx[2] = target->dynamic->cellIndex[2];
    oldIdx[3] = target->dynamic->cellIndex[3];
    
    for (n = 1; n <= 3; n++) {
        if (unlikely(target->dynamic->cellIndex[n] == 0 ||
            target->dynamic->cellIndex[n] == cellnum[n] - 1)) {
            PBC("part", target, thisList);
            break;
        }
    }
    
    LinkList("part", target, oldIdx, NULL);
    if (unlikely(oldIdx[1] == target->dynamic->cellIndex[1] &&
                 oldIdx[2] == target->dynamic->cellIndex[2] &&
                 oldIdx[3] == target->dynamic->cellIndex[3] &&
                 target->dynamic->event.time >= ZERO * 100)) {
        printf("check PCC calculation! atom: %i time: %lf\n", target->property->sequence.atomNum, target->dynamic->event.time);
    }
    
    return;
}


void HBEvent(struct AtomStr *HB_i, struct AtomStr *HB_j, struct AtomStr *neighbor_i, struct AtomStr *neighbor_j) {
    int temp;
    int atom_i = HB_i->property->num;
    int atom_j = HB_j->property->num;
    
    if (HB_i->dynamic->HB.interactionType == HBForm) {
        
        HB_i->dynamic->HB.bondConnection = atom_j;
        HB_j->dynamic->HB.bondConnection = atom_i;
        
        HB_i->property->typeofAtom = AtomTypeChange(HB_i->property->typeofAtom, 1);
        HB_j->property->typeofAtom = AtomTypeChange(HB_j->property->typeofAtom, 1);
        
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
        neighbor_i->dynamic->HBNeighbor.neighborPartner[temp] = atom_j;
    
        temp = ++neighbor_j->dynamic->HBNeighbor.neighborStatus;
        neighbor_j->dynamic->HBNeighbor.neighborPartner[temp] = atom_i;
        
        if (CheckAlphaHB(HB_i, HB_j))
            alphaHBformed ++;
        
        if (unlikely(connectionMap[atom_i][atom_j] & HB_CONNECT ||
                     connectionMap[atom_j][atom_i] & HB_CONNECT ||
                     connectionMap[atom_i][HB_j->dynamic->HB.neighbor] & NEIGHBOR_CONNECT ||
                     connectionMap[HB_j->dynamic->HB.neighbor][atom_i] & NEIGHBOR_CONNECT ||
                     connectionMap[atom_j][HB_i->dynamic->HB.neighbor] & NEIGHBOR_CONNECT ||
                     connectionMap[HB_i->dynamic->HB.neighbor][atom_j] & NEIGHBOR_CONNECT)) {
            printf("!!ERROR!!: HB connection map has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        //adjust the connection map
        connectionMap[atom_i][atom_j] ^= HB_CONNECT;
        connectionMap[atom_j][atom_i] ^= HB_CONNECT;
        connectionMap[atom_i][HB_j->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
        connectionMap[HB_j->dynamic->HB.neighbor][atom_i] ^= NEIGHBOR_CONNECT;
        connectionMap[atom_j][HB_i->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
        connectionMap[HB_i->dynamic->HB.neighbor][atom_j] ^= NEIGHBOR_CONNECT;
        
        HBnumformed ++;
        
    } else if (HB_i->dynamic->HB.interactionType == HBBreak) {
        
        HB_i->dynamic->HB.bondConnection = 0;
        HB_j->dynamic->HB.bondConnection = 0;
        
        HB_i->property->typeofAtom = AtomTypeChange(HB_i->property->typeofAtom, 0);
        HB_j->property->typeofAtom = AtomTypeChange(HB_j->property->typeofAtom, 0);
        
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
        ListRefresh(atom_j, neighbor_i->dynamic->HBNeighbor.neighborPartner, 1, 4);
        
        --neighbor_j->dynamic->HBNeighbor.neighborStatus;
        ListRefresh(atom_i, neighbor_j->dynamic->HBNeighbor.neighborPartner, 1, 4);
        
        if (CheckAlphaHB(HB_i, HB_j))
            alphaHBformed --;
        
        if (unlikely(!(connectionMap[atom_i][atom_j] & HB_CONNECT) ||
                     !(connectionMap[atom_j][atom_i] & HB_CONNECT) ||
                     !(connectionMap[atom_i][HB_j->dynamic->HB.neighbor] & NEIGHBOR_CONNECT) ||
                     !(connectionMap[HB_j->dynamic->HB.neighbor][atom_i] & NEIGHBOR_CONNECT) ||
                     !(connectionMap[atom_j][HB_i->dynamic->HB.neighbor] & NEIGHBOR_CONNECT) ||
                     !(connectionMap[HB_i->dynamic->HB.neighbor][atom_j] & NEIGHBOR_CONNECT))) {
            printf("!!ERROR!!: HB connection map has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        //adjuct the connection map
        connectionMap[atom_i][atom_j] ^= HB_CONNECT;
        connectionMap[atom_j][atom_i] ^= HB_CONNECT;
        connectionMap[atom_i][HB_j->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
        connectionMap[HB_j->dynamic->HB.neighbor][atom_i] ^= NEIGHBOR_CONNECT;
        connectionMap[atom_j][HB_i->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
        connectionMap[HB_i->dynamic->HB.neighbor][atom_j] ^= NEIGHBOR_CONNECT;
        
        HBnumformed --;
        
        if (unlikely(HBnumformed < alphaHBformed)) {
            printf("!WARNING!: total HB number is smaller than the alpha HB number! %s:%i\n", __FILE__, __LINE__);
            warningsum ++;
        }
    }
    
    return;
}


void ThermostatEvent(struct AtomStr *target) {
    double speed[4];
    double netV[4] = {0};

    if (chargeAA[COM] && flow.mark == 3) { //!!not parallizable!! will use the raw data
        CalCOMV(target->property->sequence.proteinNum, netV);
    }
    
    if (strncmp(thermostatType, "Andersen", 1) == 0) {
        for (int n = 1; n <= 3; n++) {
        regenerate:
            speed[n] = RandomVelocity(target->property->mass);
            speed[n] += flow.constV.v[n] + netV[n] +
                        flow.force.a[target->property->num][n] *
                        (currenttime - flow.force.timeRec[target->property->num]);
            
            if (unlikely(speed[n] == 0)) {
                goto regenerate;
            }
        }
        
        TRANSFER_VECTOR(target->dynamic->velocity, speed);
        flow.force.timeRec[target->property->num] = currenttime;
    }
}


void WallEvent(struct AtomStr *target) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *targetSpeed = target->dynamic->velocity;
    double speed_i[4] = {0}, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisWall = NULL;
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(speed_i, target->dynamic->velocity);
    
    if (strcmp(wallExist, "smooth") == 0) {
        thisWall = &wall[0];
        TRANSFER_VECTOR(position_j, thisWall->dynamic->coordinate);
        TRANSFER_VECTOR(speed_j, thisWall->dynamic->velocity);
        
        if (strncmp(wallType, "parallel", 1) == 0) {
            position_i[1] = position_i[3] = 0;
               speed_i[1] =    speed_i[3] = 0;
        } else if (strncmp(wallType, "cylinder", 1) == 0) {
            position_i[1] = 0;
               speed_i[1] = 0;
        }
    }
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass;
    potential = target->dynamic->event.potential;
    
    switch (target->dynamic->event.subEventType) {
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
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / target->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(targetSpeed, term, target->dynamic->velocity);

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


void ObstEvent(struct AtomStr *target) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *targetSpeed = target->dynamic->velocity;
    double speed_i[4] = {0}, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisWall = &obstObj.obst[target->dynamic->event.partner / -10];
    
    int dim = (target->dynamic->event.partner * -1) % 10;
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(speed_i,    target->dynamic->velocity);

    TRANSFER_VECTOR(position_j,   thisWall->dynamic->coordinate);
    TRANSFER_VECTOR(speed_j,      thisWall->dynamic->velocity);
    
    for (int i = 1; i <= 3; i ++) {
        if (dim != i) {
            position_i[i] = 0; speed_i[i] = 0;
            position_j[i] = 0; speed_j[i] = 0;
        }
    }
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(   speed_i,    speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij  = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass;
    potential = target->dynamic->event.potential;
    
    switch (target->dynamic->event.subEventType) {
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
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / target->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(targetSpeed, term, target->dynamic->velocity);
    
    return;
}


void TunnelEvent(struct AtomStr *target) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *targetSpeed = target->dynamic->velocity;
    double speed_i[4] = {0}, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisWall = &tunlObj.tunnel;
    
    int num = target->dynamic->event.partner / -10000;
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, tunlObj.position[num]);
    TRANSFER_VECTOR(speed_i, target->dynamic->velocity);
    TRANSFER_VECTOR(speed_j,   thisWall->dynamic->velocity);

    position_i[1] = position_j[1]; //end position
       speed_i[1] =    speed_j[1] = 0;

    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass;
    potential = target->dynamic->event.potential;
    
    switch (target->dynamic->event.subEventType) {
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
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / target->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(targetSpeed, term, target->dynamic->velocity);
    
    return;
}


void ChargeEvent(struct AtomStr *target) {
    int s_vel = 0;
    int chargeNum = target->dynamic->event.partner / -100;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, *speed_j = NULL;
    double position_i[4] = {0}, position_j[4] = {0};
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, flow.charge.position[chargeNum]);
    speed_i = target->dynamic->velocity;
    speed_j = flow.charge.velocity[chargeNum];
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass;
    potential = target->dynamic->event.potential;
    
    switch (target->dynamic->event.subEventType) {
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
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / target->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, target->dynamic->velocity);
    
    return;
}


void SphObstEvent(struct AtomStr *target) {
    int s_vel = 0;
    int thisObst = target->dynamic->event.partner / -1000;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass;
    double *speed_i = NULL, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, SphObstObj.position[thisObst]);
    
    speed_i = target->dynamic->velocity;
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(   speed_i,    speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij  = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass;
    
    switch (target->dynamic->event.subEventType) {
        case CoreColli:
            s_vel = 1;
            potential = 0;
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / target->property->mass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(speed_i, term, target->dynamic->velocity);
    
    return;
}


void InteractionEvent (struct AtomStr *target, struct AtomStr *partner, struct AtomStr *HBNeighbor_i, struct AtomStr *HBNeighbor_j, struct AtomListStr *thisList)
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
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, partner->dynamic->coordinate);
    speed_i = target->dynamic->velocity;
    speed_j = partner->dynamic->velocity;
    
    for (int n = 1; n <= 3; n++) {
        if (target->dynamic->cellIndex[n] != partner->dynamic->cellIndex[n]) {
            if (target->dynamic->cellIndex[n] == 0 && partner->dynamic->cellIndex[n] == cellnum[n] - 1) {
                position_j[n] -= boxDimension[n];
            } else if (target->dynamic->cellIndex[n] == cellnum[n] - 1 && partner->dynamic->cellIndex[n] == 0) {
                position_j[n] += boxDimension[n];
            }
        }
    } //apply pbc condition
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass * partner->property->mass / (target->property->mass + partner->property->mass);
    potential = target->dynamic->event.potential;
    
    if (unlikely(target->dynamic->event.eventType    == HB_Event && // HB event only
                (target->dynamic->event.subEventType == TBDHBForm || target->dynamic->event.subEventType == TBDHBBrik))) { //may form or break HB
        
        double position_neighbor_i[4];
        double position_neighbor_j[4];
        double r_c12_neighbor_i[4], r_neighbor_i2;
        double r_c12_neighbor_j[4], r_neighbor_j2;
        double accumPoten[6], tmp = 0;
        
        //the HB neighbors cannot be occupied before the HB forms
        if ((connectionMap[target->property->num][HBNeighbor_j->property->num] & NEIGHBOR_CONNECT ||
             connectionMap[   partner->property->num][HBNeighbor_i->property->num] & NEIGHBOR_CONNECT) &&
            target->dynamic->event.subEventType == TBDHBForm) {
            target->dynamic->HB.interactionType = HBSkip;
            return;
        }
        
        //-------------------------
        //adjust neighbor PBC
        TRANSFER_VECTOR(position_neighbor_i, HBNeighbor_i->dynamic->coordinate);
        TRANSFER_VECTOR(position_neighbor_j, HBNeighbor_j->dynamic->coordinate);
        
        for (int n = 1; n <= 3; n++) {
            
            //------------------
            //for HB_i
            if (target->dynamic->cellIndex[n] != HBNeighbor_i->dynamic->cellIndex[n]) {
                if (target->dynamic->cellIndex[n] == 0 && HBNeighbor_i->dynamic->cellIndex[n] == cellnum[n] - 1) {
                    position_neighbor_i[n] -= boxDimension[n];
                } else if (target->dynamic->cellIndex[n] == cellnum[n] - 1 && HBNeighbor_i->dynamic->cellIndex[n] == 0) {
                    position_neighbor_i[n] += boxDimension[n];
                }
            }
            //------------------
            
            
            //------------------
            //for HB_j
            //all based on the coordinates of HB_i
            if (target->dynamic->cellIndex[n] != HBNeighbor_j->dynamic->cellIndex[n]) {
                if (target->dynamic->cellIndex[n] == 0 && HBNeighbor_j->dynamic->cellIndex[n] == cellnum[n] - 1) {
                    position_neighbor_j[n] -= boxDimension[n];
                } else if (target->dynamic->cellIndex[n] == cellnum[n] - 1 && HBNeighbor_j->dynamic->cellIndex[n] == 0) {
                    position_neighbor_j[n] += boxDimension[n];
                }
            }
            //------------------
        }
        
        //measure the distances from the neighbors and calculate the total energy change
        DOT_MINUS(position_i, position_neighbor_j, r_c12_neighbor_i);
        DOT_MINUS(position_j, position_neighbor_i, r_c12_neighbor_j);
        r_neighbor_i2 = DOT_PROD(r_c12_neighbor_i, r_c12_neighbor_i);   r_neighbor_j2 = DOT_PROD(r_c12_neighbor_j, r_c12_neighbor_j);
        
        accumPoten[0] = HBEnergyChange(target,    partner->property->num, "_nonHB", thisList);
        accumPoten[1] = HBEnergyChange(partner,    target->property->num, "_nonHB", thisList);
        
        accumPoten[2] = HBEnergyChange(target,    partner->property->num, "_HB",    thisList);
        accumPoten[3] = HBEnergyChange(partner,    target->property->num, "_HB",    thisList);

        FindPair(HBNeighbor_j, target,  "neighbor", tmp, r_neighbor_i2, &temp, &temp, &neighborLowerPotential, &neighborUpperPotential, &accumPoten[4], partner, 0);
        FindPair(HBNeighbor_i, partner, "neighbor", tmp, r_neighbor_j2, &temp, &temp, &neighborLowerPotential, &neighborUpperPotential, &accumPoten[5], target, 0);
        
        if (target->dynamic->event.subEventType == TBDHBForm) { //form
			potential -= HBBarrier(target, partner); //HB barrier
            potential += (accumPoten[2] - accumPoten[0] + accumPoten[3] - accumPoten[1] + accumPoten[4] + accumPoten[5]);
        } else if (target->dynamic->event.subEventType == TBDHBBrik) { //break
            potential += HBBarrier(target, partner); //HB barrier
            potential += (accumPoten[0] - accumPoten[2] + accumPoten[1] - accumPoten[3] - accumPoten[4] - accumPoten[5]);
        } else {
            printf("!ERROR!: HB event type has something wrong!\n");
            exit(EXIT_FAILURE);
        }
        
        if (potential < -0.5 * INFINIT && target->dynamic->event.subEventType == TBDHBBrik) { //means this neighbortime must be shorter, i.e. there will be a neighbor event before the current HB event
            target->dynamic->HB.interactionType = HBSkip;
            return;
        } else if (potential > 0.5 * INFINIT && target->dynamic->event.subEventType == TBDHBForm) {
            //try to form HB, but there are neighbor atoms that not at the right positions
            target->dynamic->HB.interactionType = HBSkip;
            return;
        }
        
        if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass > 0) {
            
            if (target->dynamic->event.subEventType == TBDHBForm) {
                target->dynamic->event.subEventType = WellCaptr;
                target->dynamic->HB.interactionType = HBForm;
            } else if (target->dynamic->event.subEventType == TBDHBBrik) {
                target->dynamic->event.subEventType = WellEscap;
                target->dynamic->HB.interactionType = HBBreak;
            } else {
                printf("!!ERROR!!: should not reach here! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
        } else {
            
            //not enough:
            if (target->dynamic->event.subEventType == TBDHBForm) {
                target->dynamic->HB.interactionType = HBSkip;
                return;
            } else if (target->dynamic->event.subEventType == TBDHBBrik) {
                potential = 0;
                target->dynamic->event.subEventType = WellBounc;
                target->dynamic->HB.interactionType = NoEvent;
            } else {
                printf("!!ERROR!!: should not reach here! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    switch (target->dynamic->event.subEventType) {
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
    ((target->property->mass + partner->property->mass) * d_ij2);
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (strcmp(target->property->extraProperty[1], "_CONS")) {
        FACTOR_PROD((partner->property->mass * phi), r_ij, term);
        DOT_PLUS(speed_i, term, target->dynamic->velocity);
    }
    
    if (strcmp(partner->property->extraProperty[1], "_CONS")) {
        FACTOR_PROD((target->property->mass * phi), r_ij, term);
        DOT_MINUS(speed_j, term, partner->dynamic->velocity);
    }

    return;
}


double HBEnergyChange(struct AtomStr *target, int skipAtom, char *type, struct AtomListStr *thisList) {
    int atom_i, atom_j;
    int cell_neighbor[4], cellIndex;
    int *sourCellList = celllist;
    double cutoffRange = cutoffr * cutoffr;
    double energyChange = 0;
    double position_i[4], position_j[4];
    double r_temp[4], positionshift[4] = {0};
    double r_ij[4], r_2;
    double tmp = 0, accumPoten = 0;
    struct AtomStr *partner;
    
    atom_i = target->property->num;
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    
    for (int i = 0; i < 27; ++i) {
        //scan the neighborhood 27 subcells, include the target subcell itself
        cell_neighbor[3] = target->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = target->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = target->dynamic->cellIndex[1] + thisList->x[i];
        
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
        
        atom_j = sourCellList[atomnum + cellIndex];
        while (atom_j) {
            partner = thisList->ptr[atom_j];
            
            if (atom_i != atom_j &&
                !(connectionMap[atom_i][atom_j] & BOND_CONNECT) &&
                !(connectionMap[atom_i][atom_j] & CONSTRAINT_CONNECT) &&
                atom_j != skipAtom) {
                
                TRANSFER_VECTOR(position_j, partner->dynamic->coordinate);
                PBCShift(target, partner, positionshift);
                
                DOT_PLUS(position_j, positionshift, r_temp);
                DOT_MINUS(position_i, r_temp, r_ij); //here r_ij is before the collision
                r_2 = DOT_PROD(r_ij, r_ij);
                
                if (r_2 <= cutoffRange) {
                    if (strcmp(type, "_nonHB") == 0) {
                        FindPair(target, partner, "collision", tmp, r_2, &tmp, &tmp, &tmp, &tmp, &accumPoten, NULL, -1);
                    } else if (strcmp(type, "_HB") == 0) {
                        if (atom_j != thisList->ptr[skipAtom]->dynamic->HB.neighbor) {
                            FindPair(target, partner, "collision", tmp, r_2, &tmp, &tmp, &tmp, &tmp, &accumPoten, NULL, 1);
                        }
                    } else {
                        printf("!!ERROR!!: type is invalid! %s:%i\n", __FILE__, __LINE__);
                    }
                    if (accumPoten == INFINIT) accumPoten = 0;
                    energyChange += accumPoten;
                }
            }
            atom_j = sourCellList[atom_j];
        }        
    }
    
    return energyChange;
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

int CheckAlphaHB(struct AtomStr* HB_i, struct AtomStr* HB_j) {
    if (HB_i->property->sequence.proteinNum == HB_j->property->sequence.proteinNum &&
        ((strcmp(HB_i->property->name, "H") == 0 && strcmp(HB_j->property->name, "O") == 0 && HB_j->property->sequence.aminoacidNum == HB_i->property->sequence.aminoacidNum - 4) ||
         (strcmp(HB_i->property->name, "O") == 0 && strcmp(HB_j->property->name, "H") == 0 && HB_i->property->sequence.aminoacidNum == HB_j->property->sequence.aminoacidNum - 4)))
        return TRUE;
    
    return FALSE;
}
