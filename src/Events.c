//
//  Events.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

void CGSurfaceEvent(struct AtomStr *target, struct AtomListStr *thisThread);

static int CheckAlphaHB(struct AtomStr* HB_i, struct AtomStr* HB_j);
static int CheckHB_NeighborConnection(int atom_i, int atom_j, int *neighborList_i, int *neighborList_j, int flag);
static void HBEvent(struct AtomStr *HB_i, struct AtomStr *HB_j, int *neighborList_i, int *neighborList_j);
static void ThermostatEvent(struct AtomStr *target);
static void PBCandCrossCellEvent(struct AtomStr *target, struct AtomListStr *thisList);
static void WallEvent(struct AtomStr *target);
static void ObstEvent(struct AtomStr *target);
static void TunnelEvent(struct AtomStr *target);
static void ChargeEvent(struct AtomStr *target);
static void SphObstEvent(struct AtomStr *target);
static void RemoveNeighborFromHB(int *currentStatus, int *list, int *neighborList);
static void RemoveHBFromNeighbor(int *list, int HBNum);
static void InteractionEvent(struct AtomStr *target, struct AtomStr *partner, int *neighborList_i, int *neighborList_j, struct AtomListStr *thisList);
static void ShowHideAtoms(int *atomList, int flag);
static void HBEnergyChangeCombined(struct AtomStr *target, struct AtomStr *partner, double *accumPoten, struct AtomListStr *thisList, int *excl_ap2, int *excl_ap3, int keepType_i, int keepType_j);


void DoEvent(struct AtomListStr *thisList) {
    struct AtomStr *target, *partner;
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
            InteractionEvent(target,
                             partner,
                             target->dynamic->HB.neighbor,
                             partner->dynamic->HB.neighbor,
                             thisList);
            HBEvent(target, partner, target->dynamic->HB.neighbor, partner->dynamic->HB.neighbor);
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
        if (!celllist)
            celllist = (int *)calloc(atomnum + cellnum[1] * cellnum[2] * cellnum[3] + 1, sizeof(int));
        
        //-------------------------------------------
        //calculate which subcells the atoms are in and generate the list
        for (i = 1; i <= atomnum; i++) {
            for (n = 1; n <= 3; n++) {
                
                thisList->ptr[i]->dynamic->cellIndex[n] = floor(
                thisList->ptr[i]->dynamic->coordinate[n] / cellsize[n]);
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
            target->dynamic->cellIndex[n] = floor(
            target->dynamic->coordinate[n] / cellsize[n]);
            
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
    
    if (strcmp(surfInfo.style, "implicit") == 0) {
        if (strncmp(surfInfo.type, "sphere", 1) == 0) return;
        if (strncmp(surfInfo.type, "parallel", 1) == 0)
            step = 2;
        else if (strncmp(surfInfo.type, "cylinder", 1) == 0)
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


void HBEvent(struct AtomStr *HB_i, struct AtomStr *HB_j, int *neighborList_i, int *neighborList_j) {
    int atom_i = HB_i->property->num;
    int atom_j = HB_j->property->num;
    
    if (HB_i->dynamic->HB.interactionType == HBForm) {
        
        HBAddPartner(HB_i, atom_j);
        HBAddPartner(HB_j, atom_i);
        
        HB_i->property->typeofAtom[0] = 2; //change the type to _HB
        HB_j->property->typeofAtom[0] = 2;
        
#ifdef VIS
        if (visual) {
            ChangeColor(HB_i->property->type, HB_i->property->color);
            ChangeColor(HB_j->property->type, HB_j->property->color);
        }
#endif
        
        //===================================
        //for HB acceptor and donor
        AddNeighbor2HB(&HB_i->dynamic->HBNeighbor.neighborStatus,
                        HB_i->dynamic->HBNeighbor.neighborPartner,
                        HB_j->dynamic->HB.neighbor);
        AddNeighbor2HB(&HB_j->dynamic->HBNeighbor.neighborStatus,
                        HB_j->dynamic->HBNeighbor.neighborPartner,
                        HB_i->dynamic->HB.neighbor);
        
        //for neighbors
        AddHB2Neighbor(neighborList_i, atom_j);
        AddHB2Neighbor(neighborList_j, atom_i);
        
        if (CheckAlphaHB(HB_i, HB_j))
            alphaHBformed ++;

        if (CheckHB_NeighborConnection(atom_i, atom_j, neighborList_i, neighborList_j, 0)) {
            printf("!!ERROR!!: HB connection map has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        //adjust the connection map
        SwitchConnection(atom_i, atom_j, neighborList_i, neighborList_j);

        //!!for crosslinking simulation, need further optimization!!
        //hide the listed atoms
        //ShowHideAtoms(HB_i->dynamic->HBNeighbor.removedAtom, 0);
        //ShowHideAtoms(HB_j->dynamic->HBNeighbor.removedAtom, 0);

        HBnumformed ++;
        
    } else if (HB_i->dynamic->HB.interactionType == HBBreak) {
        
        HBRemovePartner(HB_i, atom_j);
        HBRemovePartner(HB_j, atom_i);
        
        //change the type back to original once the atom has no HB left
        if (HBPartnerNum(HB_i) == 0) HB_i->property->typeofAtom[0] = 1;
        if (HBPartnerNum(HB_j) == 0) HB_j->property->typeofAtom[0] = 1;
        
#ifdef VIS
        if (visual) {
            ChangeColor(HB_i->property->type, HB_i->property->color);
            ChangeColor(HB_j->property->type, HB_j->property->color);
        }
#endif
        
        //===================================
        //for HB acceptor and donor
        RemoveNeighborFromHB(&HB_i->dynamic->HBNeighbor.neighborStatus,
                              HB_i->dynamic->HBNeighbor.neighborPartner,
                              neighborList_j);
        RemoveNeighborFromHB(&HB_j->dynamic->HBNeighbor.neighborStatus,
                              HB_j->dynamic->HBNeighbor.neighborPartner,
                              neighborList_i);
        
        //for neighbor
        RemoveHBFromNeighbor(neighborList_i, atom_j);
        RemoveHBFromNeighbor(neighborList_j, atom_i);
        
        if (CheckAlphaHB(HB_i, HB_j))
            alphaHBformed --;
        
        if (CheckHB_NeighborConnection(atom_i, atom_j, neighborList_i, neighborList_j, 1)) {
            printf("!!ERROR!!: HB connection map has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        //adjuct the connection map
        SwitchConnection(atom_i, atom_j, neighborList_i, neighborList_j);

        //!!for crosslinking simulation, need further optimization!!
        //unhide the listed atoms
        //ShowHideAtoms(HB_i->dynamic->HBNeighbor.removedAtom, 1);
        //ShowHideAtoms(HB_j->dynamic->HBNeighbor.removedAtom, 1);
        
        HBnumformed --;
        
        if (unlikely(HBnumformed < alphaHBformed)) {
            printf("!WARNING![%li]: total HB number is smaller than the alpha HB number! %s:%i\n", ++warningsum, __FILE__, __LINE__);
        }
    }
    
    return;
}


void ThermostatEvent(struct AtomStr *target) {
    double speed[4] = {0};
    double netV[4] = {0};

    if (chargeRes[COM] && flow.mark == 3) { //!!not parallizable!! will use the raw data
        CalCOMV(target->property->sequence.molNum, netV);
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
    
    return;
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
    
    if (strcmp(surfInfo.style, "implicit") == 0) {
        thisWall = &wall[0];
        TRANSFER_VECTOR(position_j, thisWall->dynamic->coordinate);
        TRANSFER_VECTOR(speed_j, thisWall->dynamic->velocity);
        
        if (strncmp(surfInfo.type, "parallel", 1) == 0) {
            position_i[1] = position_i[3] = 0;
               speed_i[1] =    speed_i[3] = 0;
        } else if (strncmp(surfInfo.type, "cylinder", 1) == 0) {
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
    double realTime = REAL_TIME(currenttime); // wallDyn.curTime is in real time
    if (realTime >= wallDyn.curTime) {
        if (strcmp(wallDyn.mark, "pulse") == 0) {
            if ((wallDyn.origBoxDim[2] - boxDimension[2] >= wallDyn.size) ||
                (wallDyn.origBoxDim[2] == boxDimension[2]) ||
                (fabs(boxDimension[2] + (int)(realTime / wallDyn.rate + 1) * wallDyn.size - (realTime - 1) * wallDyn.step - wallDyn.origBoxDim[2]) <= wallDyn.step && 
                wallDyn.touch >= 5) //balance the compressing and expanding time
                ) {
                wallDyn.sign *= -1;
                wallDyn.touch = 0;
            }
        }
        
        wallDyn.curTime += 1;
        
        if (strcmp(surfInfo.type, "parallel") == 0) {
            boxDimension[2] += wallDyn.sign * wallDyn.step;
        } else if (strcmp(surfInfo.type, "cylinder") == 0) {
            boxDimension[2] += wallDyn.sign * wallDyn.step;
            boxDimension[3] += wallDyn.sign * wallDyn.step;
        } else if (strcmp(surfInfo.type, "sphere") == 0) {
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
    double phi, potential;
    double *speed_i = NULL, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, SphObstObj.position[thisObst]);
    
    speed_i = target->dynamic->velocity;
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(   speed_i,    speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij  = DOT_PROD(r_ij, v_ij);
    
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


// A failed HB formation at a distance where the pair's nonbonded potential also
// has a step is processed as that nonbonded collision (enter the step or bounce).
static int FailedHBToCollision(struct AtomStr *target, struct AtomStr *partner, double d_ij2, double *potential) {
    double lowerLimit, upperLimit, lowerPotential, upperPotential, accum;
    // positive direction = separating: the pair is classified on the outer side of the step
    FindPair(target, partner, InterColli, 1.0, d_ij2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &accum, NULL, 0);
    if (fabs(lowerLimit - d_ij2) > 1e-6 * d_ij2) return FALSE;
    *potential = lowerPotential;
    target->dynamic->event.subEventType = TBDCaptr;
    target->dynamic->HB.interactionType = NoEvent;
    return TRUE;
}

void InteractionEvent (struct AtomStr *target, struct AtomStr *partner, int *neighborList_i, int *neighborList_j, struct AtomListStr *thisList)
{
    int s_vel = 0;
    double v_ij[4];
    double b_ij, r_ij[4];
    double term[4];
    double phi, potential, reducedMass;
    double d_ij2;
    double *speed_i, *speed_j;
    double position_i[4], position_j[4];
    double lowerLimit, upperLimit, lowerPotential, upperPotential, accumPoten_val;
    
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_j, partner->dynamic->coordinate);
    speed_i = target->dynamic->velocity;
    speed_j = partner->dynamic->velocity;
    
    double pbcShift[4];
    PBCShift(target, partner, pbcShift);
    for (int n = 1; n <= 3; n++) {
        position_j[n] += pbcShift[n];
    } //apply pbc condition
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(speed_i, speed_j, v_ij);
    d_ij2 = DOT_PROD(r_ij, r_ij);
    b_ij = DOT_PROD(r_ij, v_ij);
    reducedMass = target->property->mass * partner->property->mass / (target->property->mass + partner->property->mass);
    potential = target->dynamic->event.potential;
    
    if (unlikely(target->dynamic->event.eventType    == HB_Event && // HB event only
                (target->dynamic->event.subEventType == TBDHBForm || target->dynamic->event.subEventType == TBDHBBrik))) { //may form or break HB
        
        struct AtomStr **HBNeighbor_i, **HBNeighbor_j;
        double **position_neighbor_i, **position_neighbor_j;
        double r_c12_neighbor_i[4], r_neighbor_i2;
        double r_c12_neighbor_j[4], r_neighbor_j2;
        double accumPoten[6] = {0}, potRecord = 0, tmp = 0;
        
        HBNeighbor_i = (struct AtomStr **)calloc((neighborList_i[0] + 1), sizeof(struct AtomStr *));
        HBNeighbor_j = (struct AtomStr **)calloc((neighborList_j[0] + 1), sizeof(struct AtomStr *));
        _CALLOC2(double, position_neighbor_i, (neighborList_i[0] + 1), 4);
        _CALLOC2(double, position_neighbor_j, (neighborList_j[0] + 1), 4);
        
        //the HB neighbors cannot be occupied before the HB forms
        if (CheckHB_NeighborConnection( target->property->num,
                                       partner->property->num,
                                       neighborList_i,
                                       neighborList_j,
                                       0) &&
            target->dynamic->event.subEventType == TBDHBForm) {
            free(HBNeighbor_i);
            free(HBNeighbor_j);
            _FREE2(position_neighbor_i);
            _FREE2(position_neighbor_j);
            if (FailedHBToCollision(target, partner, d_ij2, &potential)) goto do_collision;
            target->dynamic->HB.interactionType = HBSkip;
            return;
        }
        
        //-------------------------
        //adjust neighbor PBC
        NEIGHBORSCAN(neighborList_i, i) {
            HBNeighbor_i[i] = thisList->ptr[neighborList_i[i]];
            UpdateAtom(HBNeighbor_i[i], currenttime);
            TRANSFER_VECTOR(position_neighbor_i[i],
                            HBNeighbor_i[i]->dynamic->coordinate);
        }
        NEIGHBORSCAN(neighborList_j, i) {
            HBNeighbor_j[i] = thisList->ptr[neighborList_j[i]];
            UpdateAtom(HBNeighbor_j[i], currenttime);
            TRANSFER_VECTOR(position_neighbor_j[i],
                            HBNeighbor_j[i]->dynamic->coordinate);
        }

        
        //all neighbor positions are based on the coordinates of HB_i
        NEIGHBORSCAN(neighborList_i, i) {
            PBCShift(target, HBNeighbor_i[i], pbcShift);
            for (int n = 1; n <= 3; n++) position_neighbor_i[i][n] += pbcShift[n];
        }
        NEIGHBORSCAN(neighborList_j, i) {
            PBCShift(target, HBNeighbor_j[i], pbcShift);
            for (int n = 1; n <= 3; n++) position_neighbor_j[i][n] += pbcShift[n];
        }
        
        //measure the distances from the neighbors and calculate the total energy change
        //for target calculation
        NEIGHBORSCAN(neighborList_j, n) {
            DOT_MINUS(position_i, position_neighbor_j[n], r_c12_neighbor_i);
            r_neighbor_i2 = DOT_PROD(r_c12_neighbor_i, r_c12_neighbor_i);
            
            FindPair(HBNeighbor_j[n], target, InterNeighbor, tmp, r_neighbor_i2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &potRecord, partner, 0);
            accumPoten[4] += potRecord;
        }
        
        //for partner calculation
        NEIGHBORSCAN(neighborList_i, n) {
            DOT_MINUS(position_j, position_neighbor_i[n], r_c12_neighbor_j);
            r_neighbor_j2 = DOT_PROD(r_c12_neighbor_j, r_c12_neighbor_j);
            
            FindPair(HBNeighbor_i[n], partner, InterNeighbor, tmp, r_neighbor_j2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &potRecord, target, 0);
            accumPoten[5] += potRecord;
        }
        
        // neighborList_j: C (partner/O's bonded neighbors, not bonded to target/H) → excluded from H's _HB scan
        // neighborList_i: N (target/H's bonded neighbors, not bonded to partner/O) → excluded from O's _HB scan
        // An atom that keeps another HB (multi-partner acceptor) keeps its _HB type,
        // so its type-change term is zero: evaluate both sums with the _HB type.
        int forming = (target->dynamic->event.subEventType == TBDHBForm);
        int keepType_i = HBPartnerNum(target)  > (forming ? 0 : 1);
        int keepType_j = HBPartnerNum(partner) > (forming ? 0 : 1);
        HBEnergyChangeCombined(target, partner, accumPoten, thisList, neighborList_j, neighborList_i, keepType_i, keepType_j);


        if (target->dynamic->event.subEventType == TBDHBForm) { //form
			potential -= HBBarrier(target, partner); //HB barrier
            potential += (accumPoten[2] - accumPoten[0] + accumPoten[3] - accumPoten[1] + accumPoten[4] + accumPoten[5]);
        } else if (target->dynamic->event.subEventType == TBDHBBrik) { //break
            potential += HBBarrier(target, partner); //HB barrier
            potential += (accumPoten[0] - accumPoten[2] + accumPoten[1] - accumPoten[3] - accumPoten[4] - accumPoten[5]);
        } else {
            printf("!!ERROR!!: HB event type has something wrong! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

#ifdef HB_DEBUG
        {
            const char *evtName = (target->dynamic->event.subEventType == TBDHBForm) ? "FORM" :
                                  (target->dynamic->event.subEventType == TBDHBBrik) ? "BRIK" : "????";
            double disc = b_ij * b_ij - 2 * d_ij2 * potential / reducedMass;
            printf("[HB_DEBUG] t=%.4f %s H#%i(res%i) O#%i(res%i) | "
                   "pot=%.4f b2=%.4f disc=%.4f d2=%.4f m=%.4f | "
                   "ap[0-5]=%.4f %.4f %.4f %.4f %.4f %.4f | "
                   "%s\n",
                   currenttime, evtName,
                   target->property->num,  target->property->sequence.resNum,
                   partner->property->num, partner->property->sequence.resNum,
                   potential,
                   b_ij * b_ij,
                   disc,
                   d_ij2, reducedMass,
                   accumPoten[0], accumPoten[1], accumPoten[2], accumPoten[3], accumPoten[4], accumPoten[5],
                   (disc > 0) ? ((target->dynamic->event.subEventType == TBDHBForm) ? "->FORM" : "->BREAK") : "->BOUNCE/SKIP");
            fflush(stdout);
        }
#endif

        if (potential < -0.5 * INFINIT && target->dynamic->event.subEventType == TBDHBBrik) {
            printf("!WARNING![%li]: weird thing happens while calculating HB energy change, please check the situations!\n", ++warningsum);
            printf("         %s:%i\n", __FILE__, __LINE__);
            
            target->dynamic->HB.interactionType = HBSkip;

            free(HBNeighbor_i);
            free(HBNeighbor_j);
            _FREE2(position_neighbor_i);
            _FREE2(position_neighbor_j);
            return;
        } else if (potential > 0.5 * INFINIT && target->dynamic->event.subEventType == TBDHBForm) {
            //try to form HB, but there are neighbor atoms that not at the right positions
            free(HBNeighbor_i);
            free(HBNeighbor_j);
            _FREE2(position_neighbor_i);
            _FREE2(position_neighbor_j);
            if (FailedHBToCollision(target, partner, d_ij2, &potential)) goto do_collision;
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
                free(HBNeighbor_i);
                free(HBNeighbor_j);
                _FREE2(position_neighbor_i);
                _FREE2(position_neighbor_j);
                if (FailedHBToCollision(target, partner, d_ij2, &potential)) goto do_collision;
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
        
        free(HBNeighbor_i);
        free(HBNeighbor_j);
        _FREE2(position_neighbor_i);
        _FREE2(position_neighbor_j);
    }
    
do_collision:
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
    
    if (!FindStrInList(target->property->extraProperty, "_CONS", 0)) {
        FACTOR_PROD((partner->property->mass * phi), r_ij, term);
        DOT_PLUS(speed_i, term, target->dynamic->velocity);
    }
    
    if (!FindStrInList(partner->property->extraProperty, "_CONS", 0)) {
        FACTOR_PROD((target->property->mass * phi), r_ij, term);
        DOT_MINUS(speed_j, term, partner->dynamic->velocity);
    }

    return;
}


double HBEnergyChange(struct AtomStr *target, int skipAtom, char *type, struct AtomListStr *thisList) {
    int atom_i, atom_j;
    int cell_neighbor[4], cellIndex;
    int *sourCellList = celllist;
    int checkStatus = 0;
    double cutoffRange = cutoffr * cutoffr;
    double energyChange = 0;
    double position_i[4], position_j[4];
    double r_temp[4], positionshift[4] = {0};
    double r_ij[4], r_2;
    double tmp = 0, accumPoten = 0;
    double lowerLimit, upperLimit, lowerPotential, upperPotential;
    struct AtomStr *partner;
    
    atom_i = target->property->num;
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    
    for (int i = 0; i < thisList->stencilNum; ++i) {
        //scan the neighborhood 27 subcells, include the target subcell itself
        cell_neighbor[3] = target->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = target->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = target->dynamic->cellIndex[1] + thisList->x[i];
        
        for (int j = 1; j <= 3; j++) {
            if (unlikely(cell_neighbor[j] < 0)) {
                cell_neighbor[j] += cellnum[j];
            } else if (unlikely(cell_neighbor[j] >= cellnum[j])) {
                cell_neighbor[j] -= cellnum[j];
            }
        }
        
        cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2]
        + cell_neighbor[2] * cellnum[1]
        + cell_neighbor[1] + 1;
        
        atom_j = sourCellList[atomnum + cellIndex];
        while (atom_j) {
            partner = thisList->ptr[atom_j];
            if (FindStrInList(partner->property->extraProperty, "_HBRMV", 0)) {
                atom_j = sourCellList[atom_j];
                continue;;
            }
            
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
                        checkStatus = FindPair(target, partner, InterColli, tmp, r_2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &accumPoten, NULL, -1);
                    } else if (strcmp(type, "_HB") == 0) {
                        if (!FindElemInList(thisList->ptr[skipAtom]->dynamic->HB.neighbor, atom_j, 1, 2, 0)) {
                            checkStatus = FindPair(target, partner, InterColli, tmp, r_2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &accumPoten, NULL, 1);
                        }
                    } else {
                        printf("!!ERROR!!: type is invalid! %s:%i\n", __FILE__, __LINE__);
                    }
                    if (accumPoten == INFINIT && checkStatus == 1)
                        accumPoten = 0;
                    energyChange += accumPoten;
                }
            }
            atom_j = sourCellList[atom_j];
        }        
    }
    
    return energyChange;
}

int CheckAlphaHB(struct AtomStr* HB_i, struct AtomStr* HB_j) {
    if (HB_i->property->sequence.molNum == HB_j->property->sequence.molNum &&
        ((strcmp(HB_i->property->name, "H") == 0 && strcmp(HB_j->property->name, "O") == 0 && HB_j->property->sequence.resNum == HB_i->property->sequence.resNum - 4) ||
         (strcmp(HB_i->property->name, "O") == 0 && strcmp(HB_j->property->name, "H") == 0 && HB_i->property->sequence.resNum == HB_j->property->sequence.resNum - 4)))
        return TRUE;
    
    return FALSE;
}

static int InHBNeighborList(int *list, int num) {
    NEIGHBORSCAN(list, i) {
        if (list[i] == num) return TRUE;
    }
    return FALSE;
}

void HBAddPartner(struct AtomStr *HB_i, int partner) {
    int *list = HB_i->dynamic->HB.bondConnection;
    if (unlikely(list[0] >= HBCapacity(HB_i))) {
        printf("!!ERROR!!: atom #%i already has %i HB partners! %s:%i\n", HB_i->property->num, list[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    list[++list[0]] = partner;
}

void HBRemovePartner(struct AtomStr *HB_i, int partner) {
    int *list = HB_i->dynamic->HB.bondConnection;
    for (int i = 1; i <= list[0]; i++) {
        if (list[i] == partner) {
            for (; i < list[0]; i++) list[i] = list[i + 1];
            list[list[0]--] = 0;
            return;
        }
    }
    printf("!!ERROR!!: atom #%i has no HB with atom #%i! %s:%i\n", HB_i->property->num, partner, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}

static void CheckHBNeighborRoom(int status, int num) {
    if (unlikely(status >= HB_NEIGHBOR_MAX)) {
        printf("!!ERROR!!: atom #%i has too many HB neighbor partners! %s:%i\n", num, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
}

void AddNeighbor2HB(int *currentStatus, int *list, int *neighborList) {
    NEIGHBORSCAN(neighborList, i) {
        CheckHBNeighborRoom(*currentStatus, neighborList[i]);
        list[++(*currentStatus)] = neighborList[i];
    }
    
    return;
}

void RemoveNeighborFromHB(int *currentStatus, int *list, int *neighborList) {
    NEIGHBORSCAN(neighborList, i) {
        --(*currentStatus);
        ListRefresh(neighborList[i], list, 1, HB_NEIGHBOR_MAX);
    }
    
    return;
}

void AddHB2Neighbor(int *list, int HBNum) {
    int *currentStatus;
    for (int i = 1; i <= list[0]; i ++) {
        currentStatus = &atom[list[i]].dynamic->HBNeighbor.neighborStatus;
        CheckHBNeighborRoom(*currentStatus, list[i]);
        atom[list[i]].dynamic->HBNeighbor.neighborPartner[++(*currentStatus)] = HBNum;
    }
    
    return;
}

void RemoveHBFromNeighbor(int *list, int HBNum) {
    for (int i = 1; i <= list[0]; i ++) {
        --atom[list[i]].dynamic->HBNeighbor.neighborStatus;
        ListRefresh(HBNum, atom[list[i]].dynamic->HBNeighbor.neighborPartner, 1, HB_NEIGHBOR_MAX);
    }
    
    return;
}

int CheckHB_NeighborConnection(int atom_i, int atom_j, int *neighborList_i, int *neighborList_j, int flag) {
    int check = 0;

    if (!flag) { //before HB formation
        if (connectionMap[atom_i][atom_j] & HB_CONNECT ||
            connectionMap[atom_j][atom_i] & HB_CONNECT) {
            check ++;
        }
            
        NEIGHBORSCAN(neighborList_j, i) {
            if (connectionMap[atom_i][atom[atom_j].dynamic->HB.neighbor[i]] & NEIGHBOR_CONNECT ||
                connectionMap[atom[atom_j].dynamic->HB.neighbor[i]][atom_i] & NEIGHBOR_CONNECT) {
                check ++;
            }
        }
        
        NEIGHBORSCAN(neighborList_i, i) {
            if (connectionMap[atom_j][atom[atom_i].dynamic->HB.neighbor[i]] & NEIGHBOR_CONNECT ||
                connectionMap[atom[atom_i].dynamic->HB.neighbor[i]][atom_j] & NEIGHBOR_CONNECT) {
                check ++;
            }
        }
    } else { //after HB formation
        if (!(connectionMap[atom_i][atom_j] & HB_CONNECT) ||
            !(connectionMap[atom_j][atom_i] & HB_CONNECT)) {
            check ++;
        }
            
        NEIGHBORSCAN(neighborList_j, i) {
            if (!(connectionMap[atom_i][atom[atom_j].dynamic->HB.neighbor[i]] & NEIGHBOR_CONNECT) ||
                !(connectionMap[atom[atom_j].dynamic->HB.neighbor[i]][atom_i] & NEIGHBOR_CONNECT)) {
                check ++;
            }
        }
        
        NEIGHBORSCAN(neighborList_i, i) {
            if (!(connectionMap[atom_j][atom[atom_i].dynamic->HB.neighbor[i]] & NEIGHBOR_CONNECT) ||
                !(connectionMap[atom[atom_i].dynamic->HB.neighbor[i]][atom_j] & NEIGHBOR_CONNECT)) {
                check ++;
            }
        }
    }
    
    return check;
}


void SwitchConnection(int atom_i, int atom_j, int *neighborList_i, int *neighborList_j) {
    connectionMap[atom_i][atom_j] ^= HB_CONNECT;
    connectionMap[atom_j][atom_i] ^= HB_CONNECT;
    
    NEIGHBORSCAN(neighborList_j, i) {
        connectionMap[atom_i][atom[atom_j].dynamic->HB.neighbor[i]] ^= NEIGHBOR_CONNECT;
        connectionMap[atom[atom_j].dynamic->HB.neighbor[i]][atom_i] ^= NEIGHBOR_CONNECT;
    }
    
    NEIGHBORSCAN(neighborList_i, i) {
        connectionMap[atom_j][atom[atom_i].dynamic->HB.neighbor[i]] ^= NEIGHBOR_CONNECT;
        connectionMap[atom[atom_i].dynamic->HB.neighbor[i]][atom_j] ^= NEIGHBOR_CONNECT;
    }

    return;
}


void ShowHideAtoms(int *atomList, int flag) {
    struct PropertyStr *ptr;
    
    if (flag == 0) { //hide
        for (int i = 1; i <= atomList[0]; i ++) {
            int num = FindStrInList(atom[atomList[i]].property->extraProperty, "HBRMV", 0);
            if (!num) {
                
                //need further optimizations
                if (FindStrInList(atom[atomList[i]].property->extraProperty, "_HBRMV", 0)) {
                    return;
                }
                
                printf("!!ERROR!!: cannot find the atom needed to be removed. %i!\n", atomList[i]);
                printf("           %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            ptr = atom[atomList[i]].property;
            memset(ptr->extraProperty[num], '\0',
                   sizeof(ptr->extraProperty[num]));
            sprintf(ptr->extraProperty[num], "_HBRMV");
        }
    } else { //show
        for (int i = 1; i <= atomList[0]; i ++) {
            int num = FindStrInList(atom[atomList[i]].property->extraProperty, "_HBRMV", 0);
            if (!num) {
                printf("!!ERROR!!: cannot find the atom needed to be removed. %i!\n", atomList[i]);
                printf("           %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            ptr = atom[atomList[i]].property;
            memset(ptr->extraProperty[num], '\0',
                   sizeof(ptr->extraProperty[num]));
            sprintf(ptr->extraProperty[num], "HBRMV");
        }
        //after unhiding, if the listed atoms would need to be added in the renewlist?
        //what if the listed atom overlap on the other atoms after unhiding?
        //how to handle this situation?
    }

    return;
}

void HBEnergyChangeCombined(struct AtomStr *target, struct AtomStr *partner, double *accumPoten, struct AtomListStr *thisList, int *excl_ap2, int *excl_ap3, int keepType_i, int keepType_j) {
    int atom_i, atom_j, atom_p;
    int cell_neighbor[4], cellIndex;
    int *sourCellList = celllist;
    int checkStatus = 0;
    double cutoffRange = cutoffr * cutoffr;
    double position_i[4], position_j[4], position_p[4];
    double r_temp[4], positionshift[4] = {0};
    double r_ij2, r_pj2;
    double tmp = 0, pot = 0;
    double lowerLimit, upperLimit, lowerPotential, upperPotential;
    struct AtomStr *neighbor;
    
    atom_i = target->property->num;
    atom_p = partner->property->num;
    TRANSFER_VECTOR(position_i, target->dynamic->coordinate);
    TRANSFER_VECTOR(position_p, partner->dynamic->coordinate);
    
    for (int i = 0; i < thisList->stencilNum; ++i) {
        cell_neighbor[3] = target->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = target->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = target->dynamic->cellIndex[1] + thisList->x[i];
        
        for (int j = 1; j <= 3; j++) {
            if (unlikely(cell_neighbor[j] < 0)) {
                cell_neighbor[j] += cellnum[j];
            } else if (unlikely(cell_neighbor[j] >= cellnum[j])) {
                cell_neighbor[j] -= cellnum[j];
            }
        }
        
        cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2]
        + cell_neighbor[2] * cellnum[1]
        + cell_neighbor[1] + 1;
        
        atom_j = sourCellList[atomnum + cellIndex];
        while (atom_j) {
            neighbor = thisList->ptr[atom_j];
            if (FindStrInList(neighbor->property->extraProperty, "_HBRMV", 0)) {
                atom_j = sourCellList[atom_j];
                continue;
            }
            
            if (atom_j != atom_i && atom_j != atom_p) {
                UpdateAtom(neighbor, currenttime);
                TRANSFER_VECTOR(position_j, neighbor->dynamic->coordinate);

                // For target - neighbor
                PBCShift(target, neighbor, positionshift);
                DOT_PLUS(position_j, positionshift, r_temp);
                r_ij2 = CalDistance(position_i, r_temp);
                
                if (r_ij2 <= cutoffRange &&
                    !(connectionMap[atom_i][atom_j] & (BOND_CONNECT | CONSTRAINT_CONNECT))) {
                    // _nonHB for target (all non-bonded neighbors)
                    checkStatus = FindPair(target, neighbor, InterColli, tmp, r_ij2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &pot, NULL, keepType_i ? 1 : -1);
                    if (pot == INFINIT) pot = 0;
                    accumPoten[0] += pot;

                    // _HB for target — exclude HB-neighbor C (its contribution is captured by ap[4])
                    if (!InHBNeighborList(excl_ap2, atom_j)) {
                        checkStatus = FindPair(target, neighbor, InterColli, tmp, r_ij2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &pot, NULL, 1);
                        if (pot == INFINIT) pot = 0;
                        accumPoten[2] += pot;
                    }
#ifdef HB_DEBUG
                    else {
                        // Print what C would contribute to ap[2] if not excluded
                        double c_pot;
                        int cs = FindPair(target, neighbor, InterColli, tmp, r_ij2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &c_pot, NULL, 1);
                        if (c_pot == INFINIT) c_pot = 0;
                        printf("[HB_EXCL] t=%.4f H#%d O#%d: C(atom#%d) r_HC=%.3f excluded from ap[2]; C_HB_colli=%.4f\n",
                               currenttime, atom_i, atom_p, atom_j, sqrt(r_ij2), c_pot);
                        fflush(stdout);
                    }
#endif
                }

            }
            atom_j = sourCellList[atom_j];
        }
    }

    // Second loop: scan partner's (O's) 27 cells for ap[1] and ap[3].
    // Must be a separate scan so that atoms in O's cell neighborhood but NOT
    // in H's neighborhood (different adjacent cells) are not missed.
    for (int i = 0; i < thisList->stencilNum; ++i) {
        cell_neighbor[3] = partner->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = partner->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = partner->dynamic->cellIndex[1] + thisList->x[i];

        for (int j = 1; j <= 3; j++) {
            if (unlikely(cell_neighbor[j] < 0)) {
                cell_neighbor[j] += cellnum[j];
            } else if (unlikely(cell_neighbor[j] >= cellnum[j])) {
                cell_neighbor[j] -= cellnum[j];
            }
        }

        cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2]
        + cell_neighbor[2] * cellnum[1]
        + cell_neighbor[1] + 1;

        atom_j = sourCellList[atomnum + cellIndex];
        while (atom_j) {
            neighbor = thisList->ptr[atom_j];
            if (FindStrInList(neighbor->property->extraProperty, "_HBRMV", 0)) {
                atom_j = sourCellList[atom_j];
                continue;
            }

            if (atom_j != atom_i && atom_j != atom_p) {
                UpdateAtom(neighbor, currenttime);
                TRANSFER_VECTOR(position_j, neighbor->dynamic->coordinate);

                // For partner (O) - neighbor
                PBCShift(partner, neighbor, positionshift);
                DOT_PLUS(position_j, positionshift, r_temp);
                r_pj2 = CalDistance(position_p, r_temp);

                if (r_pj2 <= cutoffRange &&
                    !(connectionMap[atom_p][atom_j] & (BOND_CONNECT | CONSTRAINT_CONNECT))) {
                    // _nonHB for partner (all non-bonded neighbors)
                    checkStatus = FindPair(partner, neighbor, InterColli, tmp, r_pj2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &pot, NULL, keepType_j ? 1 : -1);
                    if (pot == INFINIT) pot = 0;
                    accumPoten[1] += pot;

                    // _HB for partner — exclude HB-neighbor N (its contribution is captured by ap[5])
                    if (!InHBNeighborList(excl_ap3, atom_j)) {
                        checkStatus = FindPair(partner, neighbor, InterColli, tmp, r_pj2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &pot, NULL, 1);
                        if (pot == INFINIT) pot = 0;
                        accumPoten[3] += pot;
                    }
#ifdef HB_DEBUG
                    else {
                        // Print what N would contribute to ap[3] if not excluded
                        double n_pot;
                        int cs = FindPair(partner, neighbor, InterColli, tmp, r_pj2, &lowerLimit, &upperLimit, &lowerPotential, &upperPotential, &n_pot, NULL, 1);
                        if (n_pot == INFINIT) n_pot = 0;
                        printf("[HB_EXCL] t=%.4f H#%d O#%d: N(atom#%d) r_ON=%.3f excluded from ap[3]; N_HB_colli=%.4f\n",
                               currenttime, atom_i, atom_p, atom_j, sqrt(r_pj2), n_pot);
                        fflush(stdout);
                    }
#endif
                }
            }
            atom_j = sourCellList[atom_j];
        }
    }
}
