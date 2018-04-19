//
//  ThreadProcess.c
//  Update
//
//  Created by Size Zheng on 3/28/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

int CheckAlphaHB(struct AtomStr* HB_i, struct AtomStr* HB_j);


int ThreadProcess(void) {
    
    thrInfo = (struct ThreadInfoStr*)malloc(threadNum * sizeof(struct ThreadInfoStr));
    for (int i = 0; i < threadNum; i ++) {
        thrInfo[i].threadID = i;
        thrInfo[i].threadRenewList = (int *)calloc(5, sizeof(int));
        thrInfo[i].threadRenewList[0] = 1;
        thrInfo[i].threadRenewList[1] = 0;
    }
    
    if (REMDInfo.flag  == 1) { //REMD is requested
        REMD();
    } else if (codeNum == 1 && threadNum == 1) { //single thread - one core
        SingleThread();
    } else if (codeNum == 2 && threadNum > 1) { //multiple threads - one master, one or more slaves
        MSThread();
    } else if (codeNum == 3) { //multiple threads - two or more masters
        FAThread();
    } else {
        printf("!!ERROR!!: the combination of sDMD method and thread number is invalid! %s:%i\n", __FILE__, __LINE__);
    }
    
    for (int i = 0; i < threadNum; ++i) {
        free(thrInfo[i].threadRenewList);
    }
    free(thrInfo);
    
    return 0;
}


struct ThreadStr* InitializeThread(int tid, struct AtomStr *atomLib) {
    struct ThreadStr *thisThread = (struct ThreadStr *)calloc(1, sizeof(struct ThreadStr));

    thisThread->oldTarget.dynamic = (struct DynamicStr *)calloc(1, sizeof(struct DynamicStr));
    thisThread->oldTarget.eventList = (struct EventListStr *)calloc(1, sizeof(struct EventListStr));
    
    thisThread->oldPartner.dynamic = (struct DynamicStr *)calloc(1, sizeof(struct DynamicStr));
    thisThread->oldPartner.eventList = (struct EventListStr *)calloc(1, sizeof(struct EventListStr));
    
    thisThread->raw = (struct AtomStr **)calloc(atomnum + 1, sizeof(struct AtomStr *));
    thisThread->listPtr = (struct AtomStr **)calloc(atomnum + 1, sizeof(struct AtomStr *));
    thisThread->tid = tid;
    thisThread->finishWork = 0;
    thisThread->getWork = 0;
    
    for (int i = 0; i <= atomnum; i++) {
        thisThread->raw[i] = &atomLib[i];
    }
    listInit(&thisThread->atomList);
    
    return  thisThread;
}


void FirstRun(struct ThreadStr *thisThread) {
    int n;
    int tmpList[3] = {1, 0, 0};
    
#ifndef GEL
    SDEnergyMin(50000, thisThread);
#endif
    
    PBC("full", NULL, thisThread); //check if any of the coordinates would satisfy PBC
    LinkList("full", NULL, thisThread); //create link list
    
    //=================================================
    //at the 1st step, calculate the time of all events
    for (n = 1; n <= atomnum; ++n) {
        thisThread->atomNum = n;
        tmpList[1] = n;
        
        AssignThread(tmpList, thisThread);
        
        thisThread->newTarget->dynamic->event.eventType = Invd_Event;
        thisThread->newTarget->dynamic->event.time = INFINIT;
        thisThread->newTarget->dynamic->event.partner = INVALID;
        
        BondTime(thisThread->newTarget, thisThread);
        CollisionTime(thisThread->newTarget, thisThread);
        PBCandCrossCellTime(thisThread->newTarget);

#ifndef GEL
        ThermostatTime(thisThread->newTarget);
#endif
        
#ifdef HYDROGEN_BOND
        HBTime(thisThread->newTarget, thisThread);
        HBNeighborTime(thisThread->newTarget, thisThread);
#endif
        
        if (strcmp(wallExist, "no")) {
            WallTime(thisThread->newTarget);
        }
        
        if (obstObj.mark) {
            ObstTime(thisThread->newTarget);
        }
        
        if (tunlObj.mark) {
            TunnelTime(thisThread->newTarget);
        }
        
        if (flow.mark == 3) {
            ChargeTime(thisThread->newTarget);
        }
        
        AtomDataCpy(thisThread->raw[n], thisThread->newTarget, 0);
    }
    //=================================================
    
    SchedulingRefresh(thisThread); //generate a single-event binary tree for the coming events
}


int ProcessEvent(int *threadRenewList, struct ThreadStr *thisThread) {
    int hazardType = 0;
    double timeIncr;
    struct AtomStr *newTarget, *newPartner;
    struct AtomStr *oldTarget, *oldPartner;
    
    newTarget = thisThread->newTarget;
    newPartner = thisThread->newPartner;
    oldTarget = &thisThread->oldTarget;
    oldPartner = &thisThread->oldPartner;
    
    //if newPartner = NULL, it means in this event there is no partner associate with the target atom
    //then all the pointers related to "partner" will be NULL
    //we can also by checking if threadRenewList[2] > 0 to find if there would be a partner in this event
#ifdef DEBUG_IT
    if (!(newPartner == NULL && threadRenewList[2] == 0) &&
        !(newPartner != NULL && threadRenewList[2] > 0)) {
        printf("!!ERROR!!: partner state is invalid! %s:%i\n", __FILE__, __LINE__);
    }
#endif
    hazardType = HazardCheck(oldTarget, (newPartner != NULL ? oldPartner : NULL),
                             thisThread->listPtr[oldTarget->dynamic->HB.neighbor],
                             (newPartner != NULL ? thisThread->listPtr[oldPartner->dynamic->HB.neighbor] : NULL),
                             thisThread);
    
    if (hazardType == 0) {
        timeIncr = thisThread->newTarget->dynamic->event.time;
        UpdateData(timeIncr, "partial", thisThread);
        TimeForward(timeIncr, "partial", thisThread);
        
        oldTarget->dynamic->HB.interactionType = DoEvent(thisThread); //calculate the after-event velocities
        
        newTarget->dynamic->event.counter ++;
        if (threadRenewList[2] > 0) {
            newPartner->dynamic->event.counter ++;
        }
        
        Predict(threadRenewList, thisThread); //pre-predict the after-event event
    }
    
    return hazardType;
}


void CommitEvent(struct AtomStr **destLibrary, struct AtomStr *newTargetAtom, struct AtomStr *newPartnerAtom, struct AtomStr *oldTargetAtom, struct AtomStr *oldPartnerAtom, struct AtomStr *oldTargetNeighbor, struct AtomStr *oldPartnerNeighbor) {
    int typeNum;
    struct AtomStr *HB_i, *HB_j;
    
    AtomDataCpy(destLibrary[newTargetAtom->property->num], newTargetAtom, 0);
    if (newPartnerAtom)
        AtomDataCpy(destLibrary[newPartnerAtom->property->num], newPartnerAtom, 0);
    
    typeNum = oldTargetAtom->dynamic->event.eventType;
    
    switch (typeNum) {
        case Coli_Event:
            collisioneventsum ++;
            break;
        case Bond_Event:
            bondeventsum ++;
            break;
        case HB_Event:
            HB_i = oldTargetAtom;
            HB_j = oldPartnerAtom;
            
            destLibrary[HB_i->dynamic->HB.neighbor]->dynamic->HBNeighbor = oldTargetNeighbor->dynamic->HBNeighbor;
            destLibrary[HB_j->dynamic->HB.neighbor]->dynamic->HBNeighbor = oldPartnerNeighbor->dynamic->HBNeighbor;
            
            int atom_i = HB_i->property->num;
            int atom_j = HB_j->property->num;
            
            if (oldTargetAtom->dynamic->HB.interactionType == HBForm) {
                if (CheckAlphaHB(HB_i, HB_j))
                    alphaHBformed ++;
                
#ifdef DEBUG_IT
                if (connectionMap[atom_i][atom_j] & HB_CONNECT ||
                    connectionMap[atom_j][atom_i] & HB_CONNECT ||
                    connectionMap[atom_i][HB_j->dynamic->HB.neighbor] & NEIGHBOR_CONNECT ||
                    connectionMap[HB_j->dynamic->HB.neighbor][atom_i] & NEIGHBOR_CONNECT ||
                    connectionMap[atom_j][HB_i->dynamic->HB.neighbor] & NEIGHBOR_CONNECT ||
                    connectionMap[HB_i->dynamic->HB.neighbor][atom_j] & NEIGHBOR_CONNECT) {
                    printf("!!ERROR!!: HB connection map has something wrong! %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
#endif
                
                //adjust the connection map
                connectionMap[atom_i][atom_j] ^= HB_CONNECT;
                connectionMap[atom_j][atom_i] ^= HB_CONNECT;
                connectionMap[atom_i][HB_j->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
                connectionMap[HB_j->dynamic->HB.neighbor][atom_i] ^= NEIGHBOR_CONNECT;
                connectionMap[atom_j][HB_i->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
                connectionMap[HB_i->dynamic->HB.neighbor][atom_j] ^= NEIGHBOR_CONNECT;
                
                HBnumformed ++;
            } else if (oldTargetAtom->dynamic->HB.interactionType == HBBreak) {
                if (CheckAlphaHB(HB_i, HB_j))
                    alphaHBformed --;
                
#ifdef DEBUG_IT
                if (!(connectionMap[atom_i][atom_j] & HB_CONNECT) ||
                    !(connectionMap[atom_j][atom_i] & HB_CONNECT) ||
                    !(connectionMap[atom_i][HB_j->dynamic->HB.neighbor] & NEIGHBOR_CONNECT) ||
                    !(connectionMap[HB_j->dynamic->HB.neighbor][atom_i] & NEIGHBOR_CONNECT) ||
                    !(connectionMap[atom_j][HB_i->dynamic->HB.neighbor] & NEIGHBOR_CONNECT) ||
                    !(connectionMap[HB_i->dynamic->HB.neighbor][atom_j] & NEIGHBOR_CONNECT)) {
                    printf("!!ERROR!!: HB connection map has something wrong! %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
#endif
                
                //adjuct the connection map
                connectionMap[atom_i][atom_j] ^= HB_CONNECT;
                connectionMap[atom_j][atom_i] ^= HB_CONNECT;
                connectionMap[atom_i][HB_j->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
                connectionMap[HB_j->dynamic->HB.neighbor][atom_i] ^= NEIGHBOR_CONNECT;
                connectionMap[atom_j][HB_i->dynamic->HB.neighbor] ^= NEIGHBOR_CONNECT;
                connectionMap[HB_i->dynamic->HB.neighbor][atom_j] ^= NEIGHBOR_CONNECT;
                
                HBnumformed --;
                
#ifdef DEBUG_IT
                if (HBnumformed < alphaHBformed) {
                    printf("!WARNING!: total HB number is smaller than the alpha HB number! %s:%i\n", __FILE__, __LINE__);
                }
#endif
            }
            HBeventsum ++;
            break;
            
        case HBNe_Event:
            HBNeighboreventsum ++;
            break;
        case Ther_Event:
            thermostateventsum ++;
            break;
        case Wall_Event:
        case Obst_Event:
        case Tunl_Event:
            walleventsum ++;
            break;
        case PCC_Event:
            RenewCellList(newTargetAtom, oldTargetAtom);
            pbcandcrosseventsum ++;
            break;
        case Lagv_Event:
            //LangevinEventSum ++;
            break;
            
        default:
            break;
    }
}


void Predict(int *renewList, struct ThreadStr *thisThread) {
    struct AtomStr *targetAtom = NULL;
    
    ResetTarget(renewList, thisThread);
    for (int n = 1; n <= renewList[0]; n ++) {
        targetAtom = thisThread->listPtr[renewList[n]];
        
        BondTime(targetAtom, thisThread);
        CollisionTime(targetAtom, thisThread);
        PBCandCrossCellTime(targetAtom);
        
#ifndef GEL
        ThermostatTime(targetAtom);
#endif
        
#ifdef HYDROGEN_BOND
        HBTime(targetAtom, thisThread);
        HBNeighborTime(targetAtom, thisThread);
#endif
        
        if (strcmp(wallExist, "no")) {
            WallTime(targetAtom);
        }
        
        if (obstObj.mark) {
            ObstTime(targetAtom);
        }
        
        if (tunlObj.mark) {
            TunnelTime(targetAtom);
        }
        
        if (flow.mark == 3) {
            ChargeTime(targetAtom);
        }
    }
}


//!this function should be called before DoEvent!
int HazardCheck(struct AtomStr *oldTargetAtom, struct AtomStr *oldPartner, struct AtomStr *oldTargetNeighbor, struct AtomStr *oldPartnerNeighbor, struct ThreadStr *thisThread) {
    int targetNum, partnerNum;
    //0: no hazard; -1: partner changed
    
    targetNum = oldTargetAtom->property->num;
    if (oldPartner != NULL) partnerNum = oldPartner->property->num;
    else partnerNum = 0;
    
    //partner changed
    if (partnerNum > 0 &&
        (oldTargetAtom->dynamic->event.partnerCounter < thisThread->raw[partnerNum]->dynamic->event.counter ||
         oldPartner->dynamic->event.counter < thisThread->raw[partnerNum]->dynamic->event.counter))
        return -1;
    
    //self changed
    if (oldTargetAtom->dynamic->event.counter < thisThread->raw[targetNum]->dynamic->event.counter ||
        oldTargetAtom->dynamic->event.partner != thisThread->raw[targetNum]->dynamic->event.partner)
        return -1;
    
#ifdef HYDROGEN_BOND
    //for HB event, the exact interaction types (HB formation or breaking) will be determined during DoEvent, before which if the neighbor atoms changed,
    //the further determination would be invalid.
    //also, the neighbor atoms changes may effect the prediction
    //So this situation would be included in partner changes case.
    if (oldTargetAtom->dynamic->event.eventType == HB_Event &&
        (thisThread->raw[oldTargetAtom->dynamic->HB.neighbor]->dynamic->event.counter > oldTargetNeighbor->dynamic->event.counter
         || thisThread->raw[oldPartner->dynamic->HB.neighbor]->dynamic->event.counter > oldPartnerNeighbor->dynamic->event.counter))
        return -1;
#endif
    
    return 0;
}


int CheckAlphaHB(struct AtomStr* HB_i, struct AtomStr* HB_j) {
    if (HB_i->property->sequence.proteinNum == HB_j->property->sequence.proteinNum &&
        ((strcmp(HB_i->property->name, "H") == 0 && strcmp(HB_j->property->name, "O") == 0 && HB_j->property->sequence.aminoacidNum == HB_i->property->sequence.aminoacidNum - 4) ||
         (strcmp(HB_i->property->name, "O") == 0 && strcmp(HB_j->property->name, "H") == 0 && HB_i->property->sequence.aminoacidNum == HB_j->property->sequence.aminoacidNum - 4)))
        return TRUE;
    
    return FALSE;
}

