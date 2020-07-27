//
//  SGThread.c
//  sDMD
//
//  Created by Size Zheng on 9/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

void SingleThreadRun(void);
void CGResidueRenew(int *renewList, struct AtomListStr *thisList);


void SingleThread(void) {
	Prepare();
	SingleThreadRun();
	return;
}

void Prepare(void) {
    char *extraName = NULL;
    if (REMDInfo.flag) {
        extraName = REMDInfo.REMD_ExtraName;
    }
    
    InitializeFiles(extraName, fileList);
    SaveData(lgf, fileList);

    atomList.ptr = (struct AtomStr **)calloc(atomnum + 1, sizeof(struct AtomStr *));
    for (int i = 1; i <= atomnum; i ++) {
        atomList.ptr[i] = &atom[i];
    }
    
    if (strcmp("new", neworcontinue) == 0) {
        SaveData(sysInfo, fileList);
        if (strcmp(wallExist, "smooth") == 0) SaveData(wallInfo, fileList);
        if (obstObj.mark) SaveData(obstInfo, fileList);
        if (CG.mark)      SaveData(CGInfo, fileList);
        
#ifndef GEL
        SDEnergyMin(20000, &atomList);
#endif
        PBC("full", NULL, &atomList); //check if any of the coordinates would satisfy PBC
        LinkList("full", NULL, NULL, &atomList); //create link list
        
        int *renewList = calloc(atomnum + 2, sizeof(int));
        for (int i = 1; i <= atomnum; i ++) {
            renewList[i] = i;
        }
        
        Predict(renewList, &atomList);
        free(renewList);
    } else if (strcmp("continue", neworcontinue) == 0) {
        outputrecord = currenttime + outputrate;
    } else {
        printf("!!ERROR!!: the simulation configuration is invalid!\n");
        printf("           it only support \"new\" or \"continue\"\n");
        printf("%s:%i\n", __FILE__, __LINE__);
    }
    
    CreateCBT();
    return;
}


void SingleThreadRun(void) {
    int hazardType;
    int renewList[256] = {0};
    double processratio, gap = 0;
    struct AtomListStr *thisList = &atomList;
    struct AtomStr  *target;
    struct AtomStr  *partner;
    
    while (currenttime <= timestep) {
        
        AssignJob(&target, &partner, renewList, thisList);
        UpdateData(thisList->timeInr);
        TimeForward(thisList->timeInr); //update the node time
        currenttime += thisList->timeInr;

        hazardType = HazardCheck(target, partner,
                                 thisList->ptr[target->dynamic->HB.neighbor],
                                 (partner != NULL ? thisList->ptr[partner->dynamic->HB.neighbor] : NULL),
                                 thisList);

#ifdef DEBUG_PRINTF
        printf("frame = %4li, target atom = %5i, partner = %5i, %6s:%6s, type = %2i",
               frame, target->property->num, target->dynamic->event.partner,
               target->property->name,
               partner == NULL ? "NULL" : partner->property->name,
               target->dynamic->event.eventType);
        fflush(stdout);
#endif
      
        if (!hazardType) {
#ifdef DEBUG_PRINTF
            printf(", executed!\n");
#endif
            frame++;
            DoEvent(thisList); //calculate the after-event velocities
            target->dynamic->event.counter ++;
            if (renewList[0] == 2 && strcmp(partner->property->extraProperty[1], "_CONS"))
                partner->dynamic->event.counter ++;
            
			if (CG.mark && target->dynamic->event.eventType == CGSu_Event) {
				//if the target atom is a fixing atom and the event is a CGEvent, then all the atoms in the neighboring cells will be renewed
				if (strcmp(target->property->extraProperty[1], "_CONS") == 0) {
					FixAssignRenewList(renewList, thisList);
				}

				//if the CG type is "residue", all the atoms in a residue will be renewed after a CG collision
				if (strcmp(CG.type[1], "residue") == 0 &&
					strcmp(target->property->name, "CA") == 0) {
					CGResidueRenew(renewList, thisList);
				}
			}

            Predict(renewList, thisList); //pre-predict the after-event event
            UpdateCBT(renewList);
            
            processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
            if (unlikely(processratio >= gap + 0.01)) {
                printf("Process=%8.2lf%%\r", processratio);
                fflush(stdout);
                gap += 0.01;
            }
            
            if (strcmp(wallDyn.mark, "no")) DoWallDyn();
            
            //save data to the output files
            if (unlikely(currenttime >= outputrecord)) {
                for (int i = 0; i < lenFileType; i ++) {
                    if (fileList[i].mark) {
                        SaveData(i, fileList);
                    }
                }
                outputrecord += outputrate;
            }
        } else {
#ifdef DEBUG_PRINTF
            printf(", denied, self/partner changed!\n");
#endif
            renewList[0] = 1;
            renewList[2] = 0;
            renewList[3] = 0;
            
            Predict(renewList, thisList);
            UpdateCBT(renewList);
            countReCal ++;
        }
    }
    
    FreeVariables();
    return;
}

void Predict(int *renewList, struct AtomListStr *thisList) {
    struct AtomStr *target = NULL;
    
    ResetTarget(renewList, thisList);
    for (int n = 1; renewList[n]; n ++) {
        target = thisList->ptr[renewList[n]];
        
        if (strcmp(target->property->extraProperty[1], "_CONS") == 0) {
            target->dynamic->event.time = INFINIT;
            continue;
        }
        
#ifndef GEL
        BondTime(target, thisList);
#endif
        
        CollisionTime(target, thisList);
        PBCandCrossCellTime(target);
        ThermostatTime(target);
        
#ifdef HYDROGEN_BOND
        HBTime(target, thisList);
        HBNeighborTime(target, thisList);
#endif
        
        if (strcmp(wallExist, "no")) {
            WallTime(target);
        }
        
        if (obstObj.mark && CG.mark != 2) { // no CG, every atom will calculate the time
            ObstTime(target);
        }
        
        if (obstObj.mark && CG.mark == 2 && // CG surface, only the CG beads will calculate the time
            (strcmp(target->property->extraProperty[0], "CG") == 0 ||
             strcmp(target->property->extraProperty[1], "CONS") == 0)) {
                CGSurfaceTime(target);
            }
        
        if (tunlObj.mark) {
            TunnelTime(target);
        }
        
        if (flow.mark == 3) {
            ChargeTime(target);
        }
        
        if (SphObstObj.mark) {
            SphObstTime(target);
        }
    }
    
    return;
}


//!this function should be called before DoEvent!
int HazardCheck(struct AtomStr *target, struct AtomStr *partner, struct AtomStr *targetNeighbor, struct AtomStr *partnerNeighbor, struct AtomListStr *thisList) {
    int targetNum, partnerNum;
    //0: no hazard; -1: partner changed
    
    targetNum = target->property->num;
    if (partner != NULL) partnerNum = partner->property->num;
    else partnerNum = 0;
    
    //partner changed
    if (partnerNum > 0 &&
        ( target->dynamic->event.partnerCounter < thisList->ptr[partnerNum]->dynamic->event.counter ||
         partner->dynamic->event.counter        < thisList->ptr[partnerNum]->dynamic->event.counter))
        return -1;
    
    //self changed
    if (target->dynamic->event.counter  < thisList->ptr[targetNum]->dynamic->event.counter ||
        target->dynamic->event.partner != thisList->ptr[targetNum]->dynamic->event.partner)
        return -1;
    
    return 0;
}

void AssignJob(struct AtomStr **target, struct AtomStr **partner, int *renewList, struct AtomListStr *thisList) {
    renewList[0] = 1;
    renewList[1] = 0;
    renewList[2] = 0;
    renewList[3] = 0;
    
    thisList->targetNum = CBT.node[1];
    thisList->target    = thisList->ptr[thisList->targetNum];
    thisList->timeInr   = thisList->target->dynamic->event.time;
    renewList[1]        = thisList->targetNum;
    
    if (thisList->target->dynamic->event.partner > 0) {
        thisList->partner = thisList->ptr[thisList->target->dynamic->event.partner];
        renewList[++renewList[0]] = thisList->target->dynamic->event.partner;
    } else {
        thisList->partner = NULL;
    }
    
    *target  = thisList->target;
    *partner = thisList->partner;
    
    return;
}

void FreeVariables(void) {
    return;
}
