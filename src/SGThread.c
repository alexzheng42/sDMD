//
//  SGThread.c
//  sDMD
//
//  Created by Size Zheng on 9/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

void CGResidueRenew(int *renewList, struct AtomListStr *thisList);


void SingleThread(void) {
	Prepare();
	SingleThreadRun(NULL, NULL);
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
    InitBackupEvents();
    for (int i = 1; i <= atomnum; i ++) {
        atomList.ptr[i] = &atom[i];
    }
    
    if (strcmp("new", neworcontinue) == 0) {
        SaveData(sysInfo, fileList);
        if (strcmp(surfInfo.style, "implicit") == 0)
            SaveData(wallInfo, fileList);
        if (obstObj.mark)
            SaveData(obstInfo, fileList);
        if (CG.mark)
            SaveData(CGInfo,   fileList);
        
#ifndef GEL
        BuildSurfaceConnection();
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
        RebaseTime();
        outputrecord = currenttime + outputrate;
    } else {
        printf("!!ERROR!!: the simulation configuration is invalid!\n");
        printf("           it only support \"new\" or \"continue\"\n");
        printf("%s:%i\n", __FILE__, __LINE__);
    }
    
    CreateCBT();
    return;
}


// The event loop of every run. REMD plugs its replica exchange in through the hooks, so
// it uses exactly the same event handling (backup events, CG renewal, stuck detection):
//   afterEvent: called after each executed event (after the output files are written)
//   beforeFree: called at the end, after all atoms are advanced to the final time
void SingleThreadRun(void (*afterEvent)(void), void (*beforeFree)(void)) {
    int hazardType;
    int renewList[256] = {0};
    double processratio, gap = 0;
    struct AtomListStr *thisList = &atomList;
    struct AtomStr  *target;
    struct AtomStr  *partner;
    int error_count = 0;
    int hazard_count = 0;
    long int last_warningsum = warningsum;
    int warning_count = 0;
    
    while (currenttime <= timestep) {
        
        AssignJob(&target, &partner, renewList, thisList);
        if (thisList->timeInr < 1E-14) {
            error_count++;
        }

        currenttime += thisList->timeInr;
        UpdateAtom(target, currenttime);
        if (partner != NULL) {
            UpdateAtom(partner, currenttime);
        }

        hazardType = HazardCheck(target, partner, thisList);

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
            if (renewList[0] == 2 && !FindStrInList(partner->property->extraProperty, "_CONS", 0))
                partner->dynamic->event.counter ++;
            
            // Reset loop-prevention counters upon successful event execution;
            // a run of executed events that never advance the time is still a stall
            if (thisList->timeInr >= 1E-14) error_count = 0;
            hazard_count = 0;
            warning_count = 0;
            last_warningsum = warningsum;

			if (CG.mark && target->dynamic->event.eventType == CGSu_Event) {
				//if the target atom is a fixing atom and the event is a CGEvent, then all the atoms in the neighboring cells will be renewed
				if (FindStrInList(target->property->extraProperty, "_CONS", 0)) {
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
            
            if (strcmp(wallDyn.mark, "no") != 0) DoWallDyn();
            
            //save data to the output files
            if (unlikely(currenttime >= outputrecord)) {
                UpdateAllAtoms(currenttime);
                for (int i = 0; i < lenFileType; i ++) {
                    if (fileList[i].mark) {
                        SaveData(i, fileList);
                    }
                }
                outputrecord += outputrate;
            }

            if (unlikely(currenttime >= TIME_BASE_UNIT)) RebaseTime();

            if (afterEvent) afterEvent();
        } else {
#ifdef DEBUG_PRINTF
            printf(", denied, self/partner changed!\n");
#endif
            hazard_count++;
            if (hazard_count > 10000) {
                printf("!!ERROR!!: too many consecutive HazardCheck failures! Potential infinite loop in lazy update.\n");
                printf("           %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }

            renewList[0] = 1;
            renewList[2] = 0;
            renewList[3] = 0;
            
            if (TakeBackupEvent(target, thisList)) {
                countBackup ++;
            } else {
                Predict(renewList, thisList);
                countReCal ++;
            }
            UpdateCBT(renewList);
        }

        if (warningsum > last_warningsum) {
            warning_count += (int)(warningsum - last_warningsum);
            last_warningsum = warningsum;
        }

        if (unlikely(error_count > 1000000)) {
            printf("!!ERROR!!: the time has not advanced for %i events at t = %.6lf; the simulation is stuck.\n", error_count, REAL_TIME(currenttime));
            printf("           %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        if (warningsum > 1000 && error_count > 10) {
            printf("!!ERROR!!: warnings (%li) exceeded 1000 and errors exceeded 10 — simulation appears stuck.\n", warningsum);
            printf("           %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
    
    UpdateAllAtoms(timestep);
    if (beforeFree) beforeFree();
    FreeVariables();
    return;
}

// Event times are absolute, so their precision drops as the run gets longer: past
// t = 2^23 events land off their step boundary by more than ZERO and the run stalls
// in zero-time events. So the time origin moves forward by whole TIME_BASE_UNITs
// (exact in floating point, and the event order is unchanged); the real time is
// REAL_TIME(currenttime) = timeBaseCount * TIME_BASE_UNIT + currenttime.
void RebaseTime(void) {
    if (currenttime < TIME_BASE_UNIT) return;
    long long units = (long long)(currenttime / TIME_BASE_UNIT);
    double shift = units * TIME_BASE_UNIT;

    UpdateAllAtoms(currenttime); // so that lastUpdateTime shifts exactly as well
    for (int i = 1; i <= atomnum; i++) {
        if (atom[i].dynamic->event.time < INFINIT) atom[i].dynamic->event.time -= shift;
        atom[i].dynamic->lastUpdateTime -= shift;
    }
    ShiftBackupEvents(shift);
    if (flow.force.timeRec) {
        for (int i = 1; i <= atomnum; i++) flow.force.timeRec[i] -= shift;
    }
    currenttime    -= shift;
    oldcurrenttime -= shift;
    timestep       -= shift;
    outputrecord   -= shift;
    timeBaseCount  += units;
}

void Predict(int *renewList, struct AtomListStr *thisList) {
    struct AtomStr *target = NULL;
    
    ResetTarget(renewList, thisList);
    for (int n = 1; renewList[n]; n ++) {
        target = thisList->ptr[renewList[n]];
        
        if (FindStrInList(target->property->extraProperty, "_CONS", 0)) {
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
        
        if (strcmp(surfInfo.style, "implicit") == 0) {
            WallTime(target);
        }
        
        if (obstObj.mark && CG.mark != 2) { // no CG, every atom will calculate the time
            ObstTime(target);
        }
        
        if (obstObj.mark && CG.mark == 2 && // CG surface, only the CG beads will calculate the time
            (FindStrInList(target->property->extraProperty, "CG", 0) ||
             FindStrInList(target->property->extraProperty, "CONS", 0))) {
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
int HazardCheck(struct AtomStr *target, struct AtomStr *partner, struct AtomListStr *thisList) {
    int partnerNum = -1;

    if (partner != NULL) partnerNum = partner->property->num;

    // Check if the partner's state has changed since this event was predicted.
    // target->dynamic->event.partnerCounter stores the partner's counter at prediction time.
    if (partnerNum > 0 &&
        target->dynamic->event.partnerCounter < partner->dynamic->event.counter)
        return -1;

    // In a single-threaded environment, target's own counter and partner identity 
    // shouldn't change between AssignJob and HazardCheck, but we keep the partner 
    // check for consistency with the prediction logic.
    if (target->dynamic->event.partner != (partner ? partner->property->num : INVALID))
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
    thisList->timeInr   = thisList->target->dynamic->event.time - currenttime;
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
