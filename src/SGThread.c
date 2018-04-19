//
//  SGThread.c
//  sDMD
//
//  Created by Size Zheng on 9/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

void SingleThreadRun(struct ThreadInfoStr *threadInfo);
void FreeVariables(void);


void SingleThread(void) {
    thread = (struct ThreadStr **)calloc(1, sizeof(struct ThreadStr *));
    thread[0] = InitializeThread(0, atom);
    thread[0]->fileList = InitializeFiles(NULL, fileList);
    SaveData(lgf, thread[0]);
    if (strncmp("new", neworcontinue, 1) == 0) {
        SaveData(sysInfo, thread[0]);
        FirstRun(thread[0]);
    } else if (strncmp("continue", neworcontinue, 1) == 0) {
        outputrecord = currenttime + outputrate;
        SchedulingRefresh(thread[0]);
    }
    thrInfo[0].threadRenewList[1] = eventToCommit = SchedulingNextEvent(thread[0]);
    thread[0]->atomNum = eventToCommit;
    if (thread[0]->raw[eventToCommit]->dynamic->event.partner > 0) {
        thrInfo[0].threadRenewList[++thrInfo[0].threadRenewList[0]] = thread[0]->raw[eventToCommit]->dynamic->event.partner;
    }
    AssignThread(thrInfo[0].threadRenewList, thread[0]);
    SingleThreadRun(&thrInfo[0]);
}


void SingleThreadRun(struct ThreadInfoStr *threadInfo) {
    int hazardType;
    int tid = threadInfo->threadID;
    int *threadRenewList = threadInfo->threadRenewList;
    double timeIncr;
    double processratio, gap = 0;
    struct ThreadStr *thisThread = thread[tid];
    struct AtomStr **newTarget = &thisThread->newTarget;
    struct AtomStr **newPartner = &thisThread->newPartner;
    struct AtomStr *oldTarget = &thisThread->oldTarget;
    struct AtomStr *oldPartner = &thisThread->oldPartner;
    
    while (currenttime <= timestep) {
        
        if (frame == 920) {
            printf("");
        }
        
        hazardType = ProcessEvent(threadRenewList, thisThread);

#ifdef DEBUG_PRINTF
        printf("frame = %4li, target atom = %5i, partner = %5i, %6s:%6s, type = %2i",
               frame, thisThread->atomNum, thisThread->oldTarget.dynamic->event.partner,
               thisThread->oldTarget.property->name,
               *newPartner == NULL ? "NULL" : thisThread->oldPartner.property->name,
               thisThread->oldTarget.dynamic->event.eventType);
        fflush(stdout);
#endif
        
        if (hazardType < 0) {
#ifdef DEBUG_PRINTF
            printf(", denied, partner changed!\n");
#endif
            threadRenewList[0] = 1;
            threadRenewList[2] = 0;
            
            AssignThread(threadRenewList, thisThread);
            Predict(threadRenewList, thisThread);

            AtomDataCpy(thisThread->raw[threadRenewList[1]], thisThread->listPtr[threadRenewList[1]], 0);
            SchedulingDelete(threadRenewList, thisThread);
            SchedulingAdd(threadRenewList, thisThread);
            
            AssignJob(threadRenewList, thisThread);
            AssignThread(threadRenewList, thisThread);
            
            countReCal ++;

        } else {
#ifdef DETAILS
            PrintData(threadRenewList);
#endif
#ifdef DEBUG_PRINTF
            printf(", executed!\n");
#endif
            timeIncr = thisThread->raw[threadRenewList[1]]->dynamic->event.time;
#ifdef DEBUG_IT
            if (timeIncr < 0) {
                printf("!!ERROR!!: time calculation is not correct!\n");
            }
#endif
            UpdateData(timeIncr, "atom", thisThread); //update the coordinates
            TimeForward(timeIncr, "atom", thisThread); //update the node time
            currenttime += timeIncr;
            frame++;

            CommitEvent(thisThread->raw,
                        *newTarget, (threadRenewList[2] > 0 ? *newPartner : NULL),
                        oldTarget, (threadRenewList[2] > 0 ?  oldPartner : NULL),
                        thisThread->listPtr[oldTarget->dynamic->HB.neighbor],
                        (threadRenewList[2] > 0 ? thisThread->listPtr[oldPartner->dynamic->HB.neighbor] : NULL));
            
            SchedulingDelete(threadRenewList, thisThread);
            SchedulingAdd(threadRenewList, thisThread);
            
            processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
            
            if (processratio >= gap + 0.01) {
                
                printf("Process=%8.2lf%%\r", processratio);
                fflush(stdout);
                gap += 0.01;
                
            }
            
            if (strcmp(wallDyn.mark, "no")) DoWallDyn();
            
            //save data to the output files
            if (currenttime >= outputrecord) {
                for (int i = 0; i < lenFileType; i ++) {
                    if (thisThread->fileList[i].mark) {
                        SaveData(i, thisThread);
                    }
                }
                outputrecord += outputrate;
            }
            
            AssignJob(threadRenewList, thisThread);
            AssignThread(threadRenewList, thisThread);
        }
    }
    
    FreeVariables();
    return;
}


void AssignJob(int *renewList, struct ThreadStr *thisThread) {
    //clear renew list
    renewList[0] = 1;
    renewList[1] = 0;
    renewList[2] = 0;
    
    thisThread->atomNum = SchedulingNextEvent(thisThread);
    renewList[1] = thisThread->atomNum;
    
    if (thisThread->raw[renewList[1]]->dynamic->event.partner > 0) {
        renewList[++renewList[0]] = thisThread->raw[renewList[1]]->dynamic->event.partner;
    }
    
    thisThread->getWorkTime = currenttime;
    thisThread->getWorkFrame = frame;
}


void AssignThread(int *renewList, struct ThreadStr *thisThread) {
    int cell_neighbor[4], cellIndex;
    int neighborAtom;
    int count = 0;
    int checkedCell[55] = {0}, num, flag;
    int *sourCellList = celllist;
    list *atomList = &thisThread->atomList;
    listElem *elem = (atomList->num_members == 0) ? NULL : listFirst(atomList);
    struct AtomStr *targetAtom = NULL;
    struct AtomStr **destLibrary = thisThread->listPtr;
    struct AtomStr **sourLibrary = thisThread->raw;
    
    signed int x[27] = {0, 1, 0, 0, 0, 1, 1, 1, -1,  0,  0,  0, -1, -1, -1,  0,  0,  1, -1,  1, -1, -1, -1,  1,  1,  1, -1};
    signed int y[27] = {0, 0, 1, 0, 1, 0, 1, 1,  0, -1,  0, -1,  0, -1, -1, -1,  1,  0,  0, -1,  1, -1,  1, -1,  1, -1,  1};
    signed int z[27] = {0, 0, 0, 1, 1, 1, 0, 1,  0,  0, -1, -1, -1,  0, -1,  1, -1, -1,  1,  0,  0,  1, -1, -1, -1,  1,  1};
    
    for (int n = 1; n <= renewList[0]; n ++) {
        targetAtom = sourLibrary[renewList[n]];
        for (int i = 0; i < 27; ++i) {
            //scan the neighborhood 27 subcells, include the target subcell itself
            cell_neighbor[3] = targetAtom->dynamic->cellIndex[3] + z[i];
            cell_neighbor[2] = targetAtom->dynamic->cellIndex[2] + y[i];
            cell_neighbor[1] = targetAtom->dynamic->cellIndex[1] + x[i];
            
            for (int j = 1; j <= 3; j++) {
                if (cell_neighbor[j] < 0) {
                    cell_neighbor[j] = cellnum[j] - 1;
                } else if (cell_neighbor[j] >= cellnum[j]) {
                    cell_neighbor[j] = 0;
                }
            }
            
            cellIndex = cell_neighbor[3] * cellnum[1] * cellnum[2]
            + cell_neighbor[2] * cellnum[1]
            + cell_neighbor[1] + 1;
            
            neighborAtom = sourCellList[atomnum + cellIndex];
            num = flag = 0;
            if (neighborAtom > 0) {
                while (checkedCell[num]) {
                    if (checkedCell[num] == cellIndex) {
                        flag = 1;
                        break;
                    }
                    num ++;
                }
                if (flag == 0) {
                    checkedCell[num] = cellIndex;
                }
            }
            
            while (neighborAtom != 0 && flag == 0) { //+partner, +check duplicated
                count ++;
                if (count > atomList->num_members && elem == NULL) {
                    struct AtomStr *newAtom = calloc(1, sizeof(struct AtomStr));
                    newAtom->dynamic = (struct DynamicStr *)calloc(1, sizeof(struct DynamicStr));
                    
                    AtomDataCpy(newAtom, sourLibrary[neighborAtom], 1);
                    listAppend(atomList, newAtom);
                    destLibrary[neighborAtom] = newAtom;
                } else {
                    AtomDataCpy(elem->obj, sourLibrary[neighborAtom], 1);
                    destLibrary[neighborAtom] = elem->obj;
                    elem = listNext(atomList, elem);
                }
                
                neighborAtom = sourCellList[neighborAtom];
            }
        }
    }
    
    thisThread->newTarget = destLibrary[renewList[1]];
    AtomDataCpy(&thisThread->oldTarget, sourLibrary[renewList[1]], 2);
    
    if (renewList[2] > 0) {
        thisThread->newPartner = destLibrary[renewList[2]];
        AtomDataCpy(&thisThread->oldPartner, sourLibrary[renewList[2]], 2);
    } else {
        thisThread->newPartner = NULL;
        thisThread->oldPartner.property = NULL;
    }
    
    atomList->num_members = count;
    
#ifdef DEBUG_IT
    if (checkedCell[46] != 0) {
        printf("!!ERROR!!: checkedCell list overflows!\n");
    }
#endif
}


void FreeVariables(void) {
    free(thread[0]->oldTarget.dynamic);
    free(thread[0]->oldTarget.eventList);
    
    free(thread[0]->oldPartner.dynamic);
    free(thread[0]->oldPartner.eventList);
    
    listElem *thisElem = thread[0]->atomList.anchor.next;
    listElem *nextElem;
    while (thisElem != &thread[0]->atomList.anchor) {
        nextElem = thisElem->next;
        if (thisElem->obj) {
            free(((struct AtomStr*)thisElem->obj)->dynamic);
            free(thisElem->obj);
        }
        free(thisElem);
        thisElem = nextElem;
    }
    
    free(thread[0]->raw);
    free(thread[0]->listPtr);
    
    free(thread[0]);
    free(thread);
}
