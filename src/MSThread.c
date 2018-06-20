//
//  MSThread.c
//  Update
//
//  Created by Size Zheng on 8/11/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <sys/time.h>
#include "DMD.h"

int masterOut = 0;
long missCount = 0;

#ifdef DETAILS
struct timeval zeroClock;
double FindTime(void);
#endif

void MasterThreadRun(struct ThreadInfoStr *threadInfo);
void SlaveThreadRun(struct ThreadInfoStr *threadInfo);
void ResetPreCalList(int thisAtom, struct PreCalObjStr *thisList);
void InitializePreCalList(struct PreCalObjStr *thisList);
void CreateAtomStr(struct AtomStr *thisAtom);
void MSAssignJob(int *renewList);
int CheckConflictAndAssign(struct PreCalObjStr **thisList);


void MSThread(void) {
    thread = (struct ThreadStr **)calloc(1, sizeof(struct ThreadStr*));
    thread[0] = InitializeThread(0, atom);
    thread[0]->fileList = InitializeFiles(NULL, fileList);
    
    preCalList = (struct PreCalObjStr *)calloc(atomnum + 1, sizeof(struct PreCalObjStr));
    InitializePreCalList(preCalList);
    
    SaveData(lgf, thread[0]);
    if (strncmp("new", neworcontinue, 1) == 0) {
        SaveData(sysInfo, thread[0]);
        FirstRun(thread[0]);
    } else if (strncmp("continue", neworcontinue, 1) == 0) {
        outputrecord = currenttime + outputrate;
        CreateCBT();
    }
    
    pthread_mutex_init(&mstThrLock, NULL);
    pthread_mutex_init(&slvThrLock, NULL);
    
#ifdef DETAILS
    gettimeofday(&zeroClock, 0);
#endif
    
    pthread_create(&thread_t[0], 0, (void *) &MasterThreadRun, &thrInfo[0]);
    for (int i = 1; i < threadNum; i++) {
        pthread_create(&thread_t[i], 0, (void *) &SlaveThreadRun, &thrInfo[i]);
    }
    
    for (int i = 0; i < threadNum; i ++) {
        pthread_join(thread_t[i], 0);
    }
    
    return;
}

void MasterThreadRun(struct ThreadInfoStr *threadInfo) {
    int hazardType;
    int targetNum;
    int *threadRenewList;
    double timeIncr;
    double processratio, gap = 0;
    struct ThreadStr *thisThread;
    struct PreCalObjStr *thisCalObj;
    struct AtomStr *newTarget;
    struct AtomStr *newPartner;
    struct AtomStr *oldTarget;
    struct AtomStr *oldPartner;
    struct AtomStr *neighbor_i = NULL, *neighbor_j = NULL;
    
    while (currenttime <= timestep) {
        
        eventToCommit = CBT.node[1]; //find the nearest event
        pthread_mutex_lock(&mstThrLock);
        {
            while (preCalList[eventToCommit].eventStatus != 1) { //if the target atom has not been processed yet, wait
#ifdef DETAILS
                printf("%012.3lfms, frame #%li, event to commit, atom #%i, is not available. master thread is sleeping\n", FindTime(), frame, eventToCommit);
                fflush(stdout);
#endif
                pthread_cond_wait(&mstThrQ, &mstThrLock);
                missCount ++;
                
#ifdef DETAILS
                printf("%012.3lfms, frame #%li, master thread is waking up\n", FindTime(), frame);
                fflush(stdout);
#endif
            }
            
            preCalList[eventToCommit].eventStatus = 2; //changed the status code, so the other threads will not use this atom any longer
        }
        pthread_mutex_unlock(&mstThrLock);
        
#ifdef DETAILS
        printf("%012.3lfms, frame #%li, event to commit, atom #%i, has been available.\n", FindTime(), frame, eventToCommit);
        fflush(stdout);
#endif
        
        targetNum       = eventToCommit;
        thisCalObj      = &preCalList[targetNum];
        thisThread      = thisCalObj->data; //point to the pre-calculated list
        threadRenewList = thisCalObj->renewList;

        newTarget  = thisThread->newTarget;
        oldTarget  = &thisThread->oldTarget;
        neighbor_i = thisThread->listPtr[oldTarget->dynamic->HB.neighbor];

        if (threadRenewList[2] > 0) {
            newPartner = thisThread->newPartner;
            oldPartner = &thisThread->oldPartner;
            neighbor_j = thisThread->listPtr[oldPartner->dynamic->HB.neighbor];
        } else {
            newPartner = NULL;
            oldPartner = NULL;
            neighbor_j = NULL;
        }
        
#ifdef DETAILS
        printf("%012.3lfms, frame #%li, master thread start to process job. target atom = #%i, partner = #%i, type = %i\n",
               FindTime(), frame, newTarget->property->num, oldTarget->dynamic->event.partner, oldTarget->dynamic->event.eventType);
        fflush(stdout);
#endif
        
        hazardType = HazardCheck(oldTarget, (oldPartner != NULL ? oldPartner : NULL), neighbor_i, neighbor_j, thisThread);
        
        if (hazardType < 0) {
#ifdef DETAILS
            printf("%012.3lfms, frame #%li, atom #%i, hazard code 2, recalculating\n", FindTime(), frame, eventToCommit);
            fflush(stdout);
#endif
            
            threadRenewList[0] = 1;
            threadRenewList[2] = 0;
            
            AssignThread(threadRenewList, thisThread);
            Predict(threadRenewList, thisThread);
            
            pthread_mutex_lock(&slvThrLock);
            {
                AtomDataCpy(thisThread->raw[targetNum], thisThread->listPtr[targetNum], 0);
                UpdateCBT(threadRenewList);
                ResetPreCalList(targetNum, thisCalObj);
                
#ifdef DETAILS
                printf("%012.3lfms, frame #%li, master thread finished recalculating, call slv thread\n", FindTime(), frame);
                fflush(stdout);
#endif
            }
            pthread_cond_signal(&slvThrQ);
            pthread_mutex_unlock(&slvThrLock);
            
            countReCal ++;
            
        } else {
#ifdef DEBUG_PRINTF
            PrintData(threadRenewList);
#endif
#ifdef DETAILS
            printf("%012.3lfms, frame #%li, atom #%i, commited\n", FindTime(), frame, eventToCommit);
            fflush(stdout);
#endif
            
            timeIncr = thisThread->raw[targetNum]->dynamic->event.time;
            if (unlikely(timeIncr < 0)) {
                printf("!!ERROR!!: time calculation is not correct! %s:%i\n", __FILE__, __LINE__);
            }
            
            pthread_mutex_lock(&slvThrLock);
            {
                TimeForward(timeIncr, thisThread); //update the node time
                UpdateData(timeIncr, "atom", thisThread); //update the coordinates
                currenttime += timeIncr;
                frame++;
                
                CommitEvent(thisThread->raw, newTarget, (threadRenewList[2] > 0 ? newPartner : NULL), oldTarget, (threadRenewList[2] > 0 ? oldPartner : NULL), neighbor_i, neighbor_j);
                UpdateCBT(threadRenewList);
                ResetPreCalList(targetNum, thisCalObj);
                
#ifdef DETAILS
                printf("%012.3lfms, frame #%li, master thread finished submitting, call slv thread\n", FindTime(), frame);
                fflush(stdout);
#endif
            }
            pthread_cond_signal(&slvThrQ);
            pthread_mutex_unlock(&slvThrLock);
            
            processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
            
            if (processratio >= gap + 0.01) {
                printf("Process=%8.2lf%%\r", processratio);
                fflush(stdout);
                gap += 0.01;
            }
            
            //save data to the output files
            if (currenttime >= outputrecord) {
                for (int i = 0; i < lenFileType; i ++) {
                    if (thread[0]->fileList[i].mark) {
                        SaveData(i, thread[0]);
                    }
                }
                outputrecord += outputrate;
            }
        }
        
        /*
        int i = 0;
        while (i < 10000000) {
            i ++;
        }
         */
    }
    
    pthread_mutex_lock(&slvThrLock);
    masterOut = -1;

    printf("miss count = %li\n", missCount);
    printf("!WARNING!: master thread out!\n");
    
    pthread_cond_broadcast(&slvThrQ);
    pthread_mutex_unlock(&slvThrLock);

    return;
}

void SlaveThreadRun(struct ThreadInfoStr *threadInfo) {
    int tid = threadInfo->threadID;
    int *threadRenewList;
    double timeIncr;
    struct AtomStr *targetAtom, *partner;
    struct AtomStr *oldTarget, *oldPartner;
    struct ThreadStr *thisThread;
    struct PreCalObjStr *thisList;
    
    while (currenttime <= timestep) {
        
        pthread_mutex_lock(&slvThrLock);
        {
            while ((CheckConflictAndAssign(&thisList)) && !masterOut) {
#ifdef DETAILS
                printf("%012.3lfms, frame #%li, slv thread #%i is sleeping\n", FindTime(), frame, tid);
                fflush(stdout);
#endif
                pthread_cond_wait(&slvThrQ, &slvThrLock);
#ifdef DETAILS
                printf("%012.3lfms, frame #%li, slv thread #%i is waking up\n", FindTime(), frame, tid);
                fflush(stdout);
#endif
            }
#ifdef DETAILS
            printf("%012.3lfms, frame #%li, slv thread #%i start to process atom #%i\n", FindTime(), frame, tid, thisList->renewList[1]);
            fflush(stdout);
#endif
        }
        pthread_mutex_unlock(&slvThrLock);
        
        if (unlikely(masterOut == -1)) {
            goto finish;
        }
        
        thisThread          = thisList->data;
        threadRenewList     = thisList->renewList;
        thisThread->atomNum = threadRenewList[1];
        
        timeIncr = thisThread->newTarget->dynamic->event.time;
        UpdateData(timeIncr, "partial", thisThread);
        
        targetAtom = thisThread->newTarget;
        partner    = thisThread->newPartner;
        oldTarget  = &thisThread->oldTarget;
        oldPartner = &thisThread->oldPartner;
        
        oldTarget->dynamic->HB.interactionType = DoEvent(thisThread); //calculate the after-event velocities
        
        targetAtom->dynamic->event.counter ++;
        if (threadRenewList[2] > 0) {
            partner->dynamic->event.counter ++;
        }
        
        Predict(threadRenewList, thisThread); //pre-predict the after-event event
        thisList->eventStatus = 1;
        
        pthread_mutex_lock(&mstThrLock);
        {
#ifdef DETAILS
            printf("%012.3lfms, frame #%li, slv thread #%i finished processing atom #%i, call master thread\n", FindTime(), frame, tid, thisThread->atomNum);
            fflush(stdout);
#endif
            pthread_cond_signal(&mstThrQ);
        }
        pthread_mutex_unlock(&mstThrLock);
    }
    
finish:
    pthread_mutex_lock(&slvThrLock);
    
    printf("!WARNING!: slave thread #%i out!\n", tid);
    
    pthread_cond_broadcast(&slvThrQ);
    pthread_mutex_unlock(&slvThrLock);
    
    return;
}


void MSAssignJob(int *renewList) {
    //clear renew list
    renewList[0] = 1;
    renewList[1] = 0;
    renewList[2] = 0;
    
    renewList[1] = CBT.node[1]; //need to modify
    
    if (atom[renewList[1]].dynamic->event.partner > 0) {
        renewList[++renewList[0]] = atom[renewList[1]].dynamic->event.partner;
    }
}


int CheckConflictAndAssign(struct PreCalObjStr **thisList) {
    int cell_neighbor[4], cellIndex;
    int neighborAtom;
    int count = 0;
    int checkedCell[55] = {0}, num, flag;
    int *sourCellList = celllist;
    int renewList[4] = {0};
    struct ThreadStr *thisThread;
    struct AtomStr *targetAtom = NULL;
    struct AtomStr **destLibrary;
    struct AtomStr **sourLibrary;
    list *atomList;
    listElem *elem;
    
    signed int x[27] = {0, 1, 0, 0, 0, 1, 1, 1, -1,  0,  0,  0, -1, -1, -1,  0,  0,  1, -1,  1, -1, -1, -1,  1,  1,  1, -1};
    signed int y[27] = {0, 0, 1, 0, 1, 0, 1, 1,  0, -1,  0, -1,  0, -1, -1, -1,  1,  0,  0, -1,  1, -1,  1, -1,  1, -1,  1};
    signed int z[27] = {0, 0, 0, 1, 1, 1, 0, 1,  0,  0, -1, -1, -1,  0, -1,  1, -1, -1,  1,  0,  0,  1, -1, -1, -1,  1,  1};
    
    MSAssignJob(renewList);
    
    if (nthCheck > 5) {
        return 1;
    }
    
    *thisList = &preCalList[renewList[1]];
    (*thisList)->eventStatus = 2;
    
    (*thisList)->renewList[0] = renewList[0];
    (*thisList)->renewList[1] = renewList[1];
    (*thisList)->renewList[2] = renewList[2];
    
    thisThread = (*thisList)->data;
    destLibrary = thisThread->listPtr;
    sourLibrary = thisThread->raw;
    atomList = &thisThread->atomList;
    elem = (atomList->num_members == 0) ? NULL : listFirst(atomList);

    for (int n = 1; n <= renewList[0]; n ++) {
        targetAtom = sourLibrary[renewList[n]];
        for (int i = 0; i < 27; i ++) {
            cell_neighbor[3] = targetAtom->dynamic->cellIndex[3] + z[i];
            cell_neighbor[2] = targetAtom->dynamic->cellIndex[2] + y[i];
            cell_neighbor[1] = targetAtom->dynamic->cellIndex[1] + x[i];
            
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
            
            while (neighborAtom != 0 && flag == 0) {
                if (neighborAtom != renewList[n] &&
                    preCalList[neighborAtom].eventStatus > 0 &&
                    atom[neighborAtom].dynamic->event.time < atom[renewList[1]].dynamic->event.time) {
                    (*thisList)->eventStatus = 0;
                    return 1;
                }
                neighborAtom = sourCellList[neighborAtom];
            }
        }        
    }
    
    for (int i = 0; i < 55 && checkedCell[i]; i ++) {
        neighborAtom = sourCellList[atomnum + checkedCell[i]];
        while (neighborAtom != 0) {
            count ++;
            if (unlikely(count > atomList->num_members && elem == NULL)) {
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
    atomList->num_members = count;

    thisThread->newTarget = destLibrary[renewList[1]];
    AtomDataCpy(&thisThread->oldTarget, sourLibrary[renewList[1]], 2);
    
    if (renewList[2] > 0) {
        thisThread->newPartner = destLibrary[renewList[2]];
        AtomDataCpy(&thisThread->oldPartner, sourLibrary[renewList[2]], 2);
    } else {
        thisThread->newPartner = NULL;
        thisThread->oldPartner.property = NULL;
    }
    return 0;
}


void ResetPreCalList(int thisAtom, struct PreCalObjStr *thisList) {
    thisList->eventStatus = 0;
    
    thisList->renewList[0] = 1;
    thisList->renewList[1] = 0;
    thisList->renewList[2] = 0;
    thisList->renewList[3] = 0;
}

void InitializePreCalList(struct PreCalObjStr *thisList) {
    for (int i = 0; i <= atomnum; i ++) {
        ResetPreCalList(i, &thisList[i]);
        thisList[i].data = InitializeThread(-1, atom);
    }
}

void CreateAtomStr(struct AtomStr *thisAtom) {
    thisAtom->property = (struct PropertyStr *)calloc(1, sizeof(struct PropertyStr));
    thisAtom->dynamic = (struct DynamicStr *)calloc(1, sizeof(struct DynamicStr));
}


#ifdef DETAILS
double FindTime(void) {
    struct timeval clock;
    double curTime;
    
    gettimeofday(&clock, 0);
    curTime = (clock.tv_sec - zeroClock.tv_sec) * 1000.0 + (clock.tv_usec - zeroClock.tv_usec) / 1000.0;
    
    return curTime;
}
#endif
