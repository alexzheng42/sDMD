//
//  FAThread.c
//  Update
//
//  Created by Size Zheng on 8/11/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <sys/time.h>
#include "DMD.h"

int threadOut = 0;

void FAThreadRun(struct ThreadInfoStr *threadInfo);
int CheckStatus(int *renewList, struct ThreadStr *thisThread);


void FAThread(void) {
    thread = (struct ThreadStr **)calloc(threadNum, sizeof(struct ThreadStr*));
    for (int i = 0; i < threadNum; i ++) {
        thread[i] = InitializeThread(i, atom);
        thread[i]->fileList = InitializeFiles(NULL, fileList);
    }
    
    SaveData(lgf, thread[0]);
    if (strncmp("new", neworcontinue, 1) == 0) {
        SaveData(sysInfo, thread[0]);
        FirstRun(thread[0]);
    } else if (strncmp("continue", neworcontinue, 1) == 0) {
        outputrecord = currenttime + outputrate;
        CreateCBT();
    }
    
    pthread_mutex_init(&mstThrLock, NULL);
    
    eventToCommit = CBT.node[1];
    for (int i = 0; i < threadNum; i++) {
        thrInfo[i].threadRenewList[1] = CBT.node[1]; //need to modify
        thread[i]->atomNum = nthNode;
        if (thread[i]->raw[nthNode]->dynamic->event.partner > 0) {
            thrInfo[i].threadRenewList[++thrInfo[i].threadRenewList[0]] = thread[i]->raw[nthNode]->dynamic->event.partner;
        }
        AssignThread(thrInfo[i].threadRenewList, thread[i]);
        
        pthread_create(&thread_t[i], 0, (void *) &FAThreadRun, &thrInfo[i]);
    }
    
    for (int i = 0; i < threadNum; i ++) {
        pthread_join(thread_t[i], 0);
    }
    
    return;
}

void FAThreadRun(struct ThreadInfoStr *threadInfo) {
    int hazardType;
    int missCount = 0;
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
        
        hazardType = ProcessEvent(threadRenewList, thisThread);
        
        if (hazardType < 0) {
            threadRenewList[0] = 1;
            threadRenewList[2] = 0;
            Predict(threadRenewList, thisThread);
        }
        
        //only the thread holding the next event will be processed,
        //otherwise will fall into sleep
        pthread_mutex_lock(&mstThrLock);
        {
            while (threadRenewList[1] != eventToCommit && missCount < 5 && !threadOut) {
                missCount ++;
                pthread_cond_wait(&mstThrQ, &mstThrLock);
            }
        }
        pthread_mutex_unlock(&mstThrLock);
        
        //avoid forever sleeping
        if (missCount >= 5) {
            missCount = 0;
            goto repeat;
        }
        
        if (hazardType < 0) {
#ifdef DETAILS
            printf(", denied, partner changed!\n");
#endif
            pthread_mutex_lock(&mstThrLock);
            AtomDataCpy(thisThread->raw[threadRenewList[1]], thisThread->listPtr[threadRenewList[1]], 0);
            UpdateCBT(threadRenewList);
            pthread_mutex_unlock(&mstThrLock);
            
            countReCal ++;
            
        } else {
            
#ifdef DETAILS
            printf("frame #%li, master thread start to process job. target atom = #%i, partner = #%i, type = %i\n",
                   frame, (*newTarget)->property->num, oldTarget->dynamic->event.partner, oldTarget->dynamic->event.eventType);
            fflush(stdout);
#endif
            
            timeIncr = thisThread->raw[threadRenewList[1]]->dynamic->event.time;
            if (unlikely(timeIncr < 0)) {
                printf("!!ERROR!!: time calculation is not correct!\n");
            }
            
            pthread_mutex_lock(&mstThrLock);
            UpdateData(timeIncr, "atom", thisThread); //update the coordinates
            TimeForward(timeIncr, "atom", thisThread); //update the node time
            currenttime += timeIncr;
            frame++;
            
            CommitEvent(thisThread->raw,
                        *newTarget, (threadRenewList[2] > 0 ? *newPartner : NULL),
                        oldTarget, (threadRenewList[2] > 0 ? oldPartner : NULL),
                        thisThread->listPtr[oldTarget->dynamic->HB.neighbor],
                        thisThread->listPtr[oldPartner->dynamic->HB.neighbor]);
            UpdateCBT(threadRenewList);
            pthread_mutex_unlock(&mstThrLock);
            
            if (strcmp(wallDyn.mark, "no")) DoWallDyn();
            
            processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
            if (unlikely(processratio >= gap + 0.01)) {
                
                printf("Process=%8.2lf%%\r", processratio);
                fflush(stdout);
                gap += 0.01;
                
            }
            
            //save data to the output files
            if (unlikely(currenttime >= outputrecord)) {
                for (int i = 0; i < lenFileType; i ++) {
                    if (thisThread->fileList[i].mark) {
                        SaveData(i, thisThread);
                    }
                }
                outputrecord += outputrate;
            }
        }
        eventToCommit = CBT.node[1];
        
    repeat:
        threadRenewList[0] = 1;
        threadRenewList[1] = threadRenewList[2] = 0;
        
        pthread_mutex_lock(&mstThrLock);
        pthread_cond_broadcast(&mstThrQ);
        pthread_mutex_unlock(&mstThrLock);
        
        //waiting to be called
        //main frame lock
        pthread_mutex_lock(&mstThrLock);
        {
            while (CheckStatus(threadRenewList, thisThread) && !threadOut) {
                threadRenewList[0] = 1;
                threadRenewList[1] = threadRenewList[2] = 0;
                pthread_cond_wait(&mstThrQ, &mstThrLock);
            }
        }
        AssignThread(threadRenewList, thisThread);
        pthread_mutex_unlock(&mstThrLock);
    }

    pthread_mutex_lock(&mstThrLock);
    threadOut = -1;
    
    printf("!WARNING!: master thread out!\n");
    
    pthread_cond_broadcast(&mstThrQ);
    pthread_mutex_unlock(&mstThrLock);

    return;
}


int CheckStatus(int *renewList, struct ThreadStr *thisThread) {
    int cell_neighbor[4], cellIndex;
    int neighborAtom;
    int checkedCell[55] = {0}, num, flag;
    int *sourCellList = celllist;
    struct AtomStr *targetAtom = NULL;
    struct AtomStr **sourLibrary = thisThread->raw;
    
    signed int x[27] = {0, 1, 0, 0, 0, 1, 1, 1, -1,  0,  0,  0, -1, -1, -1,  0,  0,  1, -1,  1, -1, -1, -1,  1,  1,  1, -1};
    signed int y[27] = {0, 0, 1, 0, 1, 0, 1, 1,  0, -1,  0, -1,  0, -1, -1, -1,  1,  0,  0, -1,  1, -1,  1, -1,  1, -1,  1};
    signed int z[27] = {0, 0, 0, 1, 1, 1, 0, 1,  0,  0, -1, -1, -1,  0, -1,  1, -1, -1,  1,  0,  0,  1, -1, -1, -1,  1,  1};
    
    //need to modify
    if (nthCheck > threadNum * 2) {
        return 1;
    }

    thisThread->atomNum = nthNode;
    renewList[1] = nthNode;
    
    if (thisThread->raw[nthNode]->dynamic->event.partner > 0) {
        renewList[++renewList[0]] = thisThread->raw[nthNode]->dynamic->event.partner;
    }
    
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
                if (neighborAtom != renewList[n]) {
                    for (int j = 0; j < threadNum; j ++) {
                        if (j != thisThread->tid)
                            if ((thrInfo[j].threadRenewList[1] == neighborAtom || thrInfo[j].threadRenewList[2] == neighborAtom) &&
                                atom[thrInfo[thisThread->tid].threadRenewList[1]].dynamic->event.time > atom[neighborAtom].dynamic->event.time)
                            return 1;
                    }
                }
                neighborAtom = sourCellList[neighborAtom];
            }
        }        
    }
    
    thisThread->finishWork = 0;
    return 0;
}
