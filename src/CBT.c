//
//  CBT.c
//  sDMD
//
//  Created by Size Zheng on 5/22/18.
//  Copyright © 2018 Size Zheng. All rights reserved.
//

#include "DMD.h"

void UpdateCBTRun(int atomNum);
void InsertCBTRun(int atomNum);
void DeleteCBTRun(int atomNum);


void CreateCBT(void) {
    CBT.node = (int *)calloc(atomnum * 2, sizeof(int));
    CBT.leaf = (int *)calloc(atomnum + 1, sizeof(int));
    CBT.time = (double **)calloc(atomnum + 1, sizeof(double *));
    
    for (int i = 1; i <= atomnum; i ++) {
        CBT.time[i] = &atom[i].dynamic->event.time;
    }
    
    for (int i = 1; i <= atomnum; i ++) {
        InsertCBTRun(i);
    }
    
    return;
}


void UpdateCBT(int *renewList) {
    for (int i = 1; i <= renewList[0]; i ++) {
        UpdateCBTRun(renewList[i]);
    }
    
    return;
}


void InsertCBT(int *renewList) {
    for (int i = 1; i <= renewList[0]; i ++) {
        InsertCBTRun(renewList[i]);
    }
}


void DeleteCBT(int *renewList) {
    for (int i = 1; i <= renewList[0]; i ++) {
        DeleteCBTRun(renewList[i]);
    }
}


void UpdateCBTRun(int atomNum) {
    int father, left, right, old;
    
    for (father = CBT.leaf[atomNum] / 2; father > 0; father /= 2) {
        if (CBT.node[father] != atomNum) break; //jumps to the next "for";
        
        left  = CBT.node[father * 2];
        right = CBT.node[father * 2 + 1];
        if (*CBT.time[left] < *CBT.time[right])
            CBT.node[father] = left;
        else
            CBT.node[father] = right;
    }
    
    for (; father > 0; father /= 2) {
        old   = CBT.node[father];
        left  = CBT.node[father * 2];
        right = CBT.node[father * 2 + 1];
        if (*CBT.time[left] < *CBT.time[right])
            CBT.node[father] = left;
        else
            CBT.node[father] = right;
        
        if (CBT.node[father] == old) return;
    }
    
    return;
}


void InsertCBTRun(int atomNum) {
    int i;
    
    if (CBT.count == 0) {
        CBT.node[1] = atomNum;
        CBT.count ++;
        return;
    }
    
    i = CBT.node[CBT.count];
    CBT.node[CBT.count * 2] = i;
    CBT.node[CBT.count * 2 + 1] = atomNum;
    CBT.leaf[i] = CBT.count * 2;
    CBT.leaf[atomNum] = CBT.count * 2 + 1;
    
    CBT.count ++;
    UpdateCBTRun(i);
    
    return;
}


void DeleteCBTRun(int atomNum) {
    int i;
    
    if (CBT.count < 2) {
        CBT.node[1] = 0;
        CBT.leaf[0] = 1;
        CBT.count --;
        return;
    }
    
    i = CBT.count * 2 - 1;
    if (CBT.node[i - 1] != atomNum) {
        CBT.leaf[CBT.node[i - 1]] = i / 2;
        CBT.node[i / 2] = CBT.node[i - 1];
        UpdateCBTRun(CBT.node[i - 1]);
    } else {
        CBT.leaf[CBT.node[i]] = i / 2;
        CBT.node[i / 2] = CBT.node[i];
        UpdateCBTRun(CBT.node[i]);
        CBT.count --;
        return;
    }
    
    if (CBT.node[i] != atomNum) {
        CBT.node[CBT.leaf[atomNum]] = CBT.node[i];
        CBT.leaf[CBT.node[i]] = CBT.leaf[atomNum];
        UpdateCBTRun(CBT.node[i]);
    }
    
    CBT.count --;
    return;
}

