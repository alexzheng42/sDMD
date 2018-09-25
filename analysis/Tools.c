//
//  Tools.c
//  Analysis
//
//  Created by Size Zheng on 10/31/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"


double absvalue_vector(double * vector)
{
    double v12, v22 ,v32;
    
    v12=vector[1]*vector[1];
    v22=vector[2]*vector[2];
    v32=vector[3]*vector[3];
    
    return sqrt(v12+v22+v32);
}


void EstCell(void) {
    for (int i = 1; i <= 3; i++) {
        cellNum[i] = floor(boxOrigDim[i] / cutoffr);
        
        if (cellNum[i] < 3) {
            cellNum[i] = 3;
        }
    }
    cellNum[0] = atomnum + cellNum[1] * cellNum[2] * cellNum[3] + 1;
    
    for (int i = 1; i <= 3; i++) {
        cellSize[i] = boxOrigDim[i] / cellNum[i];
    }
    
    if (celllist) {
        free(celllist);
    }
    celllist = (int *)calloc((atomnum + cellNum[1] * cellNum[2] * cellNum[3] + 1), sizeof (int));
}


void LinkList(void) {
    int i, n;
    int cellindex;
    double *boxsize;
    double *position;
    
    pointToStruct(boxsize, boxOrigDim);
    for (int i = 0; i < cellNum[0]; i ++) {
        celllist[i] = 0;
    }
    
    //-------------------------------------------
    //calculate which subcells the atoms are in and generate the list
    for (i = 1; i <= atomnum; i++) {
        
        pointToStruct(position, atom[i].dynamic->coordinate);
        
        for (n = 1; n <= 3; n++) {
            
            atom[i].dynamic->cellIndex[n] = position[n] / cellSize[n]; //include 0
            //atom on the edge of box size has been removed by function pbc()
            //integer numbers
            //for a 3*3*3 system, cell index on one side would be:
            //  [0)    [1)   [2)
            //  [0)    [1)   [2)
            //  [0)    [1)   [2)
            
            //just in case the floating-point arithmetic error
            //the PBC-adjusted coordinates may be exactly equal to the upper edge of the PBC box
            if (atom[i].dynamic->cellIndex[n] == cellNum[n]) {
                atom[i].dynamic->cellIndex[n]--;
            }
        }
        
        cellindex = atom[i].dynamic->cellIndex[3] * cellNum[1] * cellNum[2]
        + atom[i].dynamic->cellIndex[2] * cellNum[1]
        + atom[i].dynamic->cellIndex[1] + 1;
        
        celllist[i] = celllist[atomnum + cellindex];
        celllist[atomnum + cellindex] = i;
    }
}

void PrintProcess(long step) {
    int n = step % 4;
    
    switch (n) {
        case 0:
            printf("\b-");
            break;
        case 1:
            printf("\b\\");
            break;
        case 2:
            printf("\b|");
            break;
        case 3:
            printf("\b/");
            break;
            
        default:
            break;
    }
    
    fflush(stdout);
}

int AtomModel(char *type) {
    int num = 0;
    
    if (strcmp(type, "OZB") == 0) {
        num = 26;
    } else if (strcmp(type, "NZB") == 0) {
        num = 18;
    } else if (strcmp(type, "CRB") == 0) {
        num = 6;
    } else if (strcmp(type, "CA") == 0) {
        num = 9;
    } else if (strcmp(type, "CA1") == 0) {
        num = 10;
    } else if (strcmp(type, "CA2") == 0) {
        num = 11;
    } else if (strcmp(type, "CA3") == 0) {
        num = 12;
    } else if (strcmp(type, "CA1P") == 0) {
        num = 14;
    } else if (strcmp(type, "CA2P") == 0) {
        num = 15;
    } else if (strcmp(type, "CM1") == 0) {
        num = 13;
    } else if (strcmp(type, "CM1P") == 0) {
        num = 16;
    } else if (strcmp(type, "CR") == 0) {
        num = 5;
    } else if (strcmp(type, "CRN") == 0) {
        num = 7;
    } else if (strcmp(type, "CRNP") == 0) {
        num = 8;
    } else if (strcmp(type, "NCK") == 0) {
        num = 21;
    } else if (strcmp(type, "NCR") == 0) {
        num = 22;
    } else if (strcmp(type, "NR") == 0) {
        num = 23;
    } else if (strcmp(type, "NZNP") == 0) {
        num = 19;
    } else if (strcmp(type, "NZNQ") == 0) {
        num = 20;
    } else if (strcmp(type, "NZ") == 0) {
        num = 17;
    } else if (strcmp(type, "OC") == 0) {
        num = 28;
    } else if (strcmp(type, "OW") == 0) {
        num = 29;
    } else if (strcmp(type, "OZH") == 0) {
        num = 25;
    } else if (strcmp(type, "OZ") == 0) {
        num = 24;
    } else if (strcmp(type, "SGNP") == 0) {
        num = 30;
    } else if (strcmp(type, "SG") == 0) {
        num = 31;
    } else if (strcmp(type, "OZB_HB") == 0) {
        num = 27;
    } else if (strcmp(type, "HB") == 0) {
        num = 2;
    } else if (strcmp(type, "HO") == 0) {
        num = 1;
    } else if (strcmp(type, "HN") == 0) {
        num = 1;
    } else if (strcmp(type, "HS") == 0) { //hydrogen on sulfur
        num = 1;
    } else if (strcmp(type, "HW") == 0) {
        num = 1;
    } else if (strcmp(type, "H") == 0) {
        num = 1;
    } else if (strcmp(type, "H_HB") == 0) {
        num = 3;
    } else if (strcmp(type, "HB_HB") == 0) {
        num = 4;
    } else {
        printf("\nNo atom type is matched! Check the atom outerType assignment! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return num;
}
