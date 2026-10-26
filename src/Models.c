//
//  Models.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"


int FindResIndex(char *name) { 
    for (int num = 1; num <= FFInfo.group[0].val[0]; num ++)
        if (strcmp(name, FFInfo.group[num].name) == 0)
            return FFInfo.group[num].val[0];

    printf("\n!ERROR!:The amino acid name is not found in Library! Check the input coordinate file! %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


int FindAtomIndex(char *name) {
    for (int num = 1; num <= FFInfo.bead[0].val[0]; num ++)
        if (strcmp(name, FFInfo.bead[num].name) == 0)
            return FFInfo.bead[num].val[0];
    
    printf("\nNo atom type is matched to %s! Check the atom outerType assignment! %s:%i\n", name, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


double HBBarrier(struct AtomStr *atom1, struct AtomStr *atom2) {
    int type = HBModel(atom1, atom2);
    return (type > 0) ? HBPotential.barrier[HBPotential.classOf[type]] : 0;
}

int HBModel(struct AtomStr *atom1, struct AtomStr *atom2) {
    int atom1_OrigType = atom1->property->typeofAtom[1];
    int atom2_OrigType = atom2->property->typeofAtom[1];
    
    for (int i = 1; i <= HBPotential.typeNum; i ++) {
        if (potentialPairHB[i][atom1_OrigType][atom2_OrigType].dmin &&
            potentialPairHB[i][atom1_OrigType][atom2_OrigType].step) {
            return i;
        }
    }
    
    if (strcmp(Methodtype, "Ding") == 0 &&
        ((atom1_OrigType == 2 && (atom2->property->typeofRes == 3 || atom2->property->typeofRes == 6)) ||
         (atom2_OrigType == 2 && (atom1->property->typeofRes == 3 || atom1->property->typeofRes == 6)))
        )
        return -1; //no HB between ASN/GLN and backbone H for FF_Ding
    
    return -1;
}


int NeighborModel(char *atom1, char *atom2) {
    
    if (strncmp(atom1, "C", 1) == 0) {
        return 3;
    } else if (strcmp(atom1, "OW") == 0) {
        return 5;
    } else {
        if (strncmp(atom1, "N", 1) == 0) {
            if (strncmp(atom2, "O", 1) == 0) {
                return 2;
            } else {
                return 6;
            }
        } else {
            if (strncmp(atom2, "O", 1) == 0) {
                return 5;
            } else {
                return 2;
            }
        }
    }
    
    return -1;
}

int AAList(char *name) {
    /*
     This function is to distinguish single molecules and
     amino acids in proteins, such that the program could
     know where a new molecule begins while reading coordinate
     file, in order to properly simulate hydrogen bond formation,
     since HB only forms between amino acids separated by 4 if
     they are in the same protein, but can form freely if they
     are in different proteins. See TimePrediction.c : 424.
     */
    char AA[20][3] = {
                      "ALA",
                      "ARG",
                      "ASN",
                      "ASP",
                      "CYS",
                      "GLN",
                      "GLU",
                      "GLY",
                      "HIS",
                      "ILE",
                      "LEU",
                      "LYS",
                      "MET",
                      "PHE",
                      "PRO",
                      "SER",
                      "THR",
                      "TRP",
                      "TYR",
                      "VAL"
                     };
    
    for (int i = 0; i < 20; i ++) {
        if (strncmp(name, AA[i], 3) == 0)
            return 1;
    }

    return 0;
}

int SingleParticleList(char *name) {
    char type[20][3] = {"PHO", "PHI", "RHI"};
    
    for (int i = 0; i < 20; i ++) {
        if (strncmp(name, type[i], 3) == 0)
            return 1;
    }
    
    return 0;
}
