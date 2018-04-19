//
//  Models.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

int AAModel(char *component) {
    int num = 0;
    
    if (strcmp(component, "ALA") == 0) {
        num = 6; //bead num
    } else if (strcmp(component, "ARG") == 0) {
        num = 17; //bead num
    } else if (strcmp(component, "ASN") == 0) {
        num = 11; //bead num
    } else if (strcmp(component, "ASP") == 0) {
        num = 9; //bead num
    } else if (strcmp(component, "CYS") == 0) {
        num = 8; //bead num
    } else if (strcmp(component, "GLN") == 0) {
        num = 12; //bead num
    } else if (strcmp(component, "GLU") == 0) {
        num = 10; //bead num
    } else if (strcmp(component, "GLY") == 0) {
        num = 5; //bead num
    } else if (strcmp(component, "HIS") == 0) {
        num = 12; //bead num
    } else if (strcmp(component, "ILE") == 0) {
        num = 9; //bead num
    } else if (strcmp(component, "LEU") == 0) {
        num = 9; //bead num
    } else if (strcmp(component, "LYS") == 0) {
        num = 13; //bead num
    } else if (strcmp(component, "MET") == 0) {
        num = 9; //bead num
    } else if (strcmp(component, "PHE") == 0) {
        num = 12; //bead num
    } else if (strcmp(component, "PRO") == 0) {
        num = 7; //bead num
    } else if (strcmp(component, "SER") == 0) {
        num = 8; //bead num
    } else if (strcmp(component, "THR") == 0) {
        num = 9; //bead num
    } else if (strcmp(component, "TRP") == 0) {
        num = 16; //bead num
    } else if (strcmp(component, "TYR") == 0) {
        num = 14; //bead num
    } else if (strcmp(component, "VAL") == 0) {
        num = 8; //bead num
    } else if (strcmp(component, "SOL") == 0) {
        num = 3; //water
    } else if (strcmp(component, "GEL") == 0) {
        num = 1; //GEL ball
    } else {
        printf("\n!ERROR!:The amino acid name is not found in Library! Check the input coordinate file!\n");
        exit(EXIT_FAILURE);
    }
    
    return num;
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


int HBModel(struct AtomStr *atom1, struct AtomStr *atom2) {
    struct AtomStr *switcher;
    
    if (atom1->dynamic->HB.role == 'A') { //atom1 -> donator, atom2 -> acceptor
        switcher = atom1;
        atom1 = atom2;
        atom2 = switcher;
    }
    
    if (AtomTypeChange(atom1->property->type, 0) == 2 /*HB*/ &&
        AtomTypeChange(atom2->property->type, 0) == 26 /*OZB*/) {
        return 1;
    } else if (AtomTypeChange(atom1->property->type, 0) == 2) {
        if (strcmp(atom2->property->nameOfAA, "ASP") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLU") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLN") == 0 ||
            strcmp(atom2->property->nameOfAA, "ASN") == 0) {
            return 2;
        } else if (strcmp(atom2->property->nameOfAA, "HIS") == 0 ||
                   strcmp(atom2->property->nameOfAA, "TRP") == 0) {
            return 4;
        } else if (AtomTypeChange(atom2->property->type, 0) == 31) {
            return 10;
        } else {
            return 3;
        }
    } else if (AtomTypeChange(atom2->property->type, 0) == 26) {
        return 5;
    } else if (AtomTypeChange(atom1->property->type, 0) == 31 &&
               AtomTypeChange(atom2->property->type, 0) == 31) {
        return 9;
    } else {
        if (strcmp(atom2->property->nameOfAA, "ASP") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLU") == 0 ||
            strcmp(atom2->property->nameOfAA, "GLN") == 0 ||
            strcmp(atom2->property->nameOfAA, "ASN") == 0) {
            return 6;
        } else if (strcmp(atom2->property->nameOfAA, "HIS") == 0 ||
                   strcmp(atom2->property->nameOfAA, "TRP") == 0) {
            return 8;
        } else {
            return 7;
        }
    }
    
    return -1;
}


int NeighborModel(char *atom1, char *atom2) {
    int num = 0;
    
    if (strncmp(atom1, "C", 1) == 0) {
        num = 3;
    } else if (strcmp(atom1, "OW") == 0) {
        num = 5;
    } else {
        if (strncmp(atom1, "N", 1) == 0) {
            if (strncmp(atom2, "O", 1) == 0) {
                num = 2;
            } else {
                num = 6;
            }
        } else {
            if (strncmp(atom2, "O", 1) == 0) {
                num = 5;
            } else {
                num = 2;
            }
        }
    }
    
    return num;
}

