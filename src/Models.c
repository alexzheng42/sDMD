//
//  Models.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

int FindNumInAA(char *component) { //find the number of atoms containing in the AA
    if (strcmp(component, "ALA") == 0) {
        return 6;
    } else if (strcmp(component, "ARG") == 0) {
        return 17;
    } else if (strcmp(component, "ASN") == 0) {
        return 11;
    } else if (strcmp(component, "ASP") == 0) {
        return 9;
    } else if (strcmp(component, "CYS") == 0) {
        return 8;
    } else if (strcmp(component, "GLN") == 0) {
        return 12;
    } else if (strcmp(component, "GLU") == 0) {
        return 10;
    } else if (strcmp(component, "GLY") == 0) {
        return 5;
    } else if (strcmp(component, "HIS") == 0) {
        return 12;
    } else if (strcmp(component, "ILE") == 0) {
        return 9;
    } else if (strcmp(component, "LEU") == 0) {
        return 9;
    } else if (strcmp(component, "LYS") == 0) {
        return 13;
    } else if (strcmp(component, "MET") == 0) {
        return 9;
    } else if (strcmp(component, "PHE") == 0) {
        return 12;
    } else if (strcmp(component, "PRO") == 0) {
        return 7;
    } else if (strcmp(component, "SER") == 0) {
        return 8;
    } else if (strcmp(component, "THR") == 0) {
        return 9;
    } else if (strcmp(component, "TRP") == 0) {
        return 16;
    } else if (strcmp(component, "TYR") == 0) {
        return 14;
    } else if (strcmp(component, "VAL") == 0) {
        return 8;
    } else if (strcmp(component, "SOL") == 0) {
        return 3;
    } else if (strcmp(component, "GEL") == 0) {
        return 1;
    } else {
        printf("\n!ERROR!:The amino acid name is not found in Library! Check the input coordinate file!\n");
        exit(EXIT_FAILURE);
    }
    
    return -1;
}


int AAModel(char *component) {
    
    if (strcmp(component, "ALA") == 0) {
        return 1;
    } else if (strcmp(component, "ARG") == 0) {
        return 2;
    } else if (strcmp(component, "ASN") == 0) {
        return 3;
    } else if (strcmp(component, "ASP") == 0) {
        return 4;
    } else if (strcmp(component, "CYS") == 0) {
        return 5;
    } else if (strcmp(component, "GLN") == 0) {
        return 6;
    } else if (strcmp(component, "GLU") == 0) {
        return 7;
    } else if (strcmp(component, "GLY") == 0) {
        return 8;
    } else if (strcmp(component, "HIS") == 0) {
        return 9;
    } else if (strcmp(component, "ILE") == 0) {
        return 10;
    } else if (strcmp(component, "LEU") == 0) {
        return 11;
    } else if (strcmp(component, "LYS") == 0) {
        return 12;
    } else if (strcmp(component, "MET") == 0) {
        return 13;
    } else if (strcmp(component, "PHE") == 0) {
        return 14;
    } else if (strcmp(component, "PRO") == 0) {
        return 15;
    } else if (strcmp(component, "SER") == 0) {
        return 16;
    } else if (strcmp(component, "THR") == 0) {
        return 17;
    } else if (strcmp(component, "TRP") == 0) {
        return 18;
    } else if (strcmp(component, "TYR") == 0) {
        return 19;
    } else if (strcmp(component, "VAL") == 0) {
        return 20;
    } else if (strcmp(component, "SOL") == 0) {
        return 21;
    } else if (strcmp(component, "GEL") == 0) {
        return 22;
    } else if (strcmp(component, "PHO") == 0) {
        //hydrophobic surface
        return 23;
    } else if (strcmp(component, "PHI") == 0) {
        //hydrophilic surface
        return 24;
    } else if (strcmp(component, "RPH") == 0) {
        //relatively hydrophilic surface
        return 25;
    } else {
        printf("\n!ERROR!:The amino acid name is not found in Library! Check the input coordinate file!\n");
        exit(EXIT_FAILURE);
    }
    
    return -1;
}


int AtomModel(char *type) {
    
    if (strcmp(type, "OZB") == 0) {
        return 26;
    } else if (strcmp(type, "NZB") == 0) {
        return 18;
    } else if (strcmp(type, "CRB") == 0) {
        return 6;
    } else if (strcmp(type, "CA") == 0) {
        return 9;
    } else if (strcmp(type, "CA1") == 0) {
        return 10;
    } else if (strcmp(type, "CA2") == 0) {
        return 11;
    } else if (strcmp(type, "CA3") == 0) {
        return 12;
    } else if (strcmp(type, "CA1P") == 0) {
        return 14;
    } else if (strcmp(type, "CA2P") == 0) {
        return 15;
    } else if (strcmp(type, "CM1") == 0) {
        return 13;
    } else if (strcmp(type, "CM1P") == 0) {
        return 16;
    } else if (strcmp(type, "CR") == 0) {
        return 5;
    } else if (strcmp(type, "CRN") == 0) {
        return 7;
    } else if (strcmp(type, "CRNP") == 0) {
        return 8;
    } else if (strcmp(type, "NCK") == 0) {
        return 21;
    } else if (strcmp(type, "NCR") == 0) {
        return 22;
    } else if (strcmp(type, "NR") == 0) {
        return 23;
    } else if (strcmp(type, "NZNP") == 0) {
        return 19;
    } else if (strcmp(type, "NZNQ") == 0) {
        return 20;
    } else if (strcmp(type, "NZ") == 0) {
        return 17;
    } else if (strcmp(type, "OC") == 0) {
        return 28;
    } else if (strcmp(type, "OW") == 0) {
        return 29;
    } else if (strcmp(type, "OZH") == 0) {
        return 25;
    } else if (strcmp(type, "OZ") == 0) {
        return 24;
    } else if (strcmp(type, "SGNP") == 0) {
        return 30;
    } else if (strcmp(type, "SG") == 0) {
        return 31;
    } else if (strcmp(type, "OZB_HB") == 0) {
        return 27;
    } else if (strcmp(type, "HB") == 0) {
        return 2;
    } else if (strcmp(type, "HO") == 0) {
        return 1;
    } else if (strcmp(type, "HN") == 0) {
        return 1;
    } else if (strcmp(type, "HS") == 0) { //hydrogen on sulfur
        return 1;
    } else if (strcmp(type, "HW") == 0) {
        return 1;
    } else if (strcmp(type, "H") == 0) {
        return 1;
    } else if (strcmp(type, "H_HB") == 0) {
        return 3;
    } else if (strcmp(type, "HB_HB") == 0) {
        return 4;
    } else if (strcmp(type, "GEL") == 0) {
        return 32;
    } else {
        printf("\nNo atom type is matched! Check the atom outerType assignment! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return -1;
}


int HBModel(struct AtomStr *atom1, struct AtomStr *atom2) {
    struct AtomStr *switcher;
    
    if (atom1->dynamic->HB.role == 'A' ||
		atom1->dynamic->HB.role == 'S') { //atom1 -> donator, atom2 -> acceptor
        switcher = atom1;
        atom1 = atom2;
        atom2 = switcher;
    }
    
    int atom1_OrigType = AtomTypeChange(atom1->property->typeofAtom, 0);
    int atom2_OrigType = AtomTypeChange(atom2->property->typeofAtom, 0);
    
    if (atom1_OrigType == 2 /*HB*/ &&
        atom2_OrigType == 26 /*OZB*/) {
        return 1;
    } else if (atom1_OrigType == 2) { //backbone H, HB
        switch (atom2->property->typeofAA) {
            case 4:  //ASP
            case 7:  //GLU
                return 2;
                break;
            
            case 16: //SER
            case 17: //THR
            case 19: //TYR
                return 3;
                
            case 9:  //HIS
                return 4;
                
            case 5: //CYS
                return 10;

			case 3: //ASN
			case 6: //GLN
				    //no HB between ASN/GLN and backbone H
				return -1;
                
            default:
                break;
        }
    } else if (atom2_OrigType == 26) { //backbone O, OZB
        return 5;
    } else if (atom1_OrigType == 31 &&
               atom2_OrigType == 31) { //SG
        return 9;
    } else if (atom1_OrigType == 1) { //sidechain H
        switch (atom2->property->typeofAA) {
            case 4:  //ASP
            case 7:  //GLU
                return 6;
                break;
                
            case 16: //SER
            case 17: //THR
            case 19: //TYR
                return 7;
                
            case 9:  //HIS
                return 8;

            case 3:  //ASN
            case 6:  //GLN
                return 11;
                
            default:
                break;
        }
    }
    
    printf("!!ERROR!!: cannot find the matched HB pair for atom %2i(%2i%2s) and atom %2i(%2i%2s)!\n",
           atom1->property->num, atom1->property->sequence.aminoacidNum, atom1->property->name,
           atom2->property->num, atom2->property->sequence.aminoacidNum, atom2->property->name);
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    
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

