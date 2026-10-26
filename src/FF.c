//
//  FF.c
//  sDMD
//
//  Created by Size on 2020/7/19.
//  Copyright © 2020 Size Zheng. All rights reserved.
//

#include "DMD.h"
#include <stdio.h>

#define MAX_R       6.5    //cut-off radius
#define BDMAX_P     0.2    //max bond potential, kcal/mol
#define AGMAX_P     0.2    //max angle potential, kcal/mol
#define IMMAX_P     2.0    //max improper potential, kcal/mol
#define LJMAX_P     1.5    //max LJ potential, kcal/mol
#define DHMAX_N     5      //max dihedral step-function steps
#define FF          "CHARMM"

#define KJ2KCAL(val) \
val * 0.23900574;


static int FindNextBondList(int target, struct PropertyStr *property, double dmin, double dmax);
static int FindNextConstrList(int target, struct PropertyStr *property, int step, char *data);
static double FindMaxPotDis(char *type, double *parameter, double setP);
static double FindMinPotDis(char *type, double *parameter, double setP);
static double CalculateBond(char *type, double kb, double b0, double dis);
static double CalculateAngle(char *type, double ktheta, double theta0, double kub, double r0, double *bond, double angle);
static double CalculateDihedral(char *type, double kphi, double n, double delta, double phi);
static double CalculateImproper(char *type, double komega, double omega0, double omega);
static double CalculateDihedralDis(double *bond, double *angle, double deg);
static double CalculateLJ(char *type, double sigma, double radius, double dis);
static double ConvertDeg2Rad(double value);
static double ConvertRad2Deg(double value);
static double CalculateAngleCosineLaw(double a, double b, double c);
static void CreateBondData(char *FFType, char atomType[4][16], struct FFBondStr* bond, FILE *outputFile);
static void CreateAngleData(char *FFType, char atomType[4][16], double *bond, struct FFAngleStr *angle, FILE *outputFile);
static void CreateDihedralData(char *FFType, char atomType[4][16], double *bond, double *angle, struct FFDihedralStr *dihedral, FILE *outputFile);
static void CreateImproperData(char *FFType, char atomType[4][16], double *bond, double *angle, struct FFImproperStr *improper, FILE *outputFile);
static void PrintPotCur(char *FFMethod, char *type, char *atom1, char *atom2,
                        struct FFBondStr     *bond,
                        struct FFAngleStr    *angle,
                        struct FFDihedralStr *dihedral,
                        struct FFImproperStr *improper);
static void PrintDisPot(char *type, char *atom1, char *atom2,
                        struct FFBondStr     *bond,
                        struct FFAngleStr    *angle,
                        struct FFDihedralStr *dihedral,
                        struct FFImproperStr *improper);
static void SortList(struct DisPotStr *list);


void ReadFFInfo(char *method) {
    char buffer[PATH_LEN];
    char name[4][16];
    int num, pos, type[4];
    FILE *inputFile, *outputFile = NULL;

    sprintf(buffer, "%s/Library_%s/FF_Info.txt", datadir, FFLibrary);
    inputFile = fopen(buffer, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find FF_Info.txt file in %s/Library_%s/! %s:%i\n",
               datadir,
               FFLibrary,
               __FILE__, __LINE__);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }
    
#ifdef DEBUG_FF
    sprintf(buffer, "%s/BondedPotentialTable.txt", savedir);
    outputFile = fopen(buffer, "w");
#endif

    while(fgets(buffer, sizeof(buffer), inputFile)) {
        if (strncmp(buffer, "#Res", 4) == 0) {
            fscanf(inputFile, "%i\n", &num);

            FFInfo.group = (struct FFBeadStr*)calloc(num + 1, sizeof(struct FFBeadStr));
            FFInfo.group[0].val[0] = num;
            for (int n = 1; n <= num; n ++) {
                //val[0] = group index num
                //val[1] = bead num in group
                fscanf(inputFile, "%s%i%i\n", FFInfo.group[n].name, &FFInfo.group[n].val[0], &FFInfo.group[n].val[1]);
            }
        } else if (strncmp(buffer, "#Atom", 5) == 0) {
            fscanf(inputFile, "%i\n", &num);

            FFInfo.bead = (struct FFBeadStr*)calloc(num + 1, sizeof(struct FFBeadStr));
            FFInfo.bead[0].val[0] = num;
            for (int n = 1; n <= num; n ++) {
                //val[0] = bead index num
                fscanf(inputFile, "%s%i\n", FFInfo.bead[n].name, &FFInfo.bead[n].val[0]);
            }
        } else if (strncmp(buffer, "#Bond", 5) == 0) {
            _CALLOC2(struct FFBondStr, FFInfo.bond, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
            
            while (fgets(buffer, sizeof(buffer), inputFile) &&
                   buffer[0] == '#');
            sscanf(buffer, "%i", &num);
            
            for (int i = 0; i < num; i ++) {
                fgets(buffer, sizeof(buffer), inputFile);
                sscanf(buffer, "%s%s%n", name[0], name[1], &pos);
                type[0] = FindAtomIndex(name[0]);
                type[1] = FindAtomIndex(name[1]);
                
                sscanf(buffer + pos, "%lf%lf",
                       &FFInfo.bond[type[0]][type[1]].kb,
                       &FFInfo.bond[type[0]][type[1]].b0);
                
                sprintf(FFInfo.bond[type[0]][type[1]].name[0], "%s", name[0]);
                sprintf(FFInfo.bond[type[0]][type[1]].name[1], "%s", name[1]);
                
                CreateBondData(method, name, &FFInfo.bond[type[0]][type[1]], outputFile);
                FFInfo.bond[type[1]][type[0]] = FFInfo.bond[type[0]][type[1]];
            }
        } else if (strncmp(buffer, "#Angle", 6) == 0) {
            _CALLOC2(struct FFAngleStr, FFInfo.angle, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
            
            while (fgets(buffer, sizeof(buffer), inputFile) &&
                   buffer[0] == '#');
            sscanf(buffer, "%i", &num);
            
            for (int i = 0; i < num; i ++) {
                fgets(buffer, sizeof(buffer), inputFile);
                sscanf(buffer, "%s%s%s%n", name[0], name[1], name[2], &pos);
                type[0] = FindAtomIndex(name[0]);
                type[1] = FindAtomIndex(name[1]);
                type[2] = FindAtomIndex(name[2]);
                
                sprintf(FFInfo.angle[type[0]][type[2]].name[0],        "%s", name[0]);
                sprintf(FFInfo.angle[type[0]][type[2]].interAtom.name, "%s", name[1]);
                sprintf(FFInfo.angle[type[0]][type[2]].name[1],        "%s", name[2]);
                FFInfo.angle[type[0]][type[2]].interAtom.val[0] = type[1];
                
                sscanf(buffer + pos, "%lf%lf%lf%lf",
                       &FFInfo.angle[type[0]][type[2]].ktheta,
                       &FFInfo.angle[type[0]][type[2]].theta0,
                       &FFInfo.angle[type[0]][type[2]].r0,
                       &FFInfo.angle[type[0]][type[2]].kub);
                
                double bond[2]; //bond1, bond2
                bond[0] = FFInfo.bond[type[0]][type[1]].b0;
                bond[1] = FFInfo.bond[type[1]][type[2]].b0;
                
                CreateAngleData(method, name, bond, &FFInfo.angle[type[0]][type[2]], outputFile);
                FFInfo.angle[type[2]][type[0]] = FFInfo.angle[type[0]][type[2]];
            }
        } else if (strncmp(buffer, "#Dihedral", 9) == 0) {
            _CALLOC2(struct FFDihedralStr, FFInfo.dihedral, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
            
            while (fgets(buffer, sizeof(buffer), inputFile) &&
                   buffer[0] == '#');
            sscanf(buffer, "%i", &num);
            
            for (int i = 0; i < num; i ++) {
                fgets(buffer, sizeof(buffer), inputFile);
                sscanf(buffer, "%s%s%s%s%n", name[0], name[1], name[2], name[3], &pos);
                type[0] = FindAtomIndex(name[0]);
                type[1] = FindAtomIndex(name[1]);
                type[2] = FindAtomIndex(name[2]);
                type[3] = FindAtomIndex(name[3]);

                sprintf(FFInfo.dihedral[type[0]][type[3]].name[0], "%s", name[0]);
                sprintf(FFInfo.dihedral[type[0]][type[3]].name[1], "%s", name[3]);
                sprintf(FFInfo.dihedral[type[0]][type[3]].interAtom[0].name, "%s", name[1]);
                sprintf(FFInfo.dihedral[type[0]][type[3]].interAtom[1].name, "%s", name[2]);
                FFInfo.dihedral[type[0]][type[3]].interAtom[0].val[0] = type[1];
                FFInfo.dihedral[type[0]][type[3]].interAtom[1].val[0] = type[2];

                sscanf(buffer + pos, "%lf%lf%lf",
                       &FFInfo.dihedral[type[0]][type[3]].kphi,
                       &FFInfo.dihedral[type[0]][type[3]].n,
                       &FFInfo.dihedral[type[0]][type[3]].delta);
                
                //bond1, bond2, bond3
                FFInfo.dihedral[type[0]][type[3]].data.bond[0] = FFInfo.bond[type[0]][type[1]].b0;
                FFInfo.dihedral[type[0]][type[3]].data.bond[1] = FFInfo.bond[type[1]][type[2]].b0;
                FFInfo.dihedral[type[0]][type[3]].data.bond[2] = FFInfo.bond[type[2]][type[3]].b0;

                //angle1, angle2
                FFInfo.dihedral[type[0]][type[3]].data.angle[0] = FFInfo.angle[type[0]][type[2]].theta0;
                FFInfo.dihedral[type[0]][type[3]].data.angle[1] = FFInfo.angle[type[1]][type[3]].theta0;

                CreateDihedralData(method, name,
                                   FFInfo.dihedral[type[0]][type[3]].data.bond,
                                   FFInfo.dihedral[type[0]][type[3]].data.angle,
                                  &FFInfo.dihedral[type[0]][type[3]],
                                   outputFile);
                FFInfo.dihedral[type[3]][type[0]] = FFInfo.dihedral[type[0]][type[3]];
            }
        } else if (strncmp(buffer, "#Improper", 9) == 0) {
            _CALLOC2(struct FFImproperStr, FFInfo.improper, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
            
            while (fgets(buffer, sizeof(buffer), inputFile) &&
                   buffer[0] == '#');
            sscanf(buffer, "%i", &num);
            
            for (int i = 0; i < num; i ++) {
                fgets(buffer, sizeof(buffer), inputFile);
                sscanf(buffer, "%s%s%s%s%n", name[0], name[1], name[2], name[3], &pos);
                type[0] = FindAtomIndex(name[0]);
                type[1] = FindAtomIndex(name[1]);
                type[2] = FindAtomIndex(name[2]);
                type[3] = FindAtomIndex(name[3]);
                
                sprintf(FFInfo.improper[type[0]][type[3]].name[0], "%s", name[0]);
                sprintf(FFInfo.improper[type[0]][type[3]].name[1], "%s", name[3]);
                sprintf(FFInfo.improper[type[0]][type[3]].interAtom[0].name, "%s", name[1]);
                sprintf(FFInfo.improper[type[0]][type[3]].interAtom[1].name, "%s", name[2]);
                FFInfo.improper[type[0]][type[3]].interAtom[0].val[0] = type[1];
                FFInfo.improper[type[0]][type[3]].interAtom[1].val[0] = type[2];
                
                sscanf(buffer + pos, "%lf%lf",
                       &FFInfo.improper[type[0]][type[3]].komega,
                       &FFInfo.improper[type[0]][type[3]].omega0);
                
                //bond1, bond2, bond3
                FFInfo.improper[type[0]][type[3]].data.bond[0] = FFInfo.bond[type[0]][type[1]].b0;
                FFInfo.improper[type[0]][type[3]].data.bond[1] =
                    (FFInfo.angle[type[1]][type[2]].dmin + FFInfo.angle[type[1]][type[2]].dmax) / 2;
                FFInfo.improper[type[0]][type[3]].data.bond[2] =
                    (FFInfo.angle[type[2]][type[3]].dmin + FFInfo.angle[type[2]][type[3]].dmax) / 2;

                double bond[2];
                bond[0] = FFInfo.bond[type[0]][type[2]].b0;
                bond[1] = (FFInfo.angle[type[1]][type[3]].dmin + FFInfo.angle[type[1]][type[3]].dmax) / 2;
                
                //angle1, angle2
                FFInfo.improper[type[0]][type[3]].data.angle[0] =
                    CalculateAngleCosineLaw(FFInfo.improper[type[0]][type[3]].data.bond[0],
                                            FFInfo.improper[type[0]][type[3]].data.bond[1],
                                            bond[0]);
                FFInfo.improper[type[0]][type[3]].data.angle[1] =
                    CalculateAngleCosineLaw(FFInfo.improper[type[0]][type[3]].data.bond[1],
                                            FFInfo.improper[type[0]][type[3]].data.bond[2],
                                            bond[1]);
                
                CreateImproperData(method, name,
                                   FFInfo.improper[type[0]][type[3]].data.bond,
                                   FFInfo.improper[type[0]][type[3]].data.angle,
                                  &FFInfo.improper[type[0]][type[3]],
                                   outputFile);
                FFInfo.improper[type[2]][type[0]] = FFInfo.improper[type[0]][type[2]];
            }
        } else if (strncmp(buffer, "#LJ", 3) == 0) {
            FFInfo.LJ = (struct FFLJStr *)calloc((FFInfo.bead[0].val[0] + 1), sizeof(struct FFLJStr));
            
            while (fgets(buffer, sizeof(buffer), inputFile) &&
                   buffer[0] == '#');
            sscanf(buffer, "%i", &num);
            
            for (int i = 0; i < num; i ++) {
                fgets(buffer, sizeof(buffer), inputFile);
                sscanf(buffer, "%s%n", name[0], &pos);
                type[0] = FindAtomIndex(name[0]);
                
                sscanf(buffer + pos, "%lf%lf%lf%lf",
                       &FFInfo.LJ[type[0]].mass,
                       &FFInfo.LJ[type[0]].charge,
                       &FFInfo.LJ[type[0]].sigma,
                       &FFInfo.LJ[type[0]].epsilon);
                
                sprintf(FFInfo.LJ[type[0]].name, "%s", name[0]);
            }
        }
    }

    fclose(inputFile);
#ifdef DEBUG_FF
    fclose(outputFile);
#endif
    
    return;
}


int ReadAtomInfo(char *info, int resNum) {
    int pos, num;
    char targetAtomName[16];
    
    //*.property->type may have changed during the previous simulations
    if (strcmp(neworcontinue, "new") == 0) {
        sscanf(info, "%s%n", targetAtomName, &pos);
        num = FindTargetAtom(resNum, targetAtomName);
        
        if (num == FALSE && strncmp(targetAtomName, "H", 1) == 0) {
            printf("!WARNING![%li]: atom %s at molecule #%i residue #%i. if it is a terminal N, then no H attached! if it is not, then there could be a H missing!\n", ++warningsum, targetAtomName, residue[resNum].molNum, resNum);
            return -1;
        }
        
        sscanf(info + pos, "%s%lf%s%s%n", targetAtomName,
                                             &atom[num].property->mass,
                                             atom[num].property->extraProperty[1],
                                             atom[num].property->extraProperty[2],
                                             &pos);
        atom[num].property->typeofAtom[1] = FindAtomIndex(targetAtomName);
        sprintf(atom[num].property->typeName[1], "%s", targetAtomName);
        
    } else {
        return 1;
    }
    
    return 0;
}


int ReadBondInfo(char *info, int resNum) {
    int pos, num, flag;
    char targetAtomName[2][16];
    
    sscanf(info, "%s%s%n",
           targetAtomName[0],
           targetAtomName[1],
           &pos);
    
    if ((num = FindTargetAtom(resNum, targetAtomName[0])) > 0) { //available and legal
        tmpInt = FindTargetAtom(resNum, targetAtomName[1]);
        
        if (tmpInt == FALSE) {
            //if it is a light hydrogen
            if (strncmp(targetAtomName[1], "H", 1) == 0) {
                return -1;
            }
            printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[1], __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        } else if (tmpInt == -1) { //not available (terminal atom) but legal
            return 1;
        }
        
        int type[2];
        type[0] = FindAtomIndex(atom[num   ].property->typeName[1]);
        type[1] = FindAtomIndex(atom[tmpInt].property->typeName[1]);
        
        flag = FindNextBondList(tmpInt, atom[num].property,
                                FFInfo.bond[type[0]][type[1]].dmin,
                                FFInfo.bond[type[0]][type[1]].dmax);

        if (flag == 1) {
            if (connectionMap[num][tmpInt] & CONSTRAINT_CONNECT) {
                connectionMap[num][tmpInt] ^= CONSTRAINT_CONNECT;
            }
            connectionMap[num][tmpInt] |= BOND_CONNECT;
        }
        
    } else if (num == FALSE) {
        if (strncmp(targetAtomName[0], "H", 1) == 0) {
            return -1;
        }
        printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", targetAtomName[0], residue[resNum].startAtomNum, residue[resNum].endAtomNum, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return 0;
}


int ReadAngleInfo(char *info, int resNum) {
    int pos, num, flag;
    char targetAtomName[3][16];
    
    sscanf(info, "%s%s%s%n",
           targetAtomName[0],
           targetAtomName[1],
           targetAtomName[2],
           &pos);
    
    if ((num = FindTargetAtom(resNum, targetAtomName[0])) > 0) { //available and legal
        tmpInt = FindTargetAtom(resNum, targetAtomName[2]);
        
        if (tmpInt == FALSE) {
            //if it is a light hydrogen
            if (strncmp(targetAtomName[2], "H", 1) == 0) {
                return -1;
            }
            printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[2], __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        } else if (tmpInt == -1) { //not available (terminal atom) but legal
            return 1;
        }
        
        int interAtom = FindTargetAtom(resNum, targetAtomName[1]);
        if (interAtom == FALSE) {
            printf("!!ERROR!!: molecule #%i residue #%i cannot find the inter-atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[1], __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        int type[2];
        type[0] = FindAtomIndex(atom[num      ].property->typeName[1]);
        type[1] = FindAtomIndex(atom[tmpInt   ].property->typeName[1]);
        
        flag = FindNextBondList(tmpInt, atom[num].property,
                                FFInfo.angle[type[0]][type[1]].dmin,
                                FFInfo.angle[type[0]][type[1]].dmax);

        if (flag == 1) { //if it is a CONSTRAINT previously, cancel it
            if (connectionMap[num][tmpInt] & CONSTRAINT_CONNECT) {
                connectionMap[num][tmpInt] ^= CONSTRAINT_CONNECT;
            }
            connectionMap[num][tmpInt] |= BOND_CONNECT;
        }
        
    } else if (num == FALSE) {
        if (strncmp(targetAtomName[0], "H", 1) == 0) {
            return -1;
        }
        printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", targetAtomName[0], residue[resNum].startAtomNum, residue[resNum].endAtomNum, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return 0;
}


int ReadDihedralInfo(char *infoType, char *info, int resNum) {
    int pos, num, flag;
    char targetAtomName[4][16];
    struct DihedralDataStr *thisData = NULL;
    
    sscanf(info, "%s%s%s%s%n",
           targetAtomName[0],
           targetAtomName[1],
           targetAtomName[2],
           targetAtomName[3],
           &pos);
    
    if ((num = FindTargetAtom(resNum, targetAtomName[0])) > 0) { //available and legal
        tmpInt = FindTargetAtom(resNum, targetAtomName[3]);
        
        if (tmpInt == FALSE) {
            //if it is a light hydrogen
            if (strncmp(targetAtomName[3], "H", 1) == 0) {
                return -1;
            }
            printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[3], __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        } else if (tmpInt == -1) { //not available (terminal atom) but legal
            return 1;
        }
        
        //check if the inter-atoms are valid
        int interAtom[2];
        for (int i = 1; i < 3; i ++) {
            interAtom[i - 1] = FindTargetAtom(resNum, targetAtomName[i]);
            if (interAtom[i - 1] == FALSE) {
                printf("!!ERROR!!: molecule #%i residue #%i cannot find the inter-atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[i], __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }

        int type[2];
        type[0] = FindAtomIndex(atom[num   ].property->typeName[1]);
        type[1] = FindAtomIndex(atom[tmpInt].property->typeName[1]);
        
        if (strncmp(infoType, "CONSTR", 6) == 0) {
            thisData = &FFInfo.dihedral[type[0]][type[1]].data;
        } else if (strncmp(infoType, "IMPROP", 6) == 0) {
            thisData = &FFInfo.improper[type[0]][type[1]].data;
        }
        flag = FindNextConstrList(tmpInt, atom[num].property,
                                  thisData->step,
                                  thisData->dataString);
        
        if (flag == 1) {
            if (connectionMap[num][tmpInt] & BOND_CONNECT) {
                connectionMap[num][tmpInt] ^= BOND_CONNECT;
            }
            connectionMap[num][tmpInt] |= CONSTRAINT_CONNECT;
        }
    } else if (num == FALSE) {
        if (strncmp(targetAtomName[0], "H", 1) == 0) {
            return -1;
        }
        printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", targetAtomName[0], residue[resNum].startAtomNum, residue[resNum].endAtomNum, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return 0;
}


int ReadHBInfo(char *infoType, char *info, int resNum) {
    int curPos = 0, pos, num;
    char targetAtomName[2][16];

    if (strncmp(infoType, "HBCNT", 5) == 0) {
        int i = infoType[5] - '0';
        char changeType[16];

        sscanf(info + curPos, "%s%s%n",
               targetAtomName[0],
               changeType,
               &pos);
        curPos += pos;

        if ((num = FindTargetAtom(resNum, targetAtomName[0])) > 0) { //available and legal
            atom[num].dynamic->HB.neighbor[0] = i;
            atom[num].property->typeofAtom[2] = FindAtomIndex(changeType);
            sprintf(atom[num].property->typeName[2], "%s", changeType);
            sscanf(info + curPos, "%s%n", &atom[num].dynamic->HB.role, &pos); curPos += pos;

            int count = 0;
            for (int n = 0; n < i; n ++) {
                sscanf(info + curPos, "%s%n", targetAtomName[1], &pos); curPos += pos;
                tmpInt = FindTargetAtom(resNum, targetAtomName[1]);
                if (tmpInt == FALSE) {
                    printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[1], __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                } else if (tmpInt == -1) {
                    if (--atom[num].dynamic->HB.neighbor[0] == 0)
                        atom[num].dynamic->HB.role = '\0';
                    continue;
                }
                atom[num].dynamic->HB.neighbor[++count] = tmpInt;
            }
        } else if (strncmp(targetAtomName[0], "H", 1) == 0) {
            return -1;
        } else {
            printf("!!ERROR!!: cannot find atom %s near atom #%i! check the model library or the coordinate file! %s:%i\n", targetAtomName[0], num, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    } else if (strncmp(infoType, "HBRMV", 5) == 0) {
        int i = infoType[5] - '0';

        sscanf(info + curPos, "%s%n",
               targetAtomName[0],
               &pos);
        curPos += pos;
        
        if ((num = FindTargetAtom(resNum, targetAtomName[0])) > 0) {
            atom[num].dynamic->HBNeighbor.removedAtom[0] = i;
            for (int n = 0; n < i; n ++) {
                sscanf(info + curPos, "%s%n", targetAtomName[1], &pos); curPos += pos;
                tmpInt = FindTargetAtom(resNum, targetAtomName[1]);
                if (tmpInt == FALSE) {
                    printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, targetAtomName[1], __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                } else if (tmpInt == -1) {
                    atom[num].dynamic->HBNeighbor.removedAtom[0] --;
                    continue;
                }
                FindElemInList(atom[num].dynamic->HBNeighbor.removedAtom, tmpInt, 1, i, 1);
                FindStrInList(atom[tmpInt].property->extraProperty, "HBRMV", 1);
            }
        }
    } else if (strcmp(infoType, "HBMAX") == 0) {
        // HBMAX atom n: the atom can hold n HBs at the same time (default 1)
        int n = 0;
        sscanf(info + curPos, "%s%i%n", targetAtomName[0], &n, &pos);
        curPos += pos;

        if ((num = FindTargetAtom(resNum, targetAtomName[0])) > 0) {
            if (n < 1 || n > HB_MAX_PARTNER) {
                printf("!!ERROR!!: HBMAX of atom %s must be between 1 and %i! %s:%i\n", targetAtomName[0], HB_MAX_PARTNER, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            atom[num].dynamic->HB.maxPartner = n;
        }
    } else {
        printf("!!ERROR!!: HB info type is invalid!\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return 0;
}


void CreateBondData(char *FFType, char atomType[4][16], struct FFBondStr *bond, FILE *outputFile) {
    double maxDis = 0;
    
    if (strncmp(FFType, "CHARMM", 6) == 0 ||
        strcmp(FFType, "Crosslink") == 0) {
        maxDis = bond->b0;
        while ((tmpDouble = CalculateBond(FFType,
                                          bond->kb,
                                          bond->b0,
                                          maxDis))
                < BDMAX_P) {
            maxDis += 0.002;
        }
        bond->delta = maxDis / bond->b0 - 1;
        bond->dmin  = bond->b0 * (1 - bond->delta);
        bond->dmax  = bond->b0 * (1 + bond->delta);
        
        if (outputFile) {
            fprintf(outputFile, "   BOND %-7s %-7s", atomType[0], atomType[1]);
            fprintf(outputFile, "%-10.6lf %-10.6lf\n", bond->b0, bond->delta);
        }
    }
    
#ifdef DEBUG_FF
    PrintPotCur(FFType, "Bond", atomType[0], atomType[1], bond, NULL, NULL, NULL);
    PrintDisPot("Bond", atomType[0], atomType[1], bond, NULL, NULL, NULL);
#endif
    
    return;
}


void CreateAngleData(char *FFType, char atomType[4][16], double *bond, struct FFAngleStr *angle, FILE *outputFile) {
    
    if (strncmp(FFType, "CHARMM", 6) == 0) {
        angle->maxDeg = angle->theta0;
        while ((tmpDouble = CalculateAngle(FFType,
                                           angle->ktheta,
                                           angle->theta0,
                                           angle->kub,
                                           angle->r0,
                                           bond,
                                           angle->maxDeg))
               < AGMAX_P) {
            angle->maxDeg += 0.1;
        }
        //since the potential profile is like a parabola, the min and max degree is symetrical to theta0
        angle->minDeg = 2 * angle->theta0 - angle->maxDeg;
        
        angle->dmin = sqrt(bond[0] * bond[0] + bond[1] * bond[1] -
                           2 * bond[0] * bond[1] * cos(ConvertDeg2Rad(angle->minDeg)));
        angle->dmax = sqrt(bond[0] * bond[0] + bond[1] * bond[1] -
                           2 * bond[0] * bond[1] * cos(ConvertDeg2Rad(angle->maxDeg)));
        angle->delta = 2 * angle->dmax / (angle->dmin + angle->dmax) - 1;
                
        if (outputFile) {
            fprintf(outputFile, "  ANGLE %-7s %-7s", atomType[0], atomType[2]);
            fprintf(outputFile, "%-10.6lf %-10.6lf\n", 0.5 * (angle->dmax + angle->dmin), angle->delta);
        }
    }
    
#ifdef DEBUG_FF
    PrintPotCur(FFType, "Angle", atomType[0], atomType[2], NULL, angle, NULL, NULL);
    PrintDisPot("Angle", atomType[0], atomType[2], NULL, angle, NULL, NULL);
#endif
    
    return;
}


void CreateDihedralData(char *FFType, char atomType[4][16], double *bond, double *angle, struct FFDihedralStr *dihedral, FILE *outputFile) {
    // bond[3]: bond1,  bond2,  bond3
    //angle[2]: angle1, angle2
    
    int stepRecord = DHMAX_N;
    int pos = 0, curPos = 0;
    double midDeg = 0, degGap;
    double minD = INFINIT, maxD = 0;
    struct DisPotStr list[DHMAX_N + 2];
    
    if (strncmp(FFType, "CHARMM", 6) == 0 ||
        strcmp(FFType, "Crosslink") == 0) {
        //Find the min and max distances
        //min -> deg = 0
        minD = CalculateDihedralDis(bond, angle, 0);
        
        //max -> deg = 180
        maxD = CalculateDihedralDis(bond, angle, 180);
        
        list[0].dis = minD;
        list[0].val = 0.0;
        list[0].pot = 0;
        
        //set the max bound for the dihedral angle
        list[DHMAX_N + 1].dis = 4.0; //cut-off
        list[DHMAX_N + 1].val = 180.0;
        list[DHMAX_N + 1].pot = 0;
        
        degGap = 180.0 / DHMAX_N;
        midDeg = 0.5 * degGap;
        for (int i = 1, n = 1; i <= stepRecord; i ++, n ++) {
            list[i].val = n * degGap;
            list[i].dis = CalculateDihedralDis(bond, angle, list[i].val);
            list[i].pot = CalculateDihedral("CHARMM19", dihedral->kphi, dihedral->n, dihedral->delta, midDeg);
            midDeg += degGap;
            
            //just in case that the neighboring potential jumps had the same value
            if (fabs(list[i].pot - list[i - 1].pot) < ZERO) {
                i --;
                stepRecord --;
                continue;
            }
        }
        
        //export the FF data as the old format
        sprintf(dihedral->data.dataString, "%-10.6lf %n", list[0].dis, &pos); //write the min distance
        curPos += pos;

        for (int i = 1; i <= stepRecord; i ++) {
            sprintf(dihedral->data.dataString + curPos, "%-10.6lf %-10.6lf %n",
                    list[i].dis,
                    (list[i].pot - list[i + 1].pot),
                    &pos);
            curPos += pos;
        }
        sprintf(dihedral->data.dataString + curPos, "%-10.6lf", list[DHMAX_N + 1].dis);
        dihedral->data.step = stepRecord;
        
        dihedral->list = (struct DisPotStr *)calloc(stepRecord + 2, sizeof(struct DisPotStr));
        for (int i = 0; i <= stepRecord; i ++) {
            dihedral->list[i] = list[i];
        }
        dihedral->list[stepRecord + 1] = list[DHMAX_N + 1];
        
        if (outputFile) {
            fprintf(outputFile, "CONSTR%i %-7s %-7s", stepRecord, atomType[0], atomType[3]);
            fprintf(outputFile, "%s\n", dihedral->data.dataString);
        }
    }
    
#ifdef DEBUG_FF
    PrintPotCur(FFType, "Dihedral", atomType[0], atomType[3], NULL, NULL, dihedral, NULL);
    PrintDisPot("Dihedral", atomType[0], atomType[3], NULL, NULL, dihedral, NULL);
#endif
    
    return;
}


void CreateImproperData(char *FFType, char atomType[4][16], double *bond, double *angle, struct FFImproperStr *improper, FILE *outputFile) {
    // bond[3]: bond1,  bond2,  bond3
    //angle[2]: angle1, angle2
    
    double minDeg, maxDeg, curDeg;
    double minD = INFINIT, maxD = 0.0;
    struct DisPotStr list[2];
    
    if (strncmp(FFType, "CHARMM", 6) == 0 ||
        strcmp(FFType, "Crosslink") == 0) {
        maxDeg = improper->omega0;
        while ((tmpDouble = CalculateImproper(FFType,
                                              improper->komega,
                                              improper->omega0,
                                              maxDeg))
               < IMMAX_P) {
            maxDeg += 0.1;
        }
        minDeg = 2 * improper->omega0 - maxDeg;
        
        curDeg = minDeg;
        while (curDeg <= maxDeg) {
            tmpDouble = CalculateDihedralDis(bond, angle, curDeg);
            if (minD > tmpDouble) minD = tmpDouble;
            if (maxD < tmpDouble) maxD = tmpDouble;
            curDeg += 0.1;
        }
        
        //pay attention here: the minD/maxD is not exactly corresponding to the minDeg/maxDeg
        //since the improper dihedral degree does not conver the whole pi (rad) range, unlike
        //the regular dihedrals. Thus, the minDeg/maxDeg is the min/max degree the angle is
        //allowed, but it is not necessary to be corresponding to the minD/maxD.
        //Theoretically, the minD would be at 0 rad, while the maxD would be at pi rad.
        list[0].dis = minD;
        list[0].val = minDeg;
        list[0].pot = 0;
        
        list[1].dis = maxD;
        list[1].val = maxDeg;
        list[1].pot = 0;
        
        //export the FF data as the old format
        sprintf(improper->data.dataString, "%-10.6lf %-10.6lf", list[0].dis, list[1].dis);
        improper->data.step = 0;
        
        improper->list = (struct DisPotStr *)calloc(2, sizeof(struct DisPotStr));
        improper->list[0] = list[0];
        improper->list[1] = list[1];
        
        if (outputFile) {
            fprintf(outputFile, " IMPROP %-7s %-7s", atomType[0], atomType[3]);
            fprintf(outputFile, "%s\n", improper->data.dataString);
        }
    }
    
#ifdef DEBUG_FF
    PrintPotCur(FFType, "Improper", atomType[0], atomType[3], NULL, NULL, NULL, improper);
    PrintDisPot("Improper", atomType[0], atomType[3], NULL, NULL, NULL, improper);
#endif
    
    return;
}


void CreateLJData(FILE *outputFile) {
    int flag = 0, step;
    double gap[3] = {0.5, 2.0, 0.01}; //0: potential gap; 1: distance gap; 2: step length
    double parameter[2], minDis, minPDis;
    double dis, disPot[128][2] = {0};

#ifdef DEBUG_FF
    //PrintPotCur(FF, "LJ");
#endif
    
    for (int i = 1; i <= FFInfo.bead[0].val[0]; i ++) {
        for (int j = i; j <= FFInfo.bead[0].val[0]; j ++) {
            parameter[0] = sqrt(FFInfo.LJ[i].epsilon * FFInfo.LJ[j].epsilon);       //epsilon
            parameter[1] = 0.5 * (FFInfo.LJ[i].sigma + FFInfo.LJ[j].sigma);         //sigma
            
            minDis  = FindMaxPotDis(FF, parameter, LJMAX_P);
            minPDis = FindMinPotDis(FF, parameter, 0);

            step = 0; flag = 1;
            dis = minDis + gap[2];
            disPot[step][0] = dis;
            disPot[step++][1] = CalculateLJ(FF, parameter[0], parameter[1], dis);

            while (dis <= MAX_R) {
                dis += gap[2];

                if (step == 128) {
                    printf("!!ERROR!!: too many steps! please increase the step gap. %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }

                disPot[step][0] = dis;
                disPot[step][1] = CalculateLJ(FF, parameter[0], parameter[1], dis);

                if (flag && fabs(disPot[step][0] - minPDis) < gap[2]) { //only do this once
                    disPot[step][0] = minPDis;
                    disPot[step][1] = CalculateLJ(FF, parameter[0], parameter[1], dis);

                    if (disPot[step][0] - disPot[step - 1][0] > gap[1] * 0.1 ||
                        fabs(disPot[step][1] - disPot[step - 1][1]) > gap[0] * 0.1) {
                        step ++;
                    }
                    flag = 0;
                } else if (disPot[step][0] - disPot[step - 1][0] > gap[1] ||
                    fabs(disPot[step][1] - disPot[step - 1][1]) > gap[0]) {
                    step ++;
                }
            }

#ifdef DEBUG_FF
#endif

            fprintf(outputFile, "%8s %8s ", FFInfo.LJ[i].name, FFInfo.LJ[j].name);
            fprintf(outputFile, "%10.6lf ", minDis);

            int stepRecord = step; step = 0;
            while (++step < stepRecord) {
                if (fabs(disPot[step - 1][1] - disPot[step][1]) > 10E-5) {
                    fprintf(outputFile, "%10.6lf %10.6lf ", (disPot[step - 1][0] + disPot[step][0]) * 0.5, (disPot[step - 1][1] - disPot[step][1]));
                }
            }
            if (fabs(disPot[step - 1][1]) > 10E-5)
                fprintf(outputFile, "%10.6lf %10.6lf", disPot[step - 1][0], disPot[step - 1][1]);
            fprintf(outputFile, "\n");
        }
    }
    
    return;
}


double CalculateBond(char *type, double kb, double b0, double dis) {
    if (strcmp(type, "CHARMM36") == 0 ||
        strcmp(type, "Crosslink") == 0) {
        return KJ2KCAL(0.5 * kb * (dis - b0) * (dis - b0) / 100);
    } else if (strcmp(type, "CHARMM19") == 0) {
        return (kb * (dis - b0) * (dis - b0));
    }
    
    printf("!!ERROR!!: type for bond calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


//Urey-Bradley potential
double CalculateAngle(char *type, double ktheta, double theta0, double kub, double r0, double *bond, double angle) {
    double dis = 0;
    double c_theta0, c_angle;
    
    c_theta0 = ConvertDeg2Rad(theta0);
    c_angle  = ConvertDeg2Rad(angle);
    
    if (bond)
        dis = sqrt(bond[0] * bond[0] + bond[1] * bond[1] -
                   2 * bond[0] * bond[1] * cos(c_angle));
    
    if (strcmp(type, "CHARMM36") == 0) {
        return KJ2KCAL(0.5 * ktheta * (c_angle - c_theta0) * (c_angle - c_theta0) +
                       0.5 * kub    * (dis     - r0      ) * (dis     - r0      ) / 100);
    } else if (strcmp(type, "CHARMM19") == 0) {
        return ktheta * (c_theta0 - c_angle) * (c_theta0 - c_angle);
    }
    
    printf("!!ERROR!!: type for angle calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


double CalculateDihedral(char *type, double kphi, double n, double delta, double phi) {
    double c_delta, c_phi;
    
    c_delta = ConvertDeg2Rad(delta);
    c_phi   = ConvertDeg2Rad(phi);
    
    if (strcmp(type, "CHARMM36") == 0) {
        return KJ2KCAL(kphi * (1 + cos(n * c_phi - c_delta)));
    } else if (strcmp(type, "CHARMM19") == 0) {
        return kphi * (1 + cos(n * c_phi - c_delta));
    }
    
    printf("!!ERROR!!: type for proper dihedral angle calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


double CalculateImproper(char *type, double komega, double omega0, double omega) {
    double c_omega0, c_omega;
    
    c_omega0 = ConvertDeg2Rad(omega0);
    c_omega  = ConvertDeg2Rad(omega);
    
    if (strcmp(type, "CHARMM19") == 0) {
        return komega * (c_omega - c_omega0) * (c_omega - c_omega0);
    }
    
    printf("!!ERROR!!: type for improper dihedral calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


double CalculateDihedralDis(double *bond, double *angle, double deg) {
    double dis;
    double point[3][4];
    
    point[1][1] = cos(ConvertDeg2Rad(180 - angle[0])) * bond[0] * -1;
    point[1][2] = 0.0;
    point[1][3] = sin(ConvertDeg2Rad(180 - angle[0])) * bond[0];

    point[2][1] = bond[1] + cos(ConvertDeg2Rad(180 - angle[1])) * bond[2];
    point[2][2] = sin(ConvertDeg2Rad(180 - angle[1])) * bond[2] * sin(ConvertDeg2Rad(180 - deg));
    point[2][3] = sin(ConvertDeg2Rad(180 - angle[1])) * bond[2] * cos(ConvertDeg2Rad(180 - deg)) * -1;

    DOT_MINUS(point[1], point[2], point[0]);
    dis = DOT_PROD(point[0], point[0]);

    return sqrt(dis);
}


double CalculateImproperDis(double *bond, double *angle, double deg) {
    double dis = 0;
    //double pseudoBond[2];
    
    return sqrt(dis);
}


double CalculateLJ(char *type, double epsilon, double sigma, double dis) {
    if (strcmp(type, "CHARMM") == 0 ||
        strcmp(type, "Crosslink") == 0) {
        double ratio = sigma / dis;
        return KJ2KCAL(4 * epsilon * (CalExponential(ratio, 12) - CalExponential(ratio, 6)));
    }
    
    printf("!!ERROR!!: FF type for LJ calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


int FindNextBondList(int target, struct PropertyStr *property, double dmin, double dmax) {
    struct ConstraintStr *thisBond;
    
    thisBond = property->bond;
    while (thisBond != NULL) {
        if (thisBond->connection == target) {
            return 0;
        }
        thisBond = thisBond->next;
    }
    
    struct ConstraintStr *newBond = calloc(1, sizeof(struct ConstraintStr));
    newBond->connection = target;
    newBond->dmin = dmin * dmin;
    newBond->dmax = dmax * dmax;
    newBond->next = (property->bond == NULL) ? NULL : property->bond;
    property->bond = newBond;
    
    return 1;
}


int FindNextConstrList(int target, struct PropertyStr *property, int step, char *data) {
    int curPos = 0, pos;
    struct ConstraintStr *thisConstr = property->constr;
    struct StepPotenStr *thisStep = NULL;
    
    while (thisConstr != NULL) {
        if (thisConstr->connection == target) {
            return 0;
        }
        thisConstr = thisConstr->next;
    }
    
    struct ConstraintStr *newConstr = calloc(1, sizeof(struct ConstraintStr));
    newConstr->connection = target;
    sscanf(data + curPos, "%lf%n", &newConstr->dmin, &pos); curPos += pos;
    newConstr->dmin *= newConstr->dmin;
    
    for (int i = 0; i < step; i ++) {
        struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
        sscanf(data + curPos, "%lf%lf%n", &newStep->d, &newStep->e, &pos); curPos += pos;
        newStep->d *= newStep->d;
        
        if (i == 0) { //append, make sure the length is from small to large
            newConstr->step = newStep;
            thisStep = newStep;
        } else {
            thisStep->next = newStep;
            thisStep = thisStep->next;
        }
    }
    
    sscanf(data + curPos, "%lf%lf%lf%n", &newConstr->dmax, &newConstr->e, &newConstr->accumulated, &pos); curPos += pos;
    newConstr->dmax *= newConstr->dmax;
    newConstr->next = (property->constr == NULL) ? NULL : property->constr;
    property->constr = newConstr;
    
    CalAccumulatedPoten(property->constr->step);
    
    return 1;
}


double FindMaxPotDis(char *type, double *parameter, double setP) {
    double dis = parameter[1] * 0.5;
    
    if (parameter[0] == 0) //epsilon = 0
        return 0;
    
    while (CalculateLJ(type, parameter[0], parameter[1], dis) > setP) {
        dis += 0.02;
    }

    return dis;
}


double FindMinPotDis(char *type, double *parameter, double setP) {
    double min = INFINIT, minDis = 0;
    double dis = 1.0, value;
    
    if (parameter[0] == 0) //epsilon = 0
        return 0;
    
    while (dis < MAX_R) {
        value = CalculateLJ(type, parameter[0], parameter[1], dis);
        
        if (value < min) {
            min    = value;
            minDis = dis;
        }
        dis += 0.001;
    }
    return minDis;
}


void PrintPotCur(char *FFMethod, char *type, char *atom1, char *atom2,
                 struct FFBondStr     *bond,
                 struct FFAngleStr    *angle,
                 struct FFDihedralStr *dihedral,
                 struct FFImproperStr *improper) {
    FILE *potCurFile;
    char name[PATH_LEN];
    double minV, maxV;
    double val, gap, pot, dis;
    
    sprintf(name, "%s/Debug%s-%s-%s_PotCur.txt",
            savedir,
            type,
            atom1,
            atom2);
    potCurFile = fopen(name, "w");
    
    if (strcmp(type, "Bond") == 0) {
        minV = bond->dmin;
        maxV = bond->dmax;
        val = minV;
        gap = 0.0001;
        
        fprintf(potCurFile, "# %-8s %-8s\n", "dis", "energy");
        while (val <= maxV) {
            pot = CalculateBond(FFMethod, bond->kb, bond->b0, val);
            fprintf(potCurFile, "%-8.4lf %-8.4lf\n", val, pot);
            val += gap;
        }
    } else if (strcmp(type, "Angle") == 0) {
        minV = angle->minDeg;
        maxV = angle->maxDeg;
        val = minV;
        gap = 0.1;
        
        fprintf(potCurFile, "# %-8s %-8s\n", "deg", "energy");
        while (val <= maxV) {
            // NULL could be assigned. For the future update.
            pot = CalculateAngle(FFMethod, angle->ktheta, angle->theta0, angle->kub, angle->r0, NULL, val);
            fprintf(potCurFile, "%-8.4lf %-8.4lf\n", val, pot);
            val += gap;
        }
    } else if (strcmp(type, "Dihedral") == 0) {
        minV = 0;
        maxV = 180;
        val = minV;
        gap = 1;

        fprintf(potCurFile, "# %-8s %-8s %-8s\n", "dis", "deg", "energy");
        while (val <= maxV) {
            dis = CalculateDihedralDis(dihedral->data.bond, 
                                       dihedral->data.angle,
                                       val);
            pot = CalculateDihedral(FFMethod, dihedral->kphi, dihedral->n, dihedral->delta, val);
            fprintf(potCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", dis, val, pot);
            val += gap;
        }
    } else if (strcmp(type, "Improper") == 0) {
        minV = improper->list[0].val;
        maxV = improper->list[1].val;
        val = minV;
        gap = 0.1;
        
        fprintf(potCurFile, "# %-8s %-8s\n", "deg", "energy");
        while (val <= maxV) {
            pot = CalculateImproper(FFMethod, improper->komega, improper->omega0, val);
            fprintf(potCurFile, "%-8.4lf %-8.4lf\n", val, pot);
            val += gap;
        }
    }
    
    /*
    else if (strcmp(type, "LJ") == 0) {
        double parameter[2];
        
        for (int i = 1; i <= FFInfo.bead[0].val[0]; i ++) {
            if (FFInfo.LJ[i].epsilon == 0) continue;
            
            for (int j = i; j <= FFInfo.bead[0].val[0]; j ++) {
                if (FFInfo.LJ[j].epsilon == 0) continue;
                
                sprintf(name, "%s/%s-%s.txt", savedir, FFInfo.LJ[i].name, FFInfo.LJ[j].name);
                potCurFile = fopen(name, "w");
                if (potCurFile == NULL) {
                    printf("!!ERROR!!: Cannot create PotCur file in the directory %s!\n", name);
                    printf("           %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                
                parameter[0] = sqrt(FFInfo.LJ[i].epsilon * FFInfo.LJ[j].epsilon);
                parameter[1] = 0.5 * (FFInfo.LJ[i].sigma + FFInfo.LJ[j].sigma);
                
                dis = FindMaxPotDis(FFMethod,
                                    parameter,
                                    0.6);

                while (dis < MAX_R) {
                    poten = CalculateLJ(FFMethod, parameter[0], parameter[1], dis);
                    fprintf(potCurFile, "%8.4lf %8.4lf\n", dis, poten);
                    dis += gap;
                }

                fclose(potCurFile);
            }
        }
    }
     */
    
    fclose(potCurFile);
    return;
}


void PrintDisPot(char *type, char *atom1, char *atom2,
                 struct FFBondStr     *bond,
                 struct FFAngleStr    *angle,
                 struct FFDihedralStr *dihedral,
                 struct FFImproperStr *improper) {
    FILE *disPotCurFile;
    char name[PATH_LEN];
    
    sprintf(name, "%s/Debug%s-%s-%s_disPot.txt",
            savedir,
            type,
            atom1,
            atom2);
    disPotCurFile = fopen(name, "w");
    
    if (strcmp(type, "Bond") == 0) {
        fprintf(disPotCurFile, "# %-8s %-8s\n", "dis", "energy");
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf\n", bond->dmin, 0.4);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf\n", bond->dmin, 0.0);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf\n", bond->dmax, 0.0);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf\n", bond->dmax, 0.4);
    } else if (strcmp(type, "Angle") == 0) {
        fprintf(disPotCurFile, "# %-8s %-8s %-8s\n", "dis", "deg", "energy");
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", angle->dmin, angle->minDeg, 0.4);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", angle->dmin, angle->minDeg, 0.0);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", angle->dmax, angle->maxDeg, 0.0);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", angle->dmax, angle->maxDeg, 0.4);
    } else if (strcmp(type, "Dihedral") == 0) {
        fprintf(disPotCurFile, "# %-8s %-8s %-8s\n", "dis", "deg", "energy");
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n",
                                dihedral->list[0].dis,
                                0.0,
                                0.6);
        for (int i = 1; i <= DHMAX_N + 1; i ++) {
            fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n",
                                    dihedral->list[i - 1].dis,
                                    dihedral->list[i].val,
                                    dihedral->list[i].pot);
            fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n",
                                    dihedral->list[i].dis,
                                    dihedral->list[i].val,
                                    dihedral->list[i].pot);
        }
    } else if (strcmp(type, "Improper") == 0) {
        fprintf(disPotCurFile, "# %-8s %-8s %-8s\n", "dis", "deg", "energy");
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", improper->list[0].dis, improper->list[0].val, 0.4);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", improper->list[0].dis, improper->list[0].val, 0.0);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", improper->list[1].dis, improper->list[1].val, 0.0);
        fprintf(disPotCurFile, "%-8.4lf %-8.4lf %-8.4lf\n", improper->list[1].dis, improper->list[1].val, 0.4);
    }
    
    fclose(disPotCurFile);

    /*
    else if (strcmp(type, "LJ") == 0) {
        sprintf(name, "%s/%s-%s_disPot.txt",
                savedir,
                atom1,
                atom2);
        disPotCurFile = fopen(name, "w");
        
        for (int j = 0; j < step; j ++) {
            fprintf(disPotCurFile, "%8.4lf %8.4lf\n", pot[j][0], pot[j][1]);
        }
        fclose(disPotCurFile);
    }
     */

    return;
}


double ConvertDeg2Rad(double value) {
    return value / 180 * M_PI;
}


double ConvertRad2Deg(double value) {
    return value * 180 / M_PI;
}


//this function can be only used for crystalized/ordered surface materials, like graphene.
void BuildSurfaceConnection(void) {
    int cellindex, pair;
    int *atom_CellAxis, cell_neighbor[4];
    int type[2], startNum = 0;
    double ratio[2] = {0.9, 1.1};
    double dis[5] = {0, INFINIT, INFINIT, INFINIT, INFINIT};
    double position[2][4] = {0};
    double positionshift[4] = {0};
    
    if (strcmp(surfInfo.style, "explicit") != 0)
        return;
    
    if (strcmp(surfInfo.type, "graphene") == 0) {
        while (strcmp(atom[++ startNum].property->resName, "GRA") != 0);
    } else {
        return;
    }
    
    PBC("full", NULL, &atomList);
    LinkList("full", NULL, NULL, &atomList);
    
    /*
    if (strcmp(surfInfo.type, "graphene") == 0) {
        dis[1] = 1.418; //    bond distance
        dis[2] = 2.456; //   angle distance
        dis[3] = 2.836; //dihedral distance
        dis[4] = 3.758; //dihedral distance
    } else {
    }
     */
    
    //find the 4 smallest distance between atoms
    for (int i = startNum; i <= atomnum; i ++) {
        POINT_TO_STRUCT(atom_CellAxis, atom[i].dynamic->cellIndex);
        
        for (int n = 0; n < atomList.stencilNum; n ++) {
            cell_neighbor[3] = atom_CellAxis[3] + atomList.z[n];
            cell_neighbor[2] = atom_CellAxis[2] + atomList.y[n];
            cell_neighbor[1] = atom_CellAxis[1] + atomList.x[n]; //scan the neighborhood 27 subcells, include the target subcell itself
            
            for (int j = 1; j <= 3; j ++) {
                positionshift[j] = 0; //PBC position shift
                if (cell_neighbor[j] < 0) {
                    cell_neighbor[j] += cellnum[j];
                    positionshift[j] = -1 * boxDimension[j];
                } else if (cell_neighbor[j] >= cellnum[j]) {
                    cell_neighbor[j] -= cellnum[j];
                    positionshift[j] = boxDimension[j];
                }
            }
            
            cellindex = cell_neighbor[3] * cellnum[1] * cellnum[2] +
                        cell_neighbor[2] * cellnum[1] +
                        cell_neighbor[1] + 1;
            
            pair = celllist[atomnum + cellindex];
            while (pair != 0) {
                if (pair > i) {
                    DOT_PLUS(atom[pair].dynamic->coordinate, positionshift, position[0]);
                    DOT_MINUS(atom[i].dynamic->coordinate,   position[0],   position[1]);
                    dis[0] = DOT_PROD(position[1], position[1]);
                    dis[0] = sqrt(dis[0]);
                    
                    if (dis[1] > dis[0] * ratio[1]) {
                        dis[2] = dis[1];
                        dis[1] = dis[0];
                    } else if (dis[0] > dis[1] * ratio[1] &&
                               dis[2] > dis[0] * ratio[1]) {
                        dis[3] = dis[2];
                        dis[2] = dis[0];
                    } else if (dis[0] > dis[2] * ratio[1] &&
                               dis[3] > dis[0] * ratio[1]) {
                        dis[4] = dis[3];
                        dis[3] = dis[0];
                    } else if (dis[0] > dis[3] * ratio[1] &&
                               dis[4] > dis[0] * ratio[1]) {
                        dis[4] = dis[0];
                    }
                }
                pair = celllist[pair];
            }
        }
    }
    
    //start to build the bonds, angles and dihedral angles
    for (int i = startNum; i <= atomnum; i ++) {
        POINT_TO_STRUCT(atom_CellAxis, atom[i].dynamic->cellIndex);
        
        for (int n = 0; n < atomList.stencilNum; n ++) {
            cell_neighbor[3] = atom_CellAxis[3] + atomList.z[n];
            cell_neighbor[2] = atom_CellAxis[2] + atomList.y[n];
            cell_neighbor[1] = atom_CellAxis[1] + atomList.x[n]; //scan the neighborhood 27 subcells, include the target subcell itself
            
            for (int j = 1; j <= 3; j ++) {
                positionshift[j] = 0; //PBC position shift
                if (cell_neighbor[j] < 0) {
                    cell_neighbor[j] += cellnum[j];
                    positionshift[j] = -1 * boxDimension[j];
                } else if (cell_neighbor[j] >= cellnum[j]) {
                    cell_neighbor[j] -= cellnum[j];
                    positionshift[j] = boxDimension[j];
                }
            }
            
            cellindex = cell_neighbor[3] * cellnum[1] * cellnum[2] +
                        cell_neighbor[2] * cellnum[1] +
                        cell_neighbor[1] + 1;
            
            pair = celllist[atomnum + cellindex];
            while (pair != 0) {
                if (pair > i) {
                    
                    DOT_PLUS(atom[pair].dynamic->coordinate, positionshift, position[0]);
                    DOT_MINUS(atom[i].dynamic->coordinate,   position[0],   position[1]);
                    dis[0] = DOT_PROD(position[1], position[1]);
                    dis[0] = sqrt(dis[0]);
                    
                    type[0] = FindAtomIndex(atom[pair].property->typeName[1]);
                    type[1] = FindAtomIndex(atom[  i ].property->typeName[1]);
                    
                    if (dis[0] > dis[1] * ratio[0] &&
                        dis[0] < dis[1] * ratio[1]) {
                        FindNextBondList(pair, atom[i].property,
                                         FFInfo.bond[type[0]][type[1]].dmin,
                                         FFInfo.bond[type[0]][type[1]].dmax);
                        connectionMap[i][pair] |= BOND_CONNECT;
                        
                    } else if (dis[0] > dis[2] * ratio[0] &&
                               dis[0] < dis[2] * ratio[1]) {
                        FindNextBondList(pair, atom[i].property,
                                         FFInfo.bond[type[0]][type[1]].dmin,
                                         FFInfo.bond[type[0]][type[1]].dmax);
                        connectionMap[i][pair] |= BOND_CONNECT;

                    } else if ((dis[0] > dis[3] * ratio[0] && dis[0] < dis[3] * ratio[1]) ||
                               (dis[0] > dis[4] * ratio[0] && dis[0] < dis[4] * ratio[1])) {
                        FindNextConstrList(pair, atom[i].property,
                                           FFInfo.dihedral[type[0]][type[1]].data.step,
                                           FFInfo.dihedral[type[0]][type[1]].data.dataString);
                        connectionMap[i][pair] |= CONSTRAINT_CONNECT;
                    }
                }
                pair = celllist[pair];
            }
        }
        
        BondDuplicate(atom[i].property->sequence.molNum);
        ConstrDuplicate(atom[i].property->sequence.molNum);
    }
    
    free(celllist); celllist = NULL; //so that the next LinkList call can be valid.
    return;
}


void SortList(struct DisPotStr *list) {
    struct DisPotStr tmp;
    for (int i = 0; i < DHMAX_N; i ++) {
        for (int j = i + 1; j < DHMAX_N; j ++) {
            if (list[i].dis > list[j].dis) {
                tmp = list[i];
                list[i] = list[j];
                list[j] = tmp;
            }
        }
    }
    
    return;
}


double CalculateAngleCosineLaw(double a, double b, double c) {
    double val = (a*a + b*b - c*c) / (2 * a * b);
    return ConvertRad2Deg(cosh(val));
}
