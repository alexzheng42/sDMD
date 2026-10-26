//
//  HBRamach.c
//  Analysis
//
//  Created by Size Zheng on 10/31/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "Analysis.h"

int **connectionMap;
char AAName[21][4] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GEL", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
struct HBType HBSum;
struct HBReadStr hb;

static long step;

static void RamachandranPlot(int atom_N, double *phi, double *psi);
static void ResetHBNum(void);
static void CheckHelixHB(int atomNum, int connectAtom, double phi, double psi);
static void CheckBetaHB(int atomNum, int connectAtom, double phi, double psi);
static void SaveHBInfo(struct SectionStr *sect, FILE *HBInfoFile);
static int FindAANum(char * name);
static int CheckAAHBStyle(struct AtomStr *atom_i, struct AtomStr *atom_j);


void HBInfo(int id) {
    int connectAtom = 0;
    int alphaRecord = 0, betaRecord = 0;
    long **markAA = NULL, totalFrame;
    double phi = 0, psi = 0;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE ***RamachFile = NULL, *HBOutputFile, *AAMarkOutputFile = NULL;
    FILE **ContactMapFile = NULL, **ContactMapAAFile = NULL;
    FILE *TrjInputFile, *CntInputFile;
   
    printf("Analyzing HB information... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outHBInfo][id].name);
    HBOutputFile = fopen(directory, "w");
    fprintf(HBOutputFile, "# %13s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s\n",
            "Step",
            "alpha",
            "3_10",
            "pi",
            "beta",
            "SS",
            "other",
            "total",
            "intraBB",
            "intraBS",
            "intraSS",
            "interBB",
            "interBS",
            "interSS");
    
    if (numofMol == 1) {
        sprintf(directory, "%s%s", path, files[outAAMark][id].name);
        AAMarkOutputFile = fopen(directory, "w");
        fprintf(AAMarkOutputFile, "@    xaxis label \"AANum\"\n");
        fprintf(AAMarkOutputFile, "@    yaxis label \"Mark\"\n");
        fprintf(AAMarkOutputFile, "@    s0 legend \"\\f{Symbol}a\"\n");
        fprintf(AAMarkOutputFile, "@    s1 legend \"\\f{Symbol}b\"\n");
        fprintf(AAMarkOutputFile, "@TYPE xy\n");
        
        int num = molecule[1].endAANum - molecule[1].startAANum + 1;
        _CALLOC2(long, markAA, (num + 1), 2);
    }
    
    if (analysisList[AnalyzeRamach]) {
        RamachFile = (FILE ***)calloc(numofMol, sizeof(FILE **));
        for (int n = 0; n < numofMol; n ++) {
            RamachFile[n] = (FILE **)calloc(21, sizeof(FILE *));
        }

        for (int m = 1; m <= group.molNum; m ++) {
            int n = group.list_Mol[m] - 1;
            for (int i = 0; i < 21; i ++) {
                sprintf(directory, "%sRamachFile_p%i_%s.txt", path, n + 1, AAName[i]);
                RamachFile[n][i] = fopen(directory, "w");
                fprintf(RamachFile[n][i], "%s\n", "@    view 0.150000, 0.150000, 0.850000, 0.850000");
                fprintf(RamachFile[n][i], "%s\n", "@    xaxis  label \"\\xF\"");
                fprintf(RamachFile[n][i], "%s\n", "@    yaxis  label \"\\xY\"");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 type xy");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 symbol 1");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 symbol size 0.100000");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 symbol color 1");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 symbol pattern 1");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 symbol fill color 1");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 symbol fill pattern 1");
                fprintf(RamachFile[n][i], "%s\n", "@    s0 line type 0");
            }
        }
    }
    
    if (analysisList[AnalyzeConMap]) {
        ContactMapFile = (FILE **)calloc(numofMol, sizeof(FILE *));
        ContactMapAAFile = (FILE **)calloc(numofMol, sizeof(FILE *));
        
        for (int m = 1; m <= group.molNum; m ++) {
            int i = group.list_Mol[m] - 1;
            sprintf(directory, "%sContactInfo_p%i.txt", path, i + 1);
            ContactMapFile[i] = fopen(directory, "w");
            
            sprintf(directory, "%sContactAAInfo_p%i.txt", path, i + 1);
            ContactMapAAFile[i] = fopen(directory, "w");
        }
    }
    
    sprintf(directory, "%s%s", path, files[outTrj][id].name);
    TrjInputFile = fopen(directory, "r");
    if (TrjInputFile == NULL) {
        printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[outTrj][id].name, path, __FILE__, __LINE__);
        printf("           You may need to execute -rPBC first to generate the required file!\n");
        exit(EXIT_FAILURE);
    }
    
    for (int sectNum = 0; sectNum < fileList[id].count; sectNum ++) {
        memset(buffer, '\0', sizeof(buffer));
        FindTargetFile(files[inCnt][id].name, fileList[id].list[sectNum + 1], buffer);
        
        sprintf(directory, "%s%s", path, buffer);
        CntInputFile = fopen(directory, "r");
        if (CntInputFile == NULL) {
            printf("!!ERROR!!: cannot find file %s in directory %s. make sure the input path is correct and the file does exist! %s:%i\n", files[inCnt][id].name, path, __FILE__, __LINE__);
            printf("           You may need to specify the connection information file by using -cnt flag!\n");
            exit(EXIT_FAILURE);
        }
        
        totalFrame = (sectInfo[sectNum].oldTime + sectInfo[sectNum].frameCount) / sectInfo[sectNum].outputRate + 1;
        for (step = sectInfo[sectNum].oldTime / sectInfo[sectNum].outputRate; step < totalFrame; step ++) {
            if (ReadGro(TrjInputFile) || ReadConnectionMap(CntInputFile)) break;
            
            PrintProcess(step);
            
            sect = &sectInfo[sectNum];
            ResetHBNum();
            if (analysisList[AnalyzeConMap]) {
                for (int m = 1; m <= group.molNum; m ++) {
                    int n = group.list_Mol[m] - 1;
                    fprintf(ContactMapFile[n],   "Frame = %li\n", step * (long)sect->outputRate);
                    fprintf(ContactMapAAFile[n], "Frame = %li\n", step * (long)sect->outputRate);
                }
            }
            
            alphaRecord = betaRecord = 0;
            
            for (int m = 1; m <= group.atomNum; m ++) {
                int i = group.list_Atom[m];
                
                if (strcmp(atom[i].property->typeName[1], "NZB") == 0) {
                    RamachandranPlot(i, &phi, &psi);
                    
                    if (phi != 0 && psi != 0 && analysisList[AnalyzeRamach]) {
                        fprintf(RamachFile[atom[i].property->sequence.molNum - 1][FindAANum(atom[i].property->resName)], "%-10.2lf%-10.2lf\n", phi, psi);
                    }
                }
                
                connectAtom = CheckHBConnection(i);
                if (connectAtom > 0 ) {
					if (connectAtom > i) { //only count once
						HBSum.total++;
						if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 4)
							HBSum.SS++;
                        
                        //calculate inter/intra BB, BS, SS HB
                        if (atom[i].property->sequence.molNum != atom[connectAtom].property->sequence.molNum) {
                            if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 1)
                                HBSum.interBB ++;
                            else if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 2)
                                HBSum.interBS ++;
                            else if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 3)
                                HBSum.interSS ++;
                        } else {
                            if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 1)
                                HBSum.intraBB ++;
                            else if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 2)
                                HBSum.intraBS ++;
                            else if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 3)
                                HBSum.intraSS ++;
                        }
					}

					if (CheckAAHBStyle(&atom[i], &atom[connectAtom]) == 1) {
						CheckHelixHB(i, connectAtom, phi, psi);
						CheckBetaHB(i, connectAtom, phi, psi);

						if (numofMol == 1 && (alphaRecord != HBSum.helix_alpha || betaRecord != HBSum.beta)) {
							if (alphaRecord != HBSum.helix_alpha) {
								markAA[atom[i].property->sequence.resNum][0] ++;
								markAA[atom[connectAtom].property->sequence.resNum][0] ++;
							} else {
								markAA[atom[i].property->sequence.resNum][1] ++;
								markAA[atom[connectAtom].property->sequence.resNum][1] ++;
							}

							alphaRecord = HBSum.helix_alpha;
							betaRecord = HBSum.beta;
						}
					}
                }
                
                if (analysisList[AnalyzeConMap] &&
                    atom[i].property->sequence.molNum == atom[connectAtom].property->sequence.molNum) {
                    fprintf(ContactMapFile[atom[i].property->sequence.molNum - 1], "%8i%8i%10s%10s\n", i, connectAtom, atom[i].property->resName, atom[connectAtom].property->resName);
                    fprintf(ContactMapAAFile[atom[i].property->sequence.molNum - 1], "%8i%8i\n", atom[i].property->sequence.resNum, atom[connectAtom].property->sequence.resNum);
                }
            }
            
            SaveHBInfo(sect, HBOutputFile);
        }
        fclose(CntInputFile);
    }
    
    if (numofMol == 1) {
        for (int i = 1; i <= molecule[1].endAANum - molecule[1].startAANum + 1; i ++) {
            fprintf(AAMarkOutputFile, "%8i%8li%8li\n", i, markAA[i][0], markAA[i][1]);
        }
    }

    fclose(TrjInputFile);
    fclose(HBOutputFile);
    
    if (numofMol == 1) {
        fclose(AAMarkOutputFile);
        _FREE2(markAA);
    }
    
    if (analysisList[AnalyzeRamach]) {
        for (int m = 1; m <= group.molNum; m ++) {
            int n = group.list_Mol[m] - 1;
            for (int i = 0; i < 21; i ++) {
                fclose(RamachFile[n][i]);
            }
        }
        
        for (int n = 0; n < numofMol; n ++) {
            free(RamachFile[n]);
        }

        free(RamachFile); RamachFile = NULL;
    }

    if (analysisList[AnalyzeConMap]) {
        for (int m = 1; m <= group.molNum; m ++) {
            int i = group.list_Mol[m] - 1;
            fclose(ContactMapFile[i]);
            fclose(ContactMapAAFile[i]);
        }
        free(ContactMapFile);   ContactMapFile = NULL;
        free(ContactMapAAFile); ContactMapAAFile = NULL;
    }
    
    printf("\bDone!\n");
    fflush(stdout);

    return;
}


void RamachandranPlot(int atom_N, double *phi, double *psi)
{
    int j;
    int flag_phi = 1, flag_psi = 1;
    int N = atom_N, Ca, C, C_1, N_1;
    double position_N[4]  = {0};
    double position_C[4]  = {0};
    double position_Ca[4] = {0};
    double position_N1[4] = {0};
    double position_C1[4] = {0};
    double b1[4] = {0}, b2[4] = {0}, b3[4] = {0};
    double tempValue1, tempValue2;
    double tempVector1[4], tempVector2[4];
    double x,y;
    double thispi=57.29578;
    
    TRANSFER_VECTOR(position_N, atom[N].dynamic->coordinate);
    
    j = N;
    while (j > 1 && strcmp(atom[--j].property->name, "C") != 0)
        ;
    if (j != N && atom[j].property->sequence.molNum == atom[N].property->sequence.molNum) {
        C_1 = j;
        TRANSFER_VECTOR(position_C1, atom[C_1].dynamic->coordinate);
    } else {
        C_1 = 0;
        flag_phi = 0;
        *phi = 0;
    }
    
    j = N;
    while (j < atomnum && strcmp(atom[++j].property->name, "CA") != 0)
        ;
    if (atom[j].property->sequence.molNum == atom[N].property->sequence.molNum) {
        Ca = j;
        TRANSFER_VECTOR(position_Ca, atom[Ca].dynamic->coordinate);
    }
    
    j = N;
    while (j < atomnum && strcmp(atom[++j].property->name, "C") != 0)
        ;
    if (atom[j].property->sequence.molNum == atom[N].property->sequence.molNum) {
        C = j;
        TRANSFER_VECTOR(position_C, atom[C].dynamic->coordinate);
    }
    
    j = N;
    while (j < atomnum && strcmp(atom[++j].property->name, "N") != 0)
        ;
    if (atom[j].property->sequence.molNum == atom[N].property->sequence.molNum) {
        N_1 = j;
        TRANSFER_VECTOR(position_N1, atom[N_1].dynamic->coordinate)
    } else {
        N_1 = 0;
        flag_psi = 0;
        *psi = 0;
    }
    
    //--------------calculate phi and psi--------------
    if (flag_phi == 1 && flag_psi == 1) {
        DOT_MINUS(position_N, position_C1, b1);
        DOT_MINUS(position_Ca, position_N, b2);
        DOT_MINUS(position_C, position_Ca, b3);
        
        tempValue1 = absvalue_vector(b2);
        tempValue1 *= tempValue1;
        tempValue2 = DOT_PROD(b1, b3);
        
        x = -1 * tempValue1 * tempValue2;
        
        tempValue1 = DOT_PROD(b1, b2);
        tempValue2 = DOT_PROD(b2, b3);
        
        x += tempValue1 * tempValue2;
        
        
        tempValue1 = absvalue_vector(b2);
        FACTOR_PROD(tempValue1, b1, tempVector1);
        CROSS_PROD(b2, b3, tempVector2);
        y = DOT_PROD(tempVector1, tempVector2);
        
        *phi = atan2(y, x) * thispi;
        
        
        DOT_MINUS(position_Ca, position_N, b1);
        DOT_MINUS(position_C, position_Ca, b2);
        DOT_MINUS(position_N1, position_C, b3);
        
        tempValue1 = absvalue_vector(b2);
        tempValue1 *= tempValue1;
        tempValue2 = DOT_PROD(b1, b3);
        
        x = -1 * tempValue1 * tempValue2;
        
        tempValue1 = DOT_PROD(b1, b2);
        tempValue2 = DOT_PROD(b2, b3);
        
        x += tempValue1 * tempValue2;
        
        
        tempValue1 = absvalue_vector(b2);
        FACTOR_PROD(tempValue1, b1, tempVector1);
        CROSS_PROD(b2, b3, tempVector2);
        y = DOT_PROD(tempVector1, tempVector2);
        
        *psi = atan2(y, x) * thispi;
        
    }
    //--------------calculate phi and psi--------------
}

int FindAANum(char * name) {
    for (int i = 0; i < 21; i ++) {
        if (strcmp(name, AAName[i]) == 0) {
            return i;
        }
    }
    
    return -1;
}


int ReadConnectionMap(FILE *inputFile) {
    int num, temp;
    char buffer[1024];
    
    if (fscanf(inputFile, "%[^\n]\n", buffer) == EOF)
        return -1;
    
    // sDMD writes one line per atom, water included (water can form HBs, e.g. Library_TIP3P);
    // every line has to be read or the following atoms get the wrong connections
    for (int i = 1; i <= atomnum; i ++) {
        for (int n = 1; n <= atomnum; n ++) {
            connectionMap[i][n] = 0;
        }
        fgets(buffer, sizeof(buffer), inputFile);
        
        int curpos = 0, pos = 0;
        while (sscanf(buffer + curpos, "%i%i%n", &num, &temp, &pos) != EOF) {
            connectionMap[i][num] = temp;
            curpos += pos;
        }
    }
    
    return 0;
}


int CheckHBConnection(int thisAtomNum) {
    for (int i = 1; i <= atomnum; i ++) {
        if (connectionMap[thisAtomNum][i] & HBConnect) {
            return i;
        }
    }
    
    return 0;
}

void CheckHelixHB(int atomNum, int connectAtom, double phi, double psi) {
    int gapResidue = 0;
    //double dihedralSum;
    
    if (CheckAAHBStyle(&atom[atomNum], &atom[connectAtom]) != 1) {
        printf("!!ERROR!!: the target atom type has to be HB (backbone H)! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (connectAtom > 0 && phi != 0 &&
        atom[atomNum].property->sequence.molNum == atom[connectAtom].property->sequence.molNum) {
        gapResidue = atom[atomNum].property->sequence.resNum - atom[connectAtom].property->sequence.resNum;
        //dihedralSum = prePsi + phi;
        
        if (phi >= -180 && phi <= 0 &&
            psi >= -100 && psi <= 45) {
            if (/*dihedralSum >= -110 && dihedralSum <= -100 &&*/ gapResidue == 4) {
                HBSum.helix_alpha ++;
            } else if (/*dihedralSum >= -80 && dihedralSum <= -70 &&*/ gapResidue == 3) {
                HBSum.helix_310 ++;
            } else if (/*dihedralSum >= -135 && dihedralSum <= -125 &&*/ gapResidue == 5) {
                HBSum.helix_pi ++;
            }
        }
    }
    return;
}

void CheckBetaHB(int atomNum, int connectAtom, double phi, double psi) {
    if (CheckAAHBStyle(&atom[atomNum], &atom[connectAtom]) != 1) {
        printf("!!ERROR!!: the target atom type has to be HB (backbone H)! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (connectAtom > 0 && phi != 0) {
        if (phi >= -175 && phi <= -70 &&
            psi >=  90  && psi <= 170) {
            HBSum.beta ++;
        }
        /*
        if (phi >= -145 && phi <= -95 &&
            psi >=  100 && psi <= 155) {
            HBSum.beta_para ++;
        } else if (phi >= -150 && phi <= -95 &&
                   psi >=  110 && psi <= 160) {
            HBSum.beta_anti ++;
        }
         */
    }
    return;
}

void ResetHBNum(void) {
    HBSum.helix_alpha = 0;
    HBSum.helix_310   = 0;
    HBSum.helix_pi    = 0;
    
    HBSum.beta        = 0;
    /*
    HBSum.beta_anti = 0;
    HBSum.beta_para = 0;
     */
    HBSum.SS          = 0;
    HBSum.other       = 0;
    HBSum.total       = 0;
    HBSum.interBB     = 0;
    HBSum.interBS     = 0;
    HBSum.interSS     = 0;
    HBSum.intraBB     = 0;
    HBSum.intraBS     = 0;
    HBSum.intraSS     = 0;
    
    return;
}

void SaveHBInfo(struct SectionStr *sect, FILE *HBInfoFile) {
    HBSum.other = HBSum.total    - HBSum.helix_alpha - HBSum.helix_310 -
                  HBSum.helix_pi - HBSum.beta        - HBSum.SS;
    
    fprintf(HBInfoFile, "%15.2lf %17i %17i %17i %17i %17i %17i %17i %17i %17i %17i %17i %17i %17i\n",
            step * sect->outputRate,
            HBSum.helix_alpha,
            HBSum.helix_310,
            HBSum.helix_pi,
            HBSum.beta,
            HBSum.SS,
            HBSum.other,
            HBSum.total,
            HBSum.intraBB,
            HBSum.intraBS,
            HBSum.intraSS,
            HBSum.interBB,
            HBSum.interBS,
            HBSum.interSS);
}

int ReadHBFile(FILE *HBInputFile, struct HBReadStr *thisHB) {
    char buffer[1024];
    
repeat:
    if (fgets(buffer, sizeof(buffer), HBInputFile) == NULL)
        return -1;
    
    if (buffer[0] == '@' || buffer[0] == '#') {
        goto repeat;
    }
    
    sscanf(buffer, "%lf%i%i%i%i%i%i%i%i%i%i%i%i%i",
           &thisHB->step,
           &thisHB->type.helix_alpha,
           &thisHB->type.helix_310,
           &thisHB->type.helix_pi,
           &thisHB->type.beta,
           &thisHB->type.SS,
           &thisHB->type.other,
           &thisHB->type.total,
           &thisHB->type.intraBB,
           &thisHB->type.intraBS,
           &thisHB->type.intraSS,
           &thisHB->type.interBB,
           &thisHB->type.interBS,
           &thisHB->type.interSS);
    
    return 0;
}


int CheckAAHBStyle(struct AtomStr *atom_i, struct AtomStr *atom_j) {
    if ((strcmp(atom_i->property->name, "H") == 0 && strcmp(atom_j->property->name, "O") != 0) ||
        (strcmp(atom_i->property->name, "O") == 0 && strcmp(atom_j->property->name, "H") != 0)) {
        return 1; //BB
    } else if (strcmp(atom_i->property->name, "H") == 0 ||
               strcmp(atom_i->property->name, "O") == 0) {
        return 2; //BS
    } else if (strcmp(atom_i->property->name, "SG") == 0 &&
               strcmp(atom_i->property->name, "SG") == 0) {
        return 4; //SG-SG
    } else {
        return 3; //SS
    }
    
    return 0;
}
