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


void HBInfo(int id) {
    int connectAtom = 0;
    int alphaRecord = 0, betaRecord = 0;
    int beginP = 1, endP = numofprotein;
    int thisAtomNum = atomnum, offSet = 0;
    long **markAA = NULL, totalFrame;
    double phi = 0, psi = 0;
    double prePhi = 0, prePsi = 0;
    char directory[1024], buffer[1024];
    struct SectionStr *sect;
    FILE ***RamachFile = NULL, *HBOutputFile, *AAMarkOutputFile = NULL;
    FILE **ContactMapFile = NULL, **ContactMapAAFile = NULL;
    FILE *TrjInputFile, *CntInputFile;

    if (nPP) { //for only analyze a specific peptide instead of all
        beginP = endP = nPP;
        thisAtomNum = protein[nPP].endAtomNum - protein[nPP].startAtomNum + 1;
        offSet = protein[nPP].startAtomNum - 1;
    }
    
    printf("Analyzing HB information... ");
    fflush(stdout);
    
    sprintf(directory, "%s%s", path, files[outHBInfo][id].name);
    HBOutputFile = fopen(directory, "w");
    fprintf(HBOutputFile, "# %13s %17s %17s %17s %17s %17s %17s %17s\n",
            "Step",
            "alpha",
            "3_10",
            "pi",
            "beta",
            "SS",
            "other",
            "total");
    
    if (numofprotein == 1) {
        sprintf(directory, "%s%s", path, files[outAAMark][id].name);
        AAMarkOutputFile = fopen(directory, "w");
        fprintf(AAMarkOutputFile, "@    xaxis label \"AANum\"\n");
        fprintf(AAMarkOutputFile, "@    yaxis label \"Mark\"\n");
        fprintf(AAMarkOutputFile, "@    s0 legend \"\\f{Symbol}a\"\n");
        fprintf(AAMarkOutputFile, "@    s1 legend \"\\f{Symbol}b\"\n");
        fprintf(AAMarkOutputFile, "@TYPE xy\n");
        
        int num = protein[1].endAANum - protein[1].startAANum + 1;
        LONG_2CALLOC(markAA, (num + 1), 2);
    }
    
    if (analysisList[AnalyzeRamach]) {
        RamachFile = (FILE ***)calloc(numofprotein, sizeof(FILE **));
        for (int n = 0; n < numofprotein; n ++) {
            RamachFile[n] = (FILE **)calloc(21, sizeof(FILE *));
        }

        for (int n = beginP - 1; n < endP; n ++) {
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
        ContactMapFile = (FILE **)calloc(numofprotein, sizeof(FILE *));
        ContactMapAAFile = (FILE **)calloc(numofprotein, sizeof(FILE *));
        
        for (int i = beginP - 1; i < endP; i ++) {
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
                for (int n = beginP - 1; n < endP; n ++) {
                    fprintf(ContactMapFile[n],   "Frame = %li\n", step * (long)sect->outputRate);
                    fprintf(ContactMapAAFile[n], "Frame = %li\n", step * (long)sect->outputRate);
                }
            }
            
            alphaRecord = betaRecord = 0;
            
            for (int i = offSet + 1; i <= offSet + thisAtomNum; i++) {
                
                if (atom[i].property->typeofAtom == 18) {
                    RamachandranPlot(i, &phi, &psi);
                    
                    if (phi != 0 && psi != 0 && analysisList[AnalyzeRamach]) {
                        fprintf(RamachFile[atom[i].property->sequence.proteinNum - 1][FindAANum(atom[i].property->nameofAA)], "%-10.2lf%-10.2lf\n", phi, psi);
                    }
                }
                
                connectAtom = CheckHBConnection(i);
                if (connectAtom > 0 ) {
					if (connectAtom > i) { //only count once
						HBSum.total++;
						if (atom[i].property->typeofAtom == 31 && atom[connectAtom].property->typeofAtom == 31)
							HBSum.SS++;
					}

					if (atom[i].property->typeofAtom == 2 || atom[i].property->typeofAtom == 4) {
						CheckHelixHB(i, connectAtom, phi, psi);
						CheckBetaHB(i, connectAtom, phi, psi);

						prePsi = psi;
						prePhi = phi;

						if (numofprotein == 1 && (alphaRecord != HBSum.helix_alpha || betaRecord != HBSum.beta)) {
							if (alphaRecord != HBSum.helix_alpha) {
								markAA[atom[i].property->sequence.aminoacidNum][0] ++;
								markAA[atom[connectAtom].property->sequence.aminoacidNum][0] ++;
							} else {
								markAA[atom[i].property->sequence.aminoacidNum][1] ++;
								markAA[atom[connectAtom].property->sequence.aminoacidNum][1] ++;
							}

							alphaRecord = HBSum.helix_alpha;
							betaRecord = HBSum.beta;
						}
					}
                }
                
                if (analysisList[AnalyzeConMap] &&
                    atom[i].property->sequence.proteinNum == atom[connectAtom].property->sequence.proteinNum) {
                    fprintf(ContactMapFile[atom[i].property->sequence.proteinNum - 1], "%8i%8i%10s%10s\n", i, connectAtom, atom[i].property->nameofAA, atom[connectAtom].property->nameofAA);
                    fprintf(ContactMapAAFile[atom[i].property->sequence.proteinNum - 1], "%8i%8i\n", atom[i].property->sequence.aminoacidNum, atom[connectAtom].property->sequence.aminoacidNum);
                }
            }
            
            SaveHBInfo(sect, HBOutputFile);
        }
        fclose(CntInputFile);
    }
    
    if (numofprotein == 1) {
        for (int i = 1; i <= protein[1].endAANum - protein[1].startAANum + 1; i ++) {
            fprintf(AAMarkOutputFile, "%8i%8li%8li\n", i, markAA[i][0], markAA[i][1]);
        }
    }

    fclose(TrjInputFile);
    fclose(HBOutputFile);
    
    if (numofprotein == 1) {
        fclose(AAMarkOutputFile);
        _2FREE(markAA);
    }
    
    if (analysisList[AnalyzeRamach]) {
        for (int n = beginP - 1; n < endP; n ++) {
            for (int i = 0; i < 21; i ++) {
                fclose(RamachFile[n][i]);
            }
        }
        
        for (int n = 0; n < numofprotein; n ++) {
            free(RamachFile[n]);
        }

        free(RamachFile); RamachFile = NULL;
    }

    if (analysisList[AnalyzeConMap]) {
        for (int i = beginP - 1; i < endP; i ++) {
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
    int flag_phi=1, flag_psi=1;
    int N = atom_N, Ca, C, C_1, N_1;
    double position_N[4]={0}, position_C[4];
    double position_Ca[4]={0};
    double position_N1[4]={0};
    double position_C1[4]={0};
    double b1[4]={0},b2[4]={0},b3[4]={0};
    double tempValue1, tempValue2;
    double tempVector1[4], tempVector2[4];
    double x,y;
    double thispi=57.29578;
    
    TRANSFER_VECTOR(position_N, atom[N].dynamic->coordinate);
    
    j = N;
    while (j > 1 && atom[--j].property->typeofAtom != 6)
        ;
    if (j != N && atom[j].property->sequence.proteinNum == atom[N].property->sequence.proteinNum) {
        C_1 = j;
        TRANSFER_VECTOR(position_C1, atom[C_1].dynamic->coordinate);
    } else {
        C_1 = 0;
        flag_phi = 0;
        *phi = 0;
    }
    
    j = N;
    while (strcmp(atom[++j].property->name, "CA"))
        ;
    if (atom[j].property->sequence.proteinNum == atom[N].property->sequence.proteinNum) {
        Ca = j;
        TRANSFER_VECTOR(position_Ca, atom[Ca].dynamic->coordinate);
    }
    
    j = N;
    while (atom[++j].property->typeofAtom != 6)
        ;
    if (atom[j].property->sequence.proteinNum == atom[N].property->sequence.proteinNum) {
        C = j;
        TRANSFER_VECTOR(position_C, atom[C].dynamic->coordinate);
    }
    
    j = N;
    while (j < atomnum && atom[++j].property->typeofAtom != 18)
        ;
    if (atom[j].property->sequence.proteinNum == atom[N].property->sequence.proteinNum) {
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
    
    for (int i = 1; i <= atomnum; i ++) {
        if (strcmp(atom[i].property->nameofAA, "SOL") == 0) //would change later
            continue;

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
    
    if (!(atom[atomNum].property->typeofAtom == 2 || atom[atomNum].property->typeofAtom == 4)) {
        printf("!!ERROR!!: the target atom type has to be HB (backbone H)! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (connectAtom > 0 && phi != 0 &&
        atom[atomNum].property->sequence.proteinNum == atom[connectAtom].property->sequence.proteinNum) {
        gapResidue = atom[atomNum].property->sequence.aminoacidNum - atom[connectAtom].property->sequence.aminoacidNum;
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
    if (!(atom[atomNum].property->typeofAtom == 2 || atom[atomNum].property->typeofAtom == 4)) {
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
    
    return;
}

void SaveHBInfo(struct SectionStr *sect, FILE *HBInfoFile) {
    HBSum.other = HBSum.total    - HBSum.helix_alpha - HBSum.helix_310 -
                  HBSum.helix_pi - HBSum.beta        - HBSum.SS;
    
    fprintf(HBInfoFile, "%15.2lf %17i %17i %17i %17i %17i %17i %17i\n",
            step * sect->outputRate,
            HBSum.helix_alpha,
            HBSum.helix_310,
            HBSum.helix_pi,
            HBSum.beta,
            HBSum.SS,
            HBSum.other,
            HBSum.total);
}

int ReadHBFile(FILE *HBInputFile, struct HBReadStr *thisHB) {
    char buffer[1024];
    
repeat:
    if (fgets(buffer, sizeof(buffer), HBInputFile) == NULL)
        return -1;
    
    if (buffer[0] == '@' || buffer[0] == '#') {
        goto repeat;
    }
    
    sscanf(buffer, "%lf%i%i%i%i%i%i%i",
           &thisHB->step,
           &thisHB->type.helix_alpha,
           &thisHB->type.helix_310,
           &thisHB->type.helix_pi,
           &thisHB->type.beta,
           &thisHB->type.SS,
           &thisHB->type.other,
           &thisHB->type.total);
    
    return 0;
}
