//
//  CG.c
//  sDMD
//
//  Created by Size on 2019/3/12.
//  Copyright © 2019 Size Zheng. All rights reserved.
//

#include "DMD.h"

//#define CG_DEBUG
#define POT_RO      1 //ro in surface potential equation
#define MAX_R       10
#define MAX_P       0.8

struct ParameterStr {
    int coreorShell;
    int s_time;
    double r_temp[4], v_ij[4], r_ij[4], b_ij, v_2, r_2, d_ij2;
    double positionshift[4];
    double lowerPotential, upperPotential;
    double shortlimit2, longlimit2;
    double calibratedcutoffr;
    double position_i[4], speed_i[4];
    double position_j[4], speed_j[4];
};

void InitializeEvent(struct InteractionEventStr *event);
void CreateCGPotentialMatrix(void);
void ReadCGPotTable(void);
void ReadGoPara(void);
void PrintPotCur_Res(void);
void PrintPotCur_Pep(void);
void PrintStepPotFun(void);
void PrintDisPot(double pot[128][2], int step, char *sur_type, int num);
char *GetAAName(int AAType);
double FindMinDisDimen(int dim);
double CalculateTime(struct ParameterStr *parameters);
double CalculateSurPoten(struct CGBeadStr sur, struct CGBeadStr res, double dis);
double FindMaxPotDis(struct CGBeadStr sur, struct CGBeadStr res, double setP);
double FindMinPotDis(struct CGBeadStr sur, struct CGBeadStr res, int print);
double FindMaxRadius(struct CGBeadStr sur, struct CGBeadStr res, double minPDis, double diff);


void InitializeCG(void) {
    if (strcmp(CG.type[1], "alphaC") == 0) {
        CG.totalNum = totalAminoAcids;
        CG.list = (int *)calloc(totalAminoAcids + 1, sizeof(int));
        CG.thisCG = (struct CGStr *)calloc(totalAminoAcids + 1, sizeof(struct CGStr));
        
        //all the atoms in the residue are mapped on the alpha-C
        for (int i = 1; i <= totalAminoAcids; i ++) {
            for (int n = aminoacid[i].startAtomNum; n <= aminoacid[i].endAtomNum; n ++) {
                if (strcmp(atom[n].property->name, "CA") == 0) {
                    CG.list[i]         = n;
                    CG.thisCG[i].start = aminoacid[i].startAtomNum;
                    CG.thisCG[i].end   = aminoacid[i].endAtomNum;
                    sprintf(atom[n].property->extraProperty[0], "CG");
                    break;
                }
            }
        }
    } else if (strcmp(CG.type[1], "residue") == 0) {
        CG.totalNum = totalAminoAcids;
        CG.list = (int *)calloc(totalAminoAcids + 1, sizeof(int));
        CG.thisCG = (struct CGStr *)calloc(totalAminoAcids + 1, sizeof(struct CGStr));
        
        //all the atoms in the residue share the same Chi of its residue
        for (int i = 1; i <= totalAminoAcids; i ++) {
            for (int n = aminoacid[i].startAtomNum; n <= aminoacid[i].endAtomNum; n ++) {
                if (strcmp(atom[n].property->name, "CA") == 0) {
                    CG.list[i]         = n;
                    CG.thisCG[i].start = aminoacid[i].startAtomNum;
                    CG.thisCG[i].end   = aminoacid[i].endAtomNum;
                }
                sprintf(atom[n].property->extraProperty[0], "CG");
            }
        }
    } else {
        printf("!!ERROR!!: the CG type is invalid. %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
	InitializeCGPotentialMatrix();
    return;
}


void InitializeCGPotentialMatrix(void) {
    char buffer[512], statement[128];
    FILE *input_file;

    if (strcmp(CG.type[1], "alphaC") == 0 ||
        strcmp(CG.type[1], "residue"  ) == 0) {
        
        sprintf(buffer, "%s/Library_Go/CGParameterTable.txt", datadir);
        input_file = fopen(buffer, "r");
        if (input_file == NULL) {
            printf("!!ERROR!!: cannot find CGParameterTable.txt file in %s/Library_Go/!\n", datadir);
            printf("           please check your directory!\n");
            exit(EXIT_FAILURE);
        }
        
        while (fgets(statement, sizeof(statement), input_file)) {
            if (statement[0] == '#') {
                continue;
            }
            
            int pos = 0;
            sscanf(statement, "%s%n", buffer, &pos);
            if (strcmp(buffer, "THETA") == 0) {
                sscanf(statement + pos, "%lf%lf%lf%lf%lf",
                       &CGMatrix.theta[0],
                       &CGMatrix.theta[1],
                       &CGMatrix.theta[2],
                       &CGMatrix.theta[3],
                       &CGMatrix.theta[4]);
            } else {
                sscanf(statement + pos, "%lf%lf%lf%lf",
                       &CGMatrix.CGBead[AAModel(buffer)].chi,
                       &CGMatrix.CGBead[AAModel(buffer)].epsilon,
                       &CGMatrix.CGBead[AAModel(buffer)].sigma,
                       &CGMatrix.CGBead[AAModel(buffer)].delta);
            }
        }
        fclose(input_file);
        
        CreateCGPotentialMatrix();
        ReadCGPotTable();
        
#ifdef CG_DEBUG
        PrintStepPotFun();
#endif
    }
    
    return;
}


void ReadGoPara(void) {
	char buffer[512], statement[128];
	double tmp, **sigmaData;
	FILE *input_file;

	//parameter file generated by
	//http://www.mmtsb.org/webservices/gomodel.html
	sprintf(buffer, "%s/GoPara.txt", datadir);
	input_file = fopen(buffer, "r");

	DOUBLE_2CALLOC(sigmaData, NATOMTYPE, 2);
	for (int i = 1; i <= totalAminoAcids; i++) {
		fgets(statement, sizeof(statement), input_file);
		sscanf(statement, "%s%s%s%lf", buffer, buffer, buffer, &tmp);
		sigmaData[AAModel(aminoacid[i].nameofAA)][0] += tmp;
		sigmaData[AAModel(aminoacid[i].nameofAA)][1] ++;
	}

	for (int i = 0; i < NATOMTYPE; i++) {
		if (sigmaData[i][1] != 0) {
			CGMatrix.CGBead[i].sigma = sigmaData[i][0] / sigmaData[i][1];
		}
		else {
			if (CGMatrix.CGBead[i].sigma < 0) {
				CGMatrix.CGBead[i].sigma = -1 * INFINIT;
			}
		}
	}
	fclose(input_file);
	FREE2(sigmaData);

	return;
}


void CreateCGPotentialMatrix(void) {
    char buffer[512];
    char *sur_type[3] = {"PHO", "PHI", "RPH"}; //these are corresponding to the surface type, see AAModel() in Model.c.
    int step, num_sur = 3, flag;
    double minPDis, minDis;
    double gap[3] = {1.0, 2.0, 0.01}; //0: potential gap; 1: distance gap; 2: step length
    double dis, disPot[128][2] = {0};
    FILE *outputPoten;
    
#ifdef CG_DEBUG
    PrintPotCur_Res();
    PrintPotCur_Pep();
#endif
    
    sprintf(buffer, "%s/Library_Go/CGPotenTable.txt", datadir);
    outputPoten = fopen(buffer, "w");
    
    for (int i = 0; i < num_sur; i ++) {
        for (int n = 0; n < NATOMTYPE; n ++) {
            
			//the initial sigma = 0, this is to void the unavailable types
            if (CGMatrix.CGBead[n].sigma <= 0 ||
                AAModel(sur_type[0]) == n ||
                AAModel(sur_type[1]) == n ||
                AAModel(sur_type[2]) == n) {
                continue;
            }
            
            //the minimum distance
            minDis  = FindMaxPotDis(CGMatrix.CGBead[AAModel(sur_type[i])],
                                    CGMatrix.CGBead[n], MAX_P);
            //the distance at which the potenial has the minimum
            minPDis = FindMinPotDis(CGMatrix.CGBead[AAModel(sur_type[i])],
                                    CGMatrix.CGBead[n], 0);
            
            fprintf(outputPoten, "%5s %5s ", sur_type[i], GetAAName(n));
            fprintf(outputPoten, "%10.6lf ", minDis);
            
            step = 0; flag = 1;
            dis = minDis + gap[2];
            disPot[step][0] = dis;
            disPot[step++][1] = CalculateSurPoten(CGMatrix.CGBead[AAModel(sur_type[i])],
                                                  CGMatrix.CGBead[n], dis);
            
            //calculate the discretized points of potential
            while (dis <= MAX_R) {
                dis += gap[2];
                
                if (step == 128) {
                    printf("!!ERROR!!: too many steps! please increase the step gap. %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                
                disPot[step][0] = dis;
                disPot[step][1] = CalculateSurPoten(CGMatrix.CGBead[AAModel(sur_type[i])],
                                                    CGMatrix.CGBead[n], dis);
                
                if (flag && fabs(disPot[step][0] - minPDis) < gap[2]) { //only do this once
                    disPot[step][0] = minPDis;
                    disPot[step][1] = CalculateSurPoten(CGMatrix.CGBead[AAModel(sur_type[i])],
                                                        CGMatrix.CGBead[n], minPDis);
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
            
#ifdef CG_DEBUG
            PrintDisPot(disPot, step, sur_type[i], n);
#endif
            
            int stepRecord = step; step = 0;
            while (++step < stepRecord) {
                fprintf(outputPoten, "%10.6lf %10.6lf ", (disPot[step - 1][0] + disPot[step][0]) * 0.5, (disPot[step - 1][1] - disPot[step][1]));
            }
            fprintf(outputPoten, "%10.6lf %10.6lf\n", disPot[step - 1][0], disPot[step - 1][1]);
        }
    }
    
    fclose(outputPoten);
    return;
}


void PrintPotCur_Res(void) {
    FILE *potCurFile;
    char name[128];
    char *sur_type[3] = {"PHO", "PHI", "RPH"}; //these are corresponding to the surface type, see AAModel() in Model.c.
    int num_sur = 3;
    double dis, gap = 0.02;
    double poten = 0;

    for (int i = 0; i < num_sur; i ++) {
        for (int n = 0; n < NATOMTYPE; n ++) {
            
            if (CGMatrix.CGBead[n].sigma <= 0 ||
                AAModel(sur_type[0]) == n ||
                AAModel(sur_type[1]) == n ||
                AAModel(sur_type[2]) == n) {
                continue;
            }
            
            sprintf(name, "%s/%s-%s.txt", savedir, sur_type[i], GetAAName(n));
            potCurFile = fopen(name, "w");
            
            dis = FindMaxPotDis(CGMatrix.CGBead[AAModel(sur_type[i])],
                                CGMatrix.CGBead[n], 0.6);
            while (dis < MAX_R) {
                poten = CalculateSurPoten(CGMatrix.CGBead[AAModel(sur_type[i])],
                                          CGMatrix.CGBead[n], dis);
                fprintf(potCurFile, "%8.4lf  %8.4lf\n", dis, poten);
                
                dis += gap;
            }
            fprintf(potCurFile, "\n");
            fclose(potCurFile);
        }
    }
    return;
}


void PrintPotCur_Pep(void) {
    FILE *potCurFile;
    char name[128];
    char *sur_type[3] = {"PHO", "PHI", "RPH"}; //these are corresponding to the surface type, see AAModel() in Model.c.
    int num_sur = 3, dim = 2;
    double dis, shift, gap = 0.1;
    double poten = 0;
    
    for (int i = 0; i < num_sur; i ++) {
        sprintf(name, "%s/%s-Pep.txt", savedir, sur_type[i]);
        potCurFile = fopen(name, "w");
        
        dis   = 2.0;
        shift = FindMinDisDimen(dim); //on y-axis

        while (dis < MAX_R) {
            poten = 0;
            for (int n = 1; n <= totalAminoAcids; n ++) {
                
                int AAHead = aminoacid[n].startAtomNum;
                int AAEnd  = aminoacid[n].endAtomNum;
                
                for (int j = AAHead; j <= AAEnd; j ++) {
                    if (strcmp(atom[j].property->extraProperty[0], "CG") == 0) {
                        poten += CalculateSurPoten(CGMatrix.CGBead[AAModel(sur_type[i])],
                                                   CGMatrix.CGBead[aminoacid[n].type],
                                                   dis + atom[j].dynamic->coordinate[dim] - shift);
                    }
                }
            }
            
            if (poten < MAX_P) fprintf(potCurFile, "%8.4lf  %8.4lf\n", dis, poten);
            dis += gap;
        }
        fprintf(potCurFile, "\n");
        fclose(potCurFile);
    }
    return;
}


double FindMinDisDimen(int dim) {
    double dis = INFINIT;
    
    for (int n = 1; n <= totalAminoAcids; n ++) {
        int AAHead = aminoacid[n].startAtomNum;
        int AAEnd  = aminoacid[n].endAtomNum;
        
        for (int i = AAHead; i <= AAEnd; i ++) {
            if (strcmp(atom[i].property->extraProperty[0], "CG") == 0) {
                if (dis > atom[i].dynamic->coordinate[dim]) {
                    dis = atom[i].dynamic->coordinate[dim];
                }
            }
        }
    }
    
    return dis;
}


void PrintDisPot(double pot[128][2], int step, char *sur_type, int num) {
    FILE *disPotCurFile;
    char name[128];
    
    sprintf(name, "%s/disPot_%s-%s.txt", savedir, sur_type, GetAAName(num));
    disPotCurFile = fopen(name, "w");
    
    for (int j = 0; j < step; j ++) {
        fprintf(disPotCurFile, "%8.4lf %8.4lf\n", pot[j][0], pot[j][1]);
    }
    fclose(disPotCurFile);
    return;
}


void PrintStepPotFun(void) {
    FILE *stepFile;
    char name[128];
    char *sur_type[3] = {"PHO", "PHI", "RPH"}; //these are corresponding to the surface type, see AAModel() in Model.c.
    int num_sur = 3;
    
    for (int i = 0; i < num_sur; i ++) {
        for (int n = 0; n < NATOMTYPE; n ++) {
            
            if (CGMatrix.CGBead[n].sigma <= 0 ||
                AAModel(sur_type[0]) == n ||
                AAModel(sur_type[1]) == n ||
                AAModel(sur_type[2]) == n) {
                continue;
            }
            
            sprintf(name, "%s/stepFun_%s-%s.txt", savedir, sur_type[i], GetAAName(n));
            stepFile = fopen(name, "w");
            
            struct StepPotenStr *thisStep = potentialPairCG[AAModel(sur_type[i])][n].step;
            while (thisStep != NULL) {
                fprintf(stepFile, "%10.4lf  %.4lf\n", sqrt(thisStep->d), thisStep->accumulated);
                thisStep = thisStep->next;
            }
            fprintf(stepFile, "\n");
            fclose(stepFile);
        }
    }
    return;
}


double CalculateSurPoten(struct CGBeadStr sur, struct CGBeadStr res, double dis) {
    double sigma, epsilon, ratio;
    
    epsilon = sur.epsilon * res.epsilon;
    if (epsilon < 0) {
        epsilon = sqrt(epsilon * -1);
    } else {
        epsilon = sqrt(epsilon);
    }
    
    sigma = (sur.sigma + res.sigma) * 0.5;
    ratio = sigma / dis;
    
    return M_PI * POT_RO * CalExponential(sigma, 3.0) * epsilon *
           (CGMatrix.theta[0] * CalExponential(ratio, 9) -
            CGMatrix.theta[1] * CalExponential(ratio, 7) +
            CGMatrix.theta[2] * CalExponential(ratio, 3) -
           (CGMatrix.theta[3] *  sur.chi +
            CGMatrix.theta[4] * (res.chi + 2 * res.delta)) *
            CalExponential(ratio, 3));
}


double FindMaxPotDis(struct CGBeadStr sur, struct CGBeadStr res, double setP) {
	double dis = 0.1;
	
	while (CalculateSurPoten(sur, res, dis) > setP) {
		dis += 0.02;
	}

	return dis;
}


double FindMinPotDis(struct CGBeadStr sur, struct CGBeadStr res, int print) {
    double min = INFINIT, minDis = 0;
    double dis = 1.0, value;
    
    while (dis < MAX_R) {
		value = CalculateSurPoten(sur, res, dis);
        
        if (print) {
            printf("%12.4lf%12.4lf\n", dis, value);
        }
        
        if (value < min) {
            min    = value;
            minDis = dis;
        }
        dis += 0.001;
    }
    return minDis;
}


double FindMaxRadius(struct CGBeadStr sur, struct CGBeadStr res, double minPDis, double diff) {
	double pot[2];
	double dis = minPDis;

	pot[0] = CalculateSurPoten(sur, res, dis);
	pot[1] = CalculateSurPoten(sur, res, dis + 0.5);
	while (pot[0] < 0 || pot[1] < 0 || fabs(pot[0] - pot[1]) > diff) {
		dis += 0.5;
		pot[0] = CalculateSurPoten(sur, res, dis);
		pot[1] = CalculateSurPoten(sur, res, dis + 0.5);
        
        if (dis > MAX_R) return MAX_R;
	}
	
	return dis;
}


void ReadCGPotTable(void) {
    //create a matrix storing the distances and energy data of inter-molecular interactions
    int num[2] = {0};
    char type[2][10];
    char buffer[1024];
    double value[2];
    struct StepPotenStr *thisStep = NULL;
    FILE *input_file;
    
    sprintf(buffer, "%s/Library_Go/CGPotenTable.txt", datadir);
    input_file = fopen(buffer, "r");
    if (input_file == NULL) {
        printf("!!ERROR!!: cannot find CGPotenTable.txt file in %s/Library_Go/!\n", datadir);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }
    
    while (fgets(buffer, sizeof(buffer), input_file)) {
        if (buffer[0] == '#') {
            continue;
        }
        
        int pos, curpos = 0; //record and move the pointer on buffer
        sscanf(buffer, "%s%s%n", type[0], type[1], &pos);
        curpos += pos;
        
        num[0] = AAModel(type[0]);
        num[1] = AAModel(type[1]);
        sscanf(buffer + curpos, "%lf%n", &potentialPairCG[num[0]][num[1]].dmin, &pos);
        curpos += pos;
        
        potentialPairCG[num[0]][num[1]].dmin *= potentialPairCG[num[0]][num[1]].dmin; //square all distance variables, convenient to distance comparison! Exclude bond length.
        
        int tmp = 1;
        while (sscanf(buffer + curpos, "%lf%lf%n", &value[0], &value[1], &pos) != EOF) {
            curpos += pos;
            
            struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
            newStep->d = value[0];
            newStep->d *= newStep->d;
            newStep->e = value[1];
            
            if (tmp == 1) {
                potentialPairCG[num[0]][num[1]].step = newStep;
                thisStep = newStep;
                tmp ++;
            } else {
                thisStep->next = newStep;
                thisStep = thisStep->next;
            }
        }
        
        /*
         calculate the potential depth/height of the interactions
         the input potential is just the potential change at the potential step,
         we need to accumulate it to get the actual potential.
         also due to the list is from small to large,
         we have to inverse the accumulated potential.
         see "CalAccumulatedPoten" in ToolFunctions.c
         */
        CalAccumulatedPoten(potentialPairCG[num[0]][num[1]].step);
        potentialPairCG[num[1]][num[0]] = potentialPairCG[num[0]][num[1]];
    }
    fclose(input_file);
    
    return;
}


void CGSurfaceTime(struct AtomStr *target) {
    int atomNum = target->property->num;
    struct ParameterStr parameters;
    struct InteractionEventStr event;
    
    for (int n = 0; n < obstObj.num; n ++) {
        InitializeEvent(&event);
        
        struct AtomStr *thisSurface = &obstObj.obst[n];
        for (int dim = 1; dim <= 3; dim ++) {
            if (obstObj.position[n][dim] < 0) {
                continue;
            }
            
            double distance2 = 0;
            TRANSFER_VECTOR(parameters.position_i, target->dynamic->coordinate);
            TRANSFER_VECTOR(parameters.speed_i, target->dynamic->velocity);
            
            TRANSFER_VECTOR(parameters.position_j, thisSurface->dynamic->coordinate);
            TRANSFER_VECTOR(parameters.speed_j, thisSurface->dynamic->velocity);
            
            for (int i = 1; i <= 3; i ++) {
                if (i != dim) {
                    parameters.position_i[i] = 0;
                    parameters.position_j[i] = 0;
                    
                    parameters.speed_i[i] = 0;
                    parameters.speed_j[i] = 0;
                }
            }
            
            DOT_MINUS(parameters.position_i, parameters.position_j, parameters.r_ij);
            DOT_MINUS(parameters.speed_i,    parameters.speed_j,    parameters.v_ij);
            parameters.r_2  = DOT_PROD(parameters.r_ij, parameters.r_ij);
            parameters.b_ij = DOT_PROD(parameters.r_ij, parameters.v_ij);
            parameters.v_2  = DOT_PROD(parameters.v_ij, parameters.v_ij);
            distance2 = parameters.r_2;
            
            if (strcmp(target->property->extraProperty[1], "CONS") == 0) {
                if (distance2 < preBond.dis) {
                    printf("!!ERROR!!: the initial distance between the fixing atom %i and the surface is too short!\n", target->property->num);
                    printf("           the current is %.4lf, while the assigned is %.4lf\n", sqrt(distance2), sqrt(preBond.dis));
                    printf("           %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                parameters.shortlimit2 = preBond.dis;
                parameters.longlimit2  = INFINIT;
                parameters.lowerPotential = 0;
                parameters.upperPotential = INFINIT;
                parameters.coreorShell = 0;
            } else {
                parameters.coreorShell = FindPair(target, thisSurface, "CG",
                                                  parameters.b_ij, distance2,
                                                  &parameters.shortlimit2, &parameters.longlimit2,
                                                  &parameters.lowerPotential, &parameters.upperPotential,
                                                  &tmpDouble, NULL, 0);
            }
            
            if (parameters.b_ij < 0) { //approaching
                if (parameters.coreorShell == 0) {
                    event.subEventType = CoreColli;
                } else {
                    event.subEventType = TBDCaptr;
                }
                
                parameters.s_time = -1;
                parameters.d_ij2  = parameters.shortlimit2;
                event.potential   = parameters.lowerPotential;
                
            } else { //leaving
                if (parameters.coreorShell == 1) {
                    continue;
                } else {
                    event.subEventType = TBDEscap;
                }
                
                parameters.s_time = 1;
                parameters.d_ij2  = parameters.longlimit2;
                event.potential   = parameters.upperPotential;
            }
            
            event.time = CalculateTime(&parameters);
            if (unlikely(event.time < 0)) {
                printf("!!ERROR!!: CG time is less than zero, atom #%i! %s:%i\n", atomNum, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            if (JobAssign(target, NULL, &event, CGSu_Event)) {
                target->dynamic->event.partner = -1 * dim - 10 * n;
            }
        }
    }
    
    return;
}


void CGSurfaceEvent(struct AtomStr *target, struct AtomListStr *thisList) {
    int s_vel = 0;
    double v_ij[4];
    double b_ij, d_ij2, r_ij[4], term[4];
    double phi, potential, reducedMass = target->property->mass;
    double speed_i[4] = {0}, speed_j[4] = {0};
    double position_i[4] = {0}, position_j[4] = {0};
    struct AtomStr *thisSurface = &obstObj.obst[target->dynamic->event.partner / -10];
    
    if (strcmp(target->property->extraProperty[1], "CONS") == 0) {
        sprintf(target->property->extraProperty[1], "_CONS");
        target->property->mass = INFINIT * INFINIT;
        target->dynamic->event.counter ++;
        target->dynamic->event.time = INFINIT;
        target->dynamic->velocity[1] = 0;
        target->dynamic->velocity[2] = 0;
        target->dynamic->velocity[3] = 0;
        
        return;
    }
    
    int dim = (target->dynamic->event.partner * -1) % 10;
    
    TRANSFER_VECTOR(position_i,  target->dynamic->coordinate);
    TRANSFER_VECTOR(speed_i,     target->dynamic->velocity);
    
    TRANSFER_VECTOR(position_j, thisSurface->dynamic->coordinate);
    TRANSFER_VECTOR(speed_j,    thisSurface->dynamic->velocity);
    
    for (int i = 1; i <= 3; i ++) {
        if (dim != i) {
            position_i[i] = 0; speed_i[i] = 0;
            position_j[i] = 0; speed_j[i] = 0;
        }
    }
    
    DOT_MINUS(position_i, position_j, r_ij);
    DOT_MINUS(   speed_i,    speed_j, v_ij);
    d_ij2     = DOT_PROD(r_ij, r_ij);
    b_ij      = DOT_PROD(r_ij, v_ij);
    
    if (strcmp(CG.type[1], "residue") == 0) {
        int AANum = protein[target->property->sequence.proteinNum - 1].endAANum + target->property->sequence.aminoacidNum;
        potential = target->dynamic->event.potential / aminoacid[AANum].mass * target->property->mass;
    } else if (strcmp(CG.type[1], "alphaC") == 0) {
        potential = target->dynamic->event.potential;
	} else {
		printf("!!ERROR!!: the CG type is invalid! %s:%i\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

    switch (target->dynamic->event.subEventType) {
        case CoreColli:
            s_vel = 1;
            potential = 0;
            break;
            
        case TBDCaptr:
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to capture?
                // >0 -> sufficient; <=0 -> no
                potential = 0;
                s_vel = 1;
            } else {
                s_vel = -1;
            }
            break;
            
        case TBDEscap: //TBD outer
            if (b_ij * b_ij - 2 * d_ij2 * potential / reducedMass <= 0) { //sufficient energy to escape?
                potential = 0;
                s_vel = -1;
            } else {
                s_vel = 1;
            }
            break;
            
        default:
            printf("!ERROR!: wall interaction type error! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            break;
    }
    
    phi = (-1 * b_ij + s_vel * sqrt(b_ij * b_ij - 2 * d_ij2 * potential / reducedMass)) / d_ij2;
    
    if (unlikely(isnan(phi))) {
        printf("!!ERROR!!: phi is NaN! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    FACTOR_PROD(phi, r_ij, term);
    DOT_PLUS(target->dynamic->velocity, term, target->dynamic->velocity);
    
    return;
}


void CGResidueRenew(int *renewList, struct AtomListStr *thisList) {
    int CGNum = 0, startNum = 0, endNum = 0;
    struct AtomStr *target = NULL, *alphaC = thisList->ptr[renewList[1]];
    
	for (int n = 1, i = renewList[0] + 1; n <= renewList[0]; n ++) {
		target = thisList->ptr[renewList[n]];
        if (strcmp(target->property->extraProperty[0], "CG")) { // non-CG atom
            continue;
        }
        
		CGNum = protein[target->property->sequence.proteinNum - 1].endAANum + target->property->sequence.aminoacidNum;
		startNum = CG.thisCG[CGNum].start;
		endNum = CG.thisCG[CGNum].end;
        
		for (int j = startNum; j <= endNum; j ++) {
			if (!FindElemInList(renewList, j, 1, renewList[0], 0)) renewList[i++] = j;
		}
		renewList[i] = 0;
    }

	for (int i = renewList[0] + 1; renewList[i]; i ++) {
		target = thisList->ptr[renewList[i]];
		target->dynamic->event.potential    = alphaC->dynamic->event.potential;
		target->dynamic->event.subEventType = alphaC->dynamic->event.subEventType;
		target->dynamic->event.eventType    = alphaC->dynamic->event.eventType;
		target->dynamic->event.partner      = alphaC->dynamic->event.partner;

		CGSurfaceEvent(target, thisList);
		target->dynamic->event.counter++;
	}
    
    return;
}


void FixAssignRenewList(int *renewList, struct AtomListStr *thisList) {
    int cell_neighbor[4], cellIndex;
    int *sourCellList = celllist;
    int partner;
    int pos = renewList[0] + 1;
    struct AtomStr *target = thisList->ptr[renewList[1]];
    
    for (int i = 0; i < 27; ++i) {
        //scan the neighborhood 27 subcells, include the target subcell itself
        cell_neighbor[3] = target->dynamic->cellIndex[3] + thisList->z[i];
        cell_neighbor[2] = target->dynamic->cellIndex[2] + thisList->y[i];
        cell_neighbor[1] = target->dynamic->cellIndex[1] + thisList->x[i];
        
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
        
        partner = sourCellList[atomnum + cellIndex];
        while (partner) {
            if (thisList->ptr[partner]->dynamic->event.partner == target->property->num &&
                !FindElemInList(renewList, partner, 1, renewList[0], 0)) renewList[pos++] = partner;
            partner = sourCellList[partner];
        }
    }
    
    return;
}


char *GetAAName(int AAType) {
    switch (AAType) {
        case 1:
            return "ALA";
        case 2:
            return "ARG";
        case 3:
            return "ASN";
        case 4:
            return "ASP";
        case 5:
            return "CYS";
        case 6:
            return "GLN";
        case 7:
            return "GLU";
        case 8:
            return "GLY";
        case 9:
            return "HIS";
        case 10:
            return "ILE";
        case 11:
            return "LEU";
        case 12:
            return "LYS";
        case 13:
            return "MET";
        case 14:
            return "PHE";
        case 15:
            return "PRO";
        case 16:
            return "SER";
        case 17:
            return "THR";
        case 18:
            return "TRP";
        case 19:
            return "TYR";
        case 20:
            return "VAL";
        case 21:
            return "SOL";
        case 22:
            return "GEL";
        case 23:
            return "PHO";
        case 24:
            return "PHI";
        case 25:
            return "RPH";
        default:
            printf("\n!ERROR!:The amino acid type cannot be found!\n");
            exit(EXIT_FAILURE);
    }
    return "Invalid";
}

