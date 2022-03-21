//
//  FF.c
//  sDMD
//
//  Created by Size on 2020/7/19.
//  Copyright © 2020 Size Zheng. All rights reserved.
//

#include "DMD.h"

#define MAX_R       10
#define MAX_P       0.8

static void CreateBondData(char *FF, char *data, int pos, FILE *outputFile);
static void CreateAngleData(char *FF, char *data, int pos, FILE *outputFile);
static void CreateLJData(char *FF, char *data, int pos, FILE *outputFile);
static double FindMaxPotDis(char *type, double *parameter, double setP);
static double FindMinPotDis(char *type, double *parameter, double setP);
static double CalculateBond(char *type, double kb, double b0, double dis);
static double CalculateAngle(char *type, double ktheta, double theta0, double angle);
static double CalculateLJ(char *type, double sigma, double radius, double dis);

void CreateFFMatrix(char *FF) {
    char buffer[512], title[16];
    int pos = 0;
    FILE *inPara, *outputBond, *outputAngle, *outputLJ;
    
    sprintf(buffer, "%s/Library_%s/FFParameter/SOL.txt", datadir, FF);
    inPara = fopen(buffer, "r");

    sprintf(buffer, "%s/Library_%s/AA/SOL.txt", datadir, FF);
    outputBond = fopen(buffer, "w");
    outputAngle = outputBond;

    sprintf(buffer, "%s/Library_%s/InteractionPotentialTable.txt", datadir, FF);
    outputLJ = fopen(buffer, "w");
    
    while (fgets(buffer, sizeof(buffer), inPara)) {
        sscanf(buffer, "%s%n", title, &pos);
        if (strcmp(title, "BOND") == 0) {
            CreateBondData(FF, buffer, pos, outputBond);
        } else if (strcmp(title, "ANGLE") == 0) {
            CreateAngleData(FF, buffer, pos, outputAngle);
        } else if (strcmp(title, "LJ") == 0) {
            CreateLJData(FF, buffer, pos, outputLJ);
        } else if (buffer[0] == '#' ||
                   strcmp(title, "ATOM") == 0) {
            fprintf(outputBond, "%s", buffer);
        } else {
            printf("!!ERROR!!: the potential function %s is not supported!\n", title);
            printf("           %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
    
    fclose(inPara);
    fclose(outputBond);
    fclose(outputLJ);
    return;
}


void CreateBondData(char *FF, char *data, int pos, FILE *outputFile) {
    char atomType[2][32];
    double parameter[2], maxDis;
    sscanf(data + pos, "%s%s%lf%lf", atomType[0], atomType[1], &parameter[0], &parameter[1]);
    
    maxDis = parameter[1];
    while (CalculateBond(FF, parameter[0], parameter[1], maxDis) < MAX_P) {
        maxDis += 0.02;
    }
    
    fprintf(outputFile, "   BOND %-7s%-7s", atomType[0], atomType[1]);
    fprintf(outputFile, "%-10.6lf %-10.6lf\n", parameter[1], maxDis / parameter[1] - 1);
    
    return;
}


void CreateAngleData(char *FF, char *data, int pos, FILE *outputFile) {
    char atomType[2][32];
    double parameter[4], maxDeg;
    sscanf(data + pos, "%s%s%lf%lf%lf%lf", atomType[0], atomType[1],
           &parameter[0],
           &parameter[1],
           &parameter[2],
           &parameter[3]);
    
    maxDeg = parameter[1];
    while (CalculateAngle(FF, parameter[0], parameter[1], maxDeg) < MAX_P) {
        maxDeg += 0.02;
    }
    
    double originalDis = sqrt(parameter[2] * parameter[2] + parameter[3] * parameter[3] -
                              2 * parameter[2] * parameter[3] * cos(parameter[1] / 180 * M_PI));
    double extendedDis = sqrt(parameter[2] * parameter[2] + parameter[3] * parameter[3] -
                              2 * parameter[2] * parameter[3] * cos(maxDeg / 180 * M_PI));
    
    fprintf(outputFile, "  ANGLE %-7s%-7s", atomType[0], atomType[1]);
    fprintf(outputFile, "%-10.6lf %-10.6lf\n", originalDis, extendedDis / originalDis - 1);
    
    return;
}


void CreateLJData(char *FF, char *data, int pos, FILE *outputFile) {
    char atomType[2][32];
    int step, flag;
    double gap[3] = {1.0, 2.0, 0.01}; //0: potential gap; 1: distance gap; 2: step length
    double parameter[2], minDis, minPDis;
    double dis, disPot[128][2] = {0};
    
    sscanf(data + pos, "%s%s%lf%lf", atomType[0], atomType[1], &parameter[0], &parameter[1]);

    minDis  = FindMaxPotDis(FF, parameter, MAX_P);
    minPDis = parameter[1];
    
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

    fprintf(outputFile, "%5s %5s ", atomType[0], atomType[1]);
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
    
    return;
}


double FindMaxPotDis(char *type, double *parameter, double setP) {
	double dis = parameter[1] * 0.5;
	
    while (CalculateLJ(type, parameter[0], parameter[1], dis) > setP) {
        dis += 0.02;
    }

	return dis;
}


double FindMinPotDis(char *type, double *parameter, double setP) {
    double min = INFINIT, minDis = 0;
    double dis = 1.0, value;
    
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


double CalculateBond(char *type, double kb, double b0, double dis) {
    if (strcmp(type, "CHARMM19") == 0) {
        return kb * CalExponential((dis - b0), 2);
    }
    
    printf("!!ERROR!!: FF type for bond calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


double CalculateAngle(char *type, double ktheta, double theta0, double angle) {
    if (strcmp(type, "CHARMM19") == 0) {
        return ktheta * CalExponential((angle - theta0), 2);
    }
    
    printf("!!ERROR!!: FF type for angle calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}


double CalculateLJ(char *type, double sigma, double radius, double dis) {
    if (strcmp(type, "CHARMM19") == 0) {
        double rmin  = radius;
        double ratio = rmin / dis;
        return -1 * sigma * (CalExponential(ratio, 12) - 2 * CalExponential(ratio, 6));
    }
    
    printf("!!ERROR!!: FF type for LJ calculation is not available!\n");
    printf("           %s:%i\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}
