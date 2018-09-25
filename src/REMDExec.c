//
//  main.c
//  REMD
//
//  Created by Size Zheng on 11/12/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/wait.h>

#define boltzmann 0.0019872041 //kcal / mol / K

int main (int argc, char *argv[]) {
    pid_t child_pid, wpid;
    int next = 0;
    int status = 0;
    int replicaNum = 0;
    int distStatus = 1;
    char portNo[20], exRate[20], inputT[20], exName[20];
    char buffer[1024], execClient[1024], execServer[1024];
    char directory[1024];
    char *command[50];
    double *T, thisT;
    FILE *configFile;
    
    if (getcwd(directory, sizeof(directory)) == NULL)
        perror("!!ERROR!!");
    
    for (int n = 0; n < 50; n ++) {
        command[n] = (char *)calloc(50, sizeof(char));
    }
    
    for (int n = 1; n < argc; n ++) {
        if (argv[n][0] == '-') {
            if (strcmp(argv[n], "-h") == 0 || strcmp(argv[n], "-help") == 0) {
            help:
                printf("Usage: %s [distribute T or not] -f [configuration file] -args [args of executable without flag -REMD]\n\n", argv[0]);
                printf("       -nodist: [yes] default is yes. will distribute temperature to each replica.\n");
                printf("            -f: provide the configuration file.\n");
                printf("         -args: provide the flags and info for sDMD without -REMD, end with \"-end\"\n");
                exit(EXIT_SUCCESS);
            } else if (strcmp(argv[n], "-nodist") == 0) {
                distStatus = 0;
            } else if (strcmp(argv[n], "-f") == 0) {
                sprintf(buffer, "%s", argv[n + 1]);
            } else if (strcmp(argv[n], "-args") == 0) {
                next = 1;
                while (n + next < argc && strcmp(argv[n + next], "-end")) {
                    sprintf(command[next - 1], "%s", argv[n + next]);
                    next ++;
                }
                n += next;
            }
        }
    }
    next --;
    
    sprintf(execClient, "%s/sDMD", directory);
    sprintf(execServer, "%s/sServer", directory);

    configFile = fopen(buffer, "r");
    if (configFile == NULL) {
        printf("!!ERROR!!: cannot find the configuration file at %s\n", buffer);
        exit(EXIT_FAILURE);
    }
    fscanf(configFile, "%[^\n]\n", buffer);
    fscanf(configFile, "%s%s%i\n", buffer, buffer, &replicaNum);
    
    int pos = 0, cur = 0;
    T = (double *)calloc(replicaNum, sizeof(double));
    fscanf(configFile, "%[^\n]\n", buffer);
    sscanf(buffer, "%s%s%n", directory, directory, &cur);
    pos += cur;
    
    for (int i = 0; i < replicaNum; i ++) {
        if (sscanf(buffer + pos, "%lf%n", &T[i], &cur) == EOF) {
            printf("!!ERROR!!: the number of temperatures provided is less than the number of replica!\n");
            exit(EXIT_FAILURE);
        }
        T[i] *= boltzmann;
        pos += cur;
    }
    if (sscanf(buffer + pos, "%lf", &thisT) != EOF) {
        printf("!!ERROR!!: the number of temperatures provided is larger than the number of replica!\n");
        exit(EXIT_FAILURE);
    }
    
    fscanf(configFile, "%s%s%s", buffer, buffer, portNo);
    fscanf(configFile, "%s%s%s", buffer, buffer, exRate);
    fclose(configFile);
    
    for (int id = 0; id < replicaNum; id++) {
        if (distStatus) {
            thisT = T[id];
        } else {
            thisT = 0;
        }
        
        sprintf(inputT, "%lf", thisT);
        sprintf(exName, "%i", id);
        
        if (next < 0) goto help;
        sprintf(command[next],     "-REMD");
        sprintf(command[next + 1], "localhost");
        sprintf(command[next + 2], "%s", portNo);
        sprintf(command[next + 3], "%s", inputT);
        sprintf(command[next + 4], "%s", exRate);
        sprintf(command[next + 5], "%s", exName);
        command[next + 6] = (char *)0;
        
        if ((child_pid = fork()) == 0) {
            
            //direct the stdout of execv() to file outputData
            char outputData[50];
            sprintf(outputData, "outputData_%s.txt", exName);
            int fd = open(outputData, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
            dup2(fd, 1);
            close(fd);
            
            sleep(1 + id);
            if (execv(execClient, command) < 0)
                perror("!!ERROR!!");
            exit(0);
        }
    }
    
    char repNum[20]; sprintf(repNum, "%i", replicaNum);
    char seed[20];   sprintf(seed, "%i", (unsigned)time(NULL));
    char *serverCommand[] = {"./sServer", portNo, repNum, seed, (char *)0};
    if (execv(execServer, serverCommand) < 0)
        perror("!!ERROR!!");
    
    while ((wpid = wait(&status)) > 0);
    
    return 0;
}


