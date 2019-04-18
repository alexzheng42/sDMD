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

#define CHAR_LENGTH 1024
#define BOLTZMANN 0.0019872041 //kcal / mol / K

#define CHAR_LOC(name) \
        name = (char *)calloc(CHAR_LENGTH, sizeof(char));

struct REMDTempStr {
    int num;
    double T;
};

int main (int argc, char *argv[]) {
    pid_t child_pid, wpid;
    int next = 0;
    int status = 0;
    int replicaNum = 0;
    int distStatus = 1;
    char portNo[32], exRate[32];
    char *buffer, *execClient, *execServer;
    char *directory;
    char *command[32];
    struct REMDTempStr thisT;
    struct REMDTempStr *REMD_T;
    FILE *configFile;
    
    CHAR_LOC(buffer);
    CHAR_LOC(execClient);
    CHAR_LOC(execServer);
    CHAR_LOC(directory);
    
    if (getcwd(directory, CHAR_LENGTH) == NULL)
        perror("!!ERROR!!");
    
    for (int n = 0; n < 32; n ++) {
        command[n] = (char *)calloc(64, sizeof(char));
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
    fgets(buffer, CHAR_LENGTH, configFile);
    fscanf(configFile, "%s%s%i\n", buffer, buffer, &replicaNum);
    
    int pos = 0, cur = 0;
    REMD_T = (struct REMDTempStr *)calloc(replicaNum, sizeof(struct REMDTempStr));
    fgets(buffer, CHAR_LENGTH, configFile);
    sscanf(buffer, "%s%s%n", directory, directory, &cur);
    pos += cur;
    
    for (int i = 0; i < replicaNum; i ++) {
        if (sscanf(buffer + pos, "%lf%n", &REMD_T[i].T, &cur) == EOF) {
            printf("!!ERROR!!: the number of temperatures provided is less than the number of replica!\n");
            exit(EXIT_FAILURE);
        }
        REMD_T[i].T *= BOLTZMANN;
        REMD_T[i].num = i + 1;
        pos += cur;
    }
    if (sscanf(buffer + pos, "%lf", &(thisT.T)) != EOF) {
        printf("!!ERROR!!: the number of temperatures provided is larger than the number of replica!\n");
        exit(EXIT_FAILURE);
    }
    
    fscanf(configFile, "%s%s%s", buffer, buffer, portNo);
    fscanf(configFile, "%s%s%s", buffer, buffer, exRate);
    fclose(configFile);
    
    for (int id = 0; id < replicaNum; id++) {
        if (distStatus) {
            thisT = REMD_T[id];
        } else {
            thisT.T = 0;
            thisT.num = 0;
        }
        
        if (next < 0) goto help;
        sprintf(command[next],     "-REMD");
        sprintf(command[next + 1], "localhost");
        sprintf(command[next + 2], "%s" , portNo);
        sprintf(command[next + 3], "%lf", thisT.T);
        sprintf(command[next + 4], "%i" , thisT.num);
        sprintf(command[next + 5], "%s" , exRate);
        sprintf(command[next + 6], "%i" , id);
        command[next + 7] = (char *)0;
        
        if ((child_pid = fork()) == 0) {
            
            //direct the stdout of execv() to file outputData
            char outputData[64];
            sprintf(outputData, "outputData_%i.txt", id);
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
    
    free(buffer);
    free(execClient);
    free(execServer);
    free(directory);
    
    for (int n = 0; n < 32; n ++) {
        free(command[n]);
    }
    
    free(REMD_T);
    return 0;
}


