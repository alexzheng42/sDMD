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
#include <signal.h>

#define CHAR_LENGTH 1024
#define BOLTZMANN 0.0019872041 //kcal / mol / K

#define CHAR_LOC(name) \
        name = (char *)calloc(CHAR_LENGTH, sizeof(char));

struct REMDTempStr {
    int num;
    double T;
};

int main (int argc, char *argv[]) {
    pid_t child_pid, server_pid;
    int next = 0;
    int status = 0;
    int replicaNum = 0;
    char portNo[32], exRate[32];
    char *buffer, *execClient, *execServer;
    char *directory;
    char *args[64];      // sDMD arguments given after -args (copied without length limits)
    int argNum = 0;
    struct REMDTempStr thisT;
    struct REMDTempStr *REMD_T;
    FILE *configFile;

    CHAR_LOC(buffer);
    CHAR_LOC(execClient);
    CHAR_LOC(execServer);
    CHAR_LOC(directory);

    if (getcwd(directory, CHAR_LENGTH) == NULL)
        perror("!!ERROR!!");
    char workDir[CHAR_LENGTH];
    snprintf(workDir, CHAR_LENGTH, "%s", directory);   // 'directory' is reused while reading the config

    for (int n = 1; n < argc; n ++) {
        if (argv[n][0] == '-') {
            if (strcmp(argv[n], "-h") == 0 || strcmp(argv[n], "-help") == 0) {
            help:
                printf("Usage: %s [distribute T or not] -f [configuration file] -args [args of executable without flag -REMD]\n\n", argv[0]);
                printf("            -f: provide the configuration file.\n");
                printf("         -args: provide the flags and info for sDMD without -REMD, end with \"-end\"\n");
                exit(EXIT_SUCCESS);
            } else if (strcmp(argv[n], "-f") == 0 && n + 1 < argc) {
                snprintf(buffer, CHAR_LENGTH, "%s", argv[n + 1]);
            } else if (strcmp(argv[n], "-args") == 0) {
                next = 1;
                while (n + next < argc && strcmp(argv[n + next], "-end") && argNum < 50) {
                    args[argNum++] = argv[n + next];
                    next ++;
                }
                n += next;
            }
        }
    }
    if (argNum == 0) goto help;

    // sDMD and sServer are expected next to the working directory, as before
    snprintf(execClient, CHAR_LENGTH, "%s/sDMD", directory);
    snprintf(execServer, CHAR_LENGTH, "%s/sServer", directory);

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

    // the server runs as a child too, so this process can supervise all of them
    char repNum[20]; sprintf(repNum, "%i", replicaNum);
    char seed[20];   sprintf(seed, "%u", (unsigned)time(NULL));
    char *serverCommand[] = {"./sServer", portNo, repNum, seed, (char *)0};
    if ((server_pid = fork()) == 0) {
        if (execv(execServer, serverCommand) < 0)
            perror("!!ERROR!!: cannot start sServer");
        exit(EXIT_FAILURE);
    }

    pid_t *replicaPid = (pid_t *)calloc(replicaNum, sizeof(pid_t));
    for (int id = 0; id < replicaNum; id++) {
        thisT = REMD_T[id];

        // argv of sDMD: program name, the user's arguments, then -REMD and its values
        char values[6][64];
        char *command[64 + 8];
        int c = 0;
        command[c++] = execClient;
        for (int a = 0; a < argNum; a++) command[c++] = args[a];
        sprintf(values[0], "localhost");
        sprintf(values[1], "%s" , portNo);
        sprintf(values[2], "%lf", thisT.T);
        sprintf(values[3], "%i" , thisT.num);
        sprintf(values[4], "%s" , exRate);
        sprintf(values[5], "%i" , id);
        command[c++] = "-REMD";
        for (int v = 0; v < 6; v++) command[c++] = values[v];
        command[c] = (char *)0;

        if ((child_pid = fork()) == 0) {
            //direct the stdout of execv() to file outputData
            char outputData[64];
            sprintf(outputData, "outputData_%i.txt", id);
            int fd = open(outputData, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
            dup2(fd, 1);
            close(fd);

            sleep(1 + id);
            if (execv(execClient, command) < 0)
                perror("!!ERROR!!: cannot start sDMD");
            exit(EXIT_FAILURE);
        }
        replicaPid[id] = child_pid;
    }

    // wait for the replicas; if they are all gone but the server is still waiting for them,
    // something went wrong (e.g. a replica failed at start-up): stop the server instead of
    // letting it wait for its connection timeout
    int running = replicaNum, failed = 0, serverDone = 0;
    while (running > 0 || !serverDone) {
        pid_t wpid = waitpid(-1, &status, running > 0 ? 0 : WNOHANG);
        if (wpid == server_pid) {
            serverDone = 1;
        } else if (wpid > 0) {
            running--;
            if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) failed++;
        } else if (running == 0) {
            int waited = 0;
            while (waited < 30 && waitpid(server_pid, &status, WNOHANG) == 0) { sleep(1); waited++; }
            if (waited >= 30) {
                printf("!!ERROR!!: all replicas have exited but sServer is still waiting; stopping it.\n"
                       "           Check outputData_*.txt in %s for the replica output.\n", workDir);
                kill(server_pid, SIGTERM);
                waitpid(server_pid, &status, 0);
                failed++;
            }
            serverDone = 1;
        }
    }
    if (failed) {
        printf("!!ERROR!!: %i replica(s) or the server did not finish normally; see outputData_*.txt.\n", failed);
    }

    free(replicaPid);
    free(buffer);
    free(execClient);
    free(execServer);
    free(directory);
    free(REMD_T);
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
