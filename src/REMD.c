//
//  REMD.c
//  sDMD
//
//  Created by Size Zheng on 11/9/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h> /* inet_addr function */
#include <netinet/in.h>
#include <unistd.h>    /* close function */
#include <netdb.h>
#include <strings.h>
#include "DMD.h"

#define h_addr h_addr_list[0] /* for backward compatibility */

int client_open(char *serverName, int port_number);
void ReplicaExchange(double *E, double *T, char *serverName, int *port_number);
void REMDRun(struct ThreadInfoStr *threadInfo, int outputRate, char *serverName, int *port_number);


void REMD(void) {
    if (REMDInfo.REMD_T > 0) {
        targetTemperature = REMDInfo.REMD_T;
    }
    
    thread = (struct ThreadStr **)calloc(1, sizeof(struct ThreadStr *));
    thread[0] = InitializeThread(0, atom);
    thread[0]->fileList = InitializeFiles(REMDInfo.REMD_ExtraName, fileList);
    SaveData(lgf, thread[0]);
    
    if (strncmp("new", neworcontinue, 1) == 0) {
        SaveData(sysInfo, thread[0]);
        FirstRun(thread[0]);
    } else if (strncmp("continue", neworcontinue, 1) == 0) {
        outputrecord = currenttime + outputrate;
        CreateCBT();
    }
    
    thrInfo[0].threadRenewList[1] = eventToCommit = CBT.node[1];
    if (thread[0]->raw[eventToCommit]->dynamic->event.partner > 0) {
        thrInfo[0].threadRenewList[++thrInfo[0].threadRenewList[0]] = thread[0]->raw[eventToCommit]->dynamic->event.partner;
    }
    AssignThread(thrInfo[0].threadRenewList, thread[0]);
    REMDRun(&thrInfo[0], REMDInfo.REMD_OutputRate, REMDInfo.REMD_ServerName, &REMDInfo.REMD_PortNum);
}

void REMDRun(struct ThreadInfoStr *threadInfo, int outputRate, char *serverName, int *port_number) {
    int hazardType;
    int tid = threadInfo->threadID;
    int *threadRenewList = threadInfo->threadRenewList;
    int rate = ((int)(currenttime / outputRate) + 1) * outputRate;
    double timeIncr;
    double processratio, gap = 0;
    struct ThreadStr *thisThread = thread[tid];
    struct AtomStr **newTarget = &thisThread->newTarget;;
    struct AtomStr **newPartner = &thisThread->newPartner;
    struct AtomStr *oldTarget = &thisThread->oldTarget;
    struct AtomStr *oldPartner = &thisThread->oldPartner;
    
    while (currenttime <= timestep) {
        
        hazardType = ProcessEvent(threadRenewList, thisThread);
        
        if (hazardType < 0) {
            threadRenewList[0] = 1;
            threadRenewList[2] = 0;
            
            AssignThread(threadRenewList, thisThread);
            Predict(threadRenewList, thisThread);
            
            AtomDataCpy(thisThread->raw[threadRenewList[1]], thisThread->listPtr[threadRenewList[1]], 0);
            UpdateCBT(threadRenewList);
            
            AssignJob(threadRenewList, thisThread);
            AssignThread(threadRenewList, thisThread);
            
        } else {

            timeIncr = thisThread->raw[threadRenewList[1]]->dynamic->event.time;
            if (unlikely(timeIncr < 0)) {
                printf("!!ERROR!!: time calculation is not correct!\n");
            }
            
            UpdateData(timeIncr, "atom", thisThread); //update the coordinates
            TimeForward(timeIncr, thisThread); //update the node time
            currenttime += timeIncr;
            frame++;
            
            CommitEvent(thisThread->raw,
                        *newTarget, (threadRenewList[2] > 0 ? *newPartner : NULL),
                        oldTarget, (threadRenewList[2] > 0 ?  oldPartner : NULL),
                        thisThread->listPtr[oldTarget->dynamic->HB.neighbor],
                        (threadRenewList[2] > 0 ? thisThread->listPtr[oldPartner->dynamic->HB.neighbor] : NULL));
            
            UpdateCBT(threadRenewList);
            processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
            
            if (processratio >= gap + 0.01) {
                
                printf("Process=%8.2lf%%\r", processratio);
                fflush(stdout);
                gap += 0.01;
                
            }
            
            //save data to the output files
            if (currenttime >= outputrecord) {
                for (int i = 0; i < lenFileType; i ++) {
                    if (thisThread->fileList[i].mark && i != RE) {
                        SaveData(i, thisThread);
                    }
                }
                outputrecord += outputrate;
            }
            
            if (currenttime >= rate) {
                double potenEnergy = CalPotenE(thisThread);
                ReplicaExchange(&potenEnergy, &targetTemperature, serverName, port_number);
                SaveData(RE, thisThread);
                
                rate += outputRate;
            }
            
            AssignJob(threadRenewList, thisThread);
            AssignThread(threadRenewList, thisThread);
        }
    }
    
    tmpDouble = 0.0;
    ReplicaExchange(&tmpDouble, &tmpDouble, serverName, port_number);
    
    return;
}


/*
 This subroutine makes an attempted exchange from
 the current replica with energy *E and temperature
 *T through the server of name *serverName
 and port number given by *port_number
 */
void ReplicaExchange(double *E, double *T, char *serverName, int *port_number) {
    long len;
    double energy, temperature;
    
    static int is_initialized = 0;
    static int sock;
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
     First time called? Then open a socket
     *- - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    if (!is_initialized) {
        sock = client_open(serverName, *port_number);
        is_initialized = 1;
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - -*
     communicate with the server: write
     *- - - - - - - - - - - - - - - - - - - - - - - - -*/
    energy = *E;
    temperature = *T;
    
    len = write(sock, &energy, sizeof(double));
    if (unlikely(len != sizeof(double))) {
        perror("send");
        exit(1);
    }
    
    len = write(sock, &temperature, sizeof(double));
    if (unlikely(len != sizeof(double))) {
        perror("send");
        exit(1);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*
     if energy=temperature=0 close socket
     *- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    if (unlikely(*E == 0 && *T == 0)) {close(sock); return;};
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*
     receive new temperature
     *- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    temperature = 0;
    len = read(sock, &temperature, sizeof(double));
    if (unlikely(len != sizeof(double))) {
        perror("recv");
        printf("%s:%i\n", __FILE__, __LINE__);
        exit(1);
    }
    
    *T = temperature;
    
    return;
}


int client_open(char *serverName, int port_number) {
    static int sock;
    static struct sockaddr_in server_addr;
    static int error_flag;
    static struct hostent *server;
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - -*
     open an IPv4 socket of type SOCK_STREAM
     *- - - - - - - - - - - - - - - - - - - - - - - - -*/
    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == -1) {
        perror("socket");
        exit(1);
    }
    server = gethostbyname(serverName);
    if (server == NULL) {
        perror("ERROR, no such host");
        exit(1);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - -*
     set the destination address
     *- - - - - - - - - - - - - - - - - - - - - - - -*/
    bzero((char *) &server_addr, sizeof(server_addr));
    server_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr, (char *)&server_addr.sin_addr.s_addr, server->h_length);
    server_addr.sin_port = htons(port_number);
    
    /*- - - - - - - - - - - - - - - - - - - - - - - -*
     connect to server
     *- - - - - - - - - - - - - - - - - - - - - - - -*/
    error_flag = connect(sock, (struct sockaddr*)&server_addr, sizeof(struct sockaddr));
    if (error_flag == -1) {
        perror("connect");
        exit(1);
    }
    
    /*- - - - - - - - - - - - - - - - - -*
     return the open socket
     *- - - - - - - - - - - - - - - - - -*/
    return sock;
}
