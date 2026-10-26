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
#include "DMD.h"

#define h_addr h_addr_list[0] /* for backward compatibility */

int client_open(char *serverName, int port_number);
void ReplicaExchange(double *E, struct REMDTempStr *T, char *serverName, int *port_number);

// replica exchange every REMD_OutputRate time units, hooked into the common event loop
static double nextExchange;

static void REMDExchange(void) {
    if (currenttime < nextExchange) return;
    UpdateAllAtoms(currenttime);
    double potenEnergy = CalPotenE();
    ReplicaExchange(&potenEnergy, &REMDInfo.REMD_Temperature, REMDInfo.REMD_ServerName, &REMDInfo.REMD_PortNum);
    targetTemperature = REMDInfo.REMD_Temperature.T;
    SaveData(RE, fileList);
    nextExchange += REMDInfo.REMD_OutputRate;
}

// tell the server this replica has finished (energy = temperature = 0 closes the socket)
static void REMDFinish(void) {
    tmpDouble = 0.0;
    REMDInfo.REMD_Temperature.T = 0.0;
    ReplicaExchange(&tmpDouble, &REMDInfo.REMD_Temperature, REMDInfo.REMD_ServerName, &REMDInfo.REMD_PortNum);
}

void REMD(void) {
    Prepare();
    int rate = REMDInfo.REMD_OutputRate;
    nextExchange = ((int)(currenttime / rate) + 1) * rate;
    SingleThreadRun(REMDExchange, REMDFinish);
}


/*
 This subroutine makes an attempted exchange from
 the current replica with energy *E and temperature
 *T through the server of name *serverName
 and port number given by *port_number
 */
void ReplicaExchange(double *E, struct REMDTempStr *T, char *serverName, int *port_number) {
    long len;
    double energy;
    struct REMDTempStr thisT;
    
    /*
     The local static variable will be initiated only at the first call;
     its value during the last call will be kept.
     */
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
    thisT = *T;
    
    len = write(sock, &energy, sizeof(double));
    if (unlikely(len != sizeof(double))) {
        perror("send");
        exit(1);
    }
    
    len = write(sock, &thisT, sizeof(struct REMDTempStr));
    if (unlikely(len != sizeof(struct REMDTempStr))) {
        perror("send");
        exit(1);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*
     if energy=temperature=0 close socket
     *- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    if (unlikely(*E == 0 && (*T).T == 0)) {close(sock); return;};
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*
     receive new temperature
     *- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    thisT.T = 0;
    len = read(sock, &thisT, sizeof(struct REMDTempStr));
    if (unlikely(len != sizeof(struct REMDTempStr) && thisT.T)) {
        perror("recv");
        printf("%s:%i\n", __FILE__, __LINE__);
        exit(1);
    }
    
    *T = thisT;
    
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
    memset((char *) &server_addr, 0, sizeof(server_addr));
    server_addr.sin_family = AF_INET;
    memcpy((char *)&server_addr.sin_addr.s_addr, (char *)server->h_addr, server->h_length);
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
