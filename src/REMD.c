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
void ReplicaExchange(double *E, struct REMDTempStr *T, char *serverName, int *port_number);
void REMDRun(int outputRate, char* serverName, int* port_number);

void REMD(void) {
	Prepare();
	REMDRun(REMDInfo.REMD_OutputRate, REMDInfo.REMD_ServerName, &REMDInfo.REMD_PortNum);
}

void REMDRun(int outputRate, char *serverName, int *port_number) {
	int hazardType;
	int renewList[256] = { 0 };
	int rate = ((int)(currenttime / outputRate) + 1) * outputRate;
	double processratio, gap = 0;
	struct AtomListStr* thisList = &atomList;
	struct AtomStr* target;
	struct AtomStr* partner;

	while (currenttime <= timestep) {

		AssignJob(&target, &partner, renewList, thisList);
		UpdateData(thisList->timeInr);
		TimeForward(thisList->timeInr); //update the node time
		currenttime += thisList->timeInr;

		hazardType = HazardCheck(target, partner,
								 thisList->ptr[target->dynamic->HB.neighbor],
								 (partner != NULL ? thisList->ptr[partner->dynamic->HB.neighbor] : NULL),
								 thisList);

#ifdef DEBUG_PRINTF
		printf("frame = %4li, target atom = %5i, partner = %5i, %6s:%6s, type = %2i",
			   frame, target->property->num, target->dynamic->event.partner,
			   target->property->name,
			   partner == NULL ? "NULL" : partner->property->name,
			   target->dynamic->event.eventType);
		fflush(stdout);
#endif

		if (!hazardType) {
#ifdef DEBUG_PRINTF
			printf(", executed!\n");
#endif
			frame++;
			DoEvent(thisList); //calculate the after-event velocities
			target->dynamic->event.counter++;
			if (renewList[0] == 2 && strcmp(partner->property->extraProperty[1], "_CONS"))
				partner->dynamic->event.counter++;

			//if the target atom is a fixing atom and the event is a CGEvent, then all the atoms in the neighboring cells will be renewed
			if (CG.mark && target->dynamic->event.eventType == CGSu_Event &&
				strcmp(target->property->extraProperty[1], "_CONS") == 0) {
				FixAssignRenewList(renewList, thisList);
			}

			Predict(renewList, thisList); //pre-predict the after-event event
			UpdateCBT(renewList);

			processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
			if (unlikely(processratio >= gap + 0.01)) {
				printf("Process=%8.2lf%%\r", processratio);
				fflush(stdout);
				gap += 0.01;
			}

			if (strcmp(wallDyn.mark, "no")) DoWallDyn();

			//save data to the output files
			if (unlikely(currenttime >= outputrecord)) {
				for (int i = 0; i < lenFileType; i++) {
					if (fileList[i].mark) {
						SaveData(i, fileList);
					}
				}
				outputrecord += outputrate;
			}

			if (currenttime >= rate) {
				double potenEnergy = CalPotenE();
				ReplicaExchange(&potenEnergy, &REMDInfo.REMD_Temperature, serverName, port_number);
				targetTemperature = REMDInfo.REMD_Temperature.T;
				SaveData(RE, fileList);
                
                rate += outputRate;
            }

		} else {
#ifdef DEBUG_PRINTF
			printf(", denied, self/partner changed!\n");
#endif
			renewList[0] = 1;
			renewList[2] = 0;
			renewList[3] = 0;

			Predict(renewList, thisList);
			UpdateCBT(renewList);
			countReCal++;
		}
	}

	tmpDouble = 0.0;
	REMDInfo.REMD_Temperature.T = 0.0;
	ReplicaExchange(&tmpDouble, &REMDInfo.REMD_Temperature, serverName, port_number);

	FreeVariables();
    return;
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
