#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <arpa/inet.h>

//#define PRTDETAIL
#define BACKLOG 15         /* buffer for waiting connections */
#define TIMEOUT 3600       /* exit if conect’s aren’t made in this #seconds */

#define min(x, y) (x < y) ? x : y

struct REMDTempStr {
    int num;
    double T;
};

double **probability;

void server(int *socket, int number_of_connects);
void collect(int *sockets, double *Epot, struct REMDTempStr *T, int number_of_connects);
void rearrange(double *Epot, struct REMDTempStr *T, int number_of_connects);
void swap(struct REMDTempStr *x, struct REMDTempStr *y);
void distribute(int *sockets, struct REMDTempStr *T, int number_of_connects);
int open_connections(int *sockets, int connects_wanted, int port_number, int timeout);
int close_connections(int *sockets, int connections_open);


/*
This is the server code which communicates through the
specified port number (essentially arbitrary except
low numbers are reserved) and expects a certain number
of clients to connect.
*/
int main(int argc, char *argv[]) {
	int *sockets;
	int port_number;
	int numconnects;
	int connects_made;
	long seed;

	if (argc != 4) {
		fprintf(stderr, "!usage: %s [port_number] [#_of_connections] [seed]\n", argv[0]);
		exit(0);
	}

	sscanf(argv[1], "%d", &port_number);
	sscanf(argv[2], "%d", &numconnects);
	sscanf(argv[3], "%ld", &seed);
	srand48(seed);

	sockets = calloc(numconnects, sizeof(int));
	if (!sockets) {
		perror("calloc");
		exit(1);
	}

	if (numconnects != (connects_made = open_connections(sockets, numconnects, port_number, TIMEOUT))) {
		perror("open_connections");
		fprintf(stderr, "! Only %d connections made instead of %d\n", connects_made, numconnects);
		exit(1);
	}
    
    probability = (double **)calloc(numconnects + 1, sizeof(double *));
    for (int i = 0; i <= numconnects; i ++) {
        probability[i] = (double *)calloc(2, sizeof(double));
    }
    
	server(sockets, numconnects);
    
    printf("\ntotal exchange rate = %.2lf%%\n", probability[numconnects][1] / probability[numconnects][0] * 100);
    for (int i = 0; i < numconnects - 1; i ++) {
        printf("%i <-> %i: %.2lf%%\n", i + 1, i + 2, probability[i][1] / probability[i][0] * 100);
    }
    
    for (int i = 0; i <= numconnects; i ++) {
        free(probability[i]);
    }
    free(probability);

	return 0;
}

void server(int *sockets, int number_of_connects) {
	double *Epot;
    struct REMDTempStr *T;

	if ((Epot = calloc(number_of_connects, sizeof(double))) == 0) {
		perror("calloc");
		exit(1);
	}

	if ((T = calloc(number_of_connects, sizeof(struct REMDTempStr))) == 0) {
		perror("calloc");
		exit(1);
	}

	collect(sockets, Epot, T, number_of_connects);
	while (T[0].T != 0.0) {
		rearrange(Epot, T, number_of_connects);
		distribute(sockets, T, number_of_connects);
		collect(sockets, Epot, T, number_of_connects);
	}
	close_connections(sockets, number_of_connects);

	free(Epot);
	free(T);
}

void collect(int *sockets, double *Epot, struct REMDTempStr *T, int number_of_connects) {
	static int i;
	static long recv_size;

	for (i = 0; i < number_of_connects; i++) {
		recv_size = read(sockets[i], &Epot[i], sizeof(double));
		if (recv_size != sizeof(double)) {
			perror("recv");
			printf("%s:%i\n", __FILE__, __LINE__);
			exit(1);
		}
	}

	for (i = 0; i < number_of_connects; i++) {
		recv_size = read(sockets[i], &T[i], sizeof(struct REMDTempStr));
		if (recv_size != sizeof(struct REMDTempStr)) {
			perror("recv");
			printf("%s:%i\n", __FILE__, __LINE__);
			exit(1);
		}
	}
}

void rearrange(double *Epot, struct REMDTempStr *T, int number_of_connects) {
	static int i;
	static int first, second;
	static double dE;
    double beta1, beta2;
    double record = 0;

	for (i = 0; i < number_of_connects * number_of_connects; i++) {
		first = drand48() * number_of_connects;
		second = drand48() * number_of_connects;
        if (first == second) {
            continue;
        }
        
        beta1 = 1 / T[first].T;
        beta2 = 1 / T[second].T;
        
		dE = (beta2 - beta1) * (Epot[first] - Epot[second]);

        if (first == second - 1 || first == second + 1) {
            probability[min(first, second)][0] ++;
        }
        probability[number_of_connects][0] ++;
        
#ifdef PRTDETAIL
        printf("first = %2i, second = %2i, Epot[%i] = %6.2lf, Eport[%i] = %6.2lf, T[%i] = %6.2lf, T[%i] = %6.2lf, dE = %8.4lf, exp = %8.4lf ",
               first, second, first, Epot[first], second, Epot[second], first, T[first], second, T[second], dE, exp(-1 * dE));
#endif
        
        if (dE <= 0.0 || (record = drand48()) < exp(-1 * dE)) {
            
#ifdef PRTDETAIL
            printf("rand = %8.4lf, accept.\n", record);
#endif
            
			swap(T + first, T + second);

            if (first == second - 1 || first == second + 1) {
                probability[min(first, second)][1] ++;
            }
            probability[number_of_connects][1] ++;
            
        } else {
#ifdef PRTDETAIL
            printf("rand = %8.4lf, reject.\n", record);
#endif
        }
	}
    
#ifdef PRTDETAIL
    printf("\n\n");
#endif
}

void swap(struct REMDTempStr *x, struct REMDTempStr *y) {
	static struct REMDTempStr z;
	z = *x;
	*x = *y;
	*y = z;
}

void distribute(int *sockets, struct REMDTempStr *T, int number_of_connects) {
	static int i;
	static long size;

	for (i = 0; i < number_of_connects; i++) {
		if ((size = write(sockets[i], &T[i], sizeof(struct REMDTempStr))) != sizeof(struct REMDTempStr)) {
			perror("send");
			exit(1);
		}
	}
}


int open_connections(int *sockets, int connects_wanted, int port_number, int timeout) {
    unsigned int sin_size;
    int sock_master;
    time_t current_time;
    int connects_made = 0;
    struct sockaddr_in my_address, *their_address;
    
    if ((their_address = calloc(connects_wanted, sizeof(struct sockaddr_in))) == 0) {
        perror("open_connections");
        exit(1);
    }
    
    if ((sock_master = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
        perror ("socket");
        exit(1);
    }
    
    bzero((char *) &my_address, sizeof(my_address));
    my_address.sin_family = AF_INET;
    my_address.sin_port = htons(port_number);
    my_address.sin_addr.s_addr = INADDR_ANY;
    
    if (bind(sock_master, (struct sockaddr *)&my_address, sizeof(struct sockaddr)) == -1) {
        perror ("bind");
        exit(1);
    }
    
    if (listen(sock_master, BACKLOG) == -1) {
        perror ("listen");
        exit(1);
    }
    
    current_time = time(0);
    sin_size = sizeof(struct sockaddr_in);
    do {
        if ((sockets[connects_made] = accept(sock_master, (struct sockaddr *)&(their_address[connects_made]), &sin_size)) != -1) {
            ++connects_made;
            printf("connection %i established\n", connects_made);
        }
    } while ((time(0) < current_time + timeout) && (connects_made < connects_wanted));
    
    free(their_address);
    return connects_made;
}

int close_connections(int *sockets, int connections_open) {
    int connects_closed = 0;
    int i;
    
    for (i = 0; i < connections_open; i++)
        connects_closed += !close(sockets[i]);
    return connects_closed = connections_open;
}
