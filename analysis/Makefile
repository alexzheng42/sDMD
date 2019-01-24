CC = gcc
CFLAGS = -I. -Os -std=c99 -w
DEPS = Analysis.h
OBJ = Analysis.o 		  \
	  Cluster.o  		  \
	  Energy.o    		  \
	  FileManage.o 		  \
	  HBRamach.o   		  \
	  PBCAdjust.o  		  \
	  REMD.o       		  \
	  RG.o 				  \
	  RMSD.o              \
	  SystemInformation.o \
	  Tools.o

all: sAnalysis

sAnalysis: $(OBJ)
	$(CC) -o $@ $^ -lm

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) $^

clean:
	rm -f *.o *.h.gch sAnalysis
