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
	  Tools.o 			  \
	  Distance.o 		  \
	  PSD.o 			  \
	  Density.o 		  \

all: sAnalysis

sAnalysis: $(OBJ)
	$(CC) -o $@ $^ -lm

# $< compiles only the .c file; listing the header as a prerequisite still rebuilds on header edits
%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o *.h.gch sAnalysis