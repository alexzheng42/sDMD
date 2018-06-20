CC = gcc
CFLAGS = -I. -O2 -std=c99
DEPS = List.h DMD.h
OBJ = CBT.o \
      DataSave.o \
      DMD.o \
      Events.o \
      Initialization.o \
      Models.o \
      ThreadProcess.o \
      TimePrediction.o \
      ToolFunctions.o \
      List.o \
      MSThread.o \
      SGThread.o \
      REMD.o

all: sDMD sServer sREMD

sServer: Server.o
	$(CC) -o $@ $^ -lm

sREMD: REMDExec.o
	$(CC) -o $@ $^

sDMD: $(OBJ)
	$(CC) -o $@ $^ -lm -pthread

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) $^

clean:
	rm -f *.o *.h.gch sDMD sServer sREMD
