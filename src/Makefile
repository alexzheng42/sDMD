CC = gcc
CFLAGS = -I. -O2 -std=c99 -w
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

VISOBJ = Visual.o
VISFLAG = -D VIS
VIS ?= no

# Linux (default)
LDFLAGS = -lGL -lGLU -lglut

# macOS
ifeq ($(shell uname -s), Darwin)
	LDFLAGS = -framework OpenGL -framework GLUT -framework Cocoa
endif

#Windows (cygwin)
ifeq ($(OS), Windows_NT)
	LDFLAGS = -lgl -lglu -lglut
endif

all: sDMD sServer sREMD

sServer: Server.o
	$(CC) -o $@ $^ -lm

sREMD: REMDExec.o
	$(CC) -o $@ $^

#enable visual
ifeq ($(VIS), yes)
sDMD: $(OBJ) $(VISOBJ)
	$(CC) -o $@ $^ -lm -pthread $(LDFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) $^ $(VISFLAG)

else

#disable visual
sDMD: $(OBJ)
	$(CC) -o $@ $^ -lm -pthread

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) $^
endif

clean:
	rm -f *.o *.h.gch sDMD sServer sREMD
