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
      REMD.o \
      Visual.o

# Linux (default)
LDFLAGS = -lGL -lGLU -lglut

# macOS
ifeq ($(shell uname -s), Darwin)
    LDFLAGS = -framework OpenGL -framework GLUT -framework Cocoa
endif

#Windows (cygwin)
ifeq ($(OS), Windows_NT)
	LDFLAGS = -lopengl32 -lglu32 -lglut32
endif

all: sDMD sServer sREMD

sServer: Server.o
	$(CC) -o $@ $^ -lm

sREMD: REMDExec.o
	$(CC) -o $@ $^

sDMD: $(OBJ)
	$(CC) -o $@ $^ -lm -pthread $(LDFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) $^

clean:
	rm -f *.o *.h.gch sDMD sServer sREMD
