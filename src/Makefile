#Defines the compiler
CC = g++
#Displays all compilation warnings
CFLAGS += -Wall
#Optimize the code for speed
CFLAGS += -O3 -ffast-math
#Enable debuging with gdb (makes the executable bigger)
#CFLAGS += -ggdb
#Enable gprof based profiling
CFLAGS += -pg
#Allow to use intrinsic functions
CFLAGS += -march=native
#On SSE2 compatible processors this will compile a faster program, in float 
#precision. Comment for a slower program in double precision.
DEFINES += -DVECTORIZE

all : findMinimum

findMinimum : octtree.o mycode.o common.o main.o
	$(CC) $(CFLAGS) $(DEFINES) $^ -o $@ && mv $@ ../bin/

octtree.o mycode.o common.o main.o : %.o : %.cpp
	$(CC) $(CFLAGS) $(DEFINES) -c $<

clean :
	rm -rf *.o
