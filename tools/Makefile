#Defines the compiler
CC = gcc
CFLAGS += -lm -std=c99
#Displays all compilation warnings
CFLAGS += -Wall -Wextra

all : depouille

depouille : depouille.o
	$(CC) $(CFLAGS) $^ -o $@ && mv $@ ../bin/
depouille.o : %.o : %.c
	$(CC) $(CFLAGS) -g -c $<

clean :
	rm -rf *.o
