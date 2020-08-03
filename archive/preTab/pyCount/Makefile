NAME=pyCount
CC=gcc
#orignially compiled with icc
CFLAGS=-Wall -fopenmp -std=gnu99 -O3 -g
SYMBS=-DBITSIZE=64 -DIMAGESIZE=64
INCS= ../config/config.c
LIBS=
LINKS=-lm -lconfig -lpari
all: pyCount.c sieve.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o $(NAME) $^ $(LIBS) $(LINKS)
lib:
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c init.c -o init.o
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c mult.c -o mult.o
	ar rcs lib$(NAME).a init.o mult.o
clean:
	rm *.o
