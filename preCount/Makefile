NAME=main.exe
CC=gcc
#orignially compiled with icc
CFLAGS= -fopenmp -std=c99 -g -O2 -I../shared/include/
#-fOpenMP: Threading
#-std=c99: loop problem on some machines

SYMBS=
INCS=
LIBS=
LINKS=-lm -L/home/gavin.guinn/flint-2.6.2 -lflint -lmpfr -lgmp -lpthread
all: main.c ../shared/src/sieve.c  ../shared/src/properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o $(NAME) $^ $(LIBS) $(LINKS)
#lib:
#	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c init.c -o init.o
#	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c mult.c -o mult.o
#	ar rcs lib$(NAME).a init.o mult.o
clean:
	rm *.o


