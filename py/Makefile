NAME=py
CC=gcc
#orignially compiled with icc
CFLAGS= -fopenmp -std=c99 
SYMBS=
INCS=
LIBS=
LINKS=-lm
all: py.c sieve.c properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o $(NAME) $^ $(LIBS) $(LINKS)
lib:
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c init.c -o init.o
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c mult.c -o mult.o
	ar rcs lib$(NAME).a init.o mult.o
clean:
	rm *.o


