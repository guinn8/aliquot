NAME=py
CC=gcc
#orignially compiled with icc
CFLAGS= -fopenmp -g 
SYMBS=
INCS=
LIBS=
LINKS=-lm -lflint -lmpfr -lgmp -lpthread
all: py.c sieve.c mmapArray.c properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o $(NAME) $^ $(LIBS) $(LINKS)
lib:
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c init.c -o init.o
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c mult.c -o mult.o
	ar rcs lib$(NAME).a init.o mult.o
clean:
	rm *.o


