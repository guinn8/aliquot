NAME=main.exe
CC=gcc
#orignially compiled with icc
CFLAGS= 
SYMBS= -std=c99 -O3 -I../shared/include/
INCS=
LIBS=
LINKS=-lm
all: main.c ../shared/src/properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o $(NAME) $^ $(LIBS) $(LINKS)
#lib:
#	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c init.c -o init.o
#	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -c mult.c -o mult.o
#	ar rcs lib$(NAME).a init.o mult.o
clean:
	rm *.o


