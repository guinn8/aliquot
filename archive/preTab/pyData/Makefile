NAME=pyData
CC=gcc
CFLAGS= -g 
SYMBS=
INCS= ../config/config.c
LIBS=
LINKS=-lm -lconfig
all: pyData.c 
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o $(NAME) $^ $(LIBS) $(LINKS)
clean:
	rm *.o
