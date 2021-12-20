
SRCDIR=src
CFLAGS=-fopenmp -std=c99  -O3 
CC=gcc
NAME=main.exe
SYMBS=
INCS=-Iinc -Iflint-2.7.1 -I.
LIBS= 
LINKS= -lm -lpthread

all: pre_count pre_ennum pre_est

pre_count: $(SRCDIR)/pre_count_main.c $(SRCDIR)/sieve.c  $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o bin/pre_count $^ $(LIBS) $(LINKS)

pre_ennum: $(SRCDIR)/pre_ennum_main.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o bin/pre_ennum $^ $(LIBS) $(LINKS)

pre_est: $(SRCDIR)/pre_est_main.c $(SRCDIR)/sieve.c  $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(SYMBS) $(INCS) -o bin/pre_est $^ $(LIBS) $(LINKS)

clean: 
	rm bin/*

.PHONY: all clean