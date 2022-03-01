CFLAGS = -g -Wall -fopenmp -lm -O3 #-ftree-vectorize -fopt-info-vec-all 
SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread

all: cli

cli: sieve_sn.c sieve_sn_cli.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)
