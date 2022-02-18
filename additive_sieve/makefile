CFLAGS = -g -Wall -fopenmp -lm -O3 #-ftree-vectorize -fopt-info-vec-all 
SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread

all: sn_test

sn_test: sn.c sn_test.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)
