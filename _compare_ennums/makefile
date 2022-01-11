SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread
CFLAGS=-fopenmp -std=c99  -O3 -fopenmp-simd -lflint -lmpfr -lgmp -lpthread

all: sn_sieve sn_accum

sn_sieve: sn_sieve.c $(SRCDIR)/sieve.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

sn_accum: sn.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

test: sn_sieve sn_accum
	./test.sh


.PHONY: clean
clean: 
	rm sn_accum
	rm sn_sieve