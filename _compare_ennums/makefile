CFLAGS = -g -Wall -fopenmp -lm -O3

all: sn_sieve sn_accum

sn_sieve: sn_sieve.c sieve.c
	gcc -o $@ $^ $(CFLAGS)

sn_accum: sn.c
	gcc -o $@ $^ $(CFLAGS)

test: sn_sieve sn_accum
	./test.sh