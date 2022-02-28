SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lgmp
CFLAGS=-g -fopenmp -Wall -Wextra -Wshadow -std=c99 -Wpedantic -O3 -Wstrict-overflow -fno-strict-aliasing

huer_model: huer_model.c ../sumdiv/sumdiv.c ../sumdiv/factor.o
	gcc $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

../sumdiv/factor.o: ../sumdiv/factor.c
	gcc  -c $(INCS) -o $@ $^ $(LINKS)

.PHONY: clean
clean: 
	rm huer_model