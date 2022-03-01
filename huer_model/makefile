SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lgmp
CFLAGS=-g -fopenmp -Wall -Wextra -Wshadow -std=c99 -Wpedantic -O3 -Wstrict-overflow -fno-strict-aliasing



huer_model: huer_model.c ../factor/factor.c
	gcc $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)


.PHONY: clean
clean: 
	rm huer_model