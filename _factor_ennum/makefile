SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread
CFLAGS=-fopenmp -std=c99  -O3 

factor_ennum: factor_ennum.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

.PHONY: clean
clean: 
	rm factor_ennum