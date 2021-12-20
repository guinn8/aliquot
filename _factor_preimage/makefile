SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread
CFLAGS=-fopenmp -std=c99  -O3 

factor_preimage: factor_preimage.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

.PHONY: clean
clean: 
	rm factor_preimage