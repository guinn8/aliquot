SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread
CFLAGS=-fopenmp -std=c99  -O3 


pom_yang_preimage: pom_yang_preimage.c $(SRCDIR)/properSumDiv.c $(SRCDIR)/sieve.c 
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

.PHONY: clean
clean: 
	rm pom_yang_preimage