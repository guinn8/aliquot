SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread
CFLAGS=-fopenmp -std=c99  -O3 


pom_yang_ennum: pom_yang_ennum.c $(SRCDIR)/properSumDiv.c $(SRCDIR)/sieve.c 
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)