SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread
CFLAGS=-fopenmp -std=c99  -O3 

huer_model: huer_model.c $(SRCDIR)/sieve.c  $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

.PHONY: clean
clean: 
	rm huer_model