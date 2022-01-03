SRCDIR=../src
INCS=-I../inc
LINKS=-lflint -lmpfr -lgmp -lpthread -lm
CFLAGS=-g -fopenmp -Wall 


main: main.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c 
	gcc $(CFLAGS) -O3 -DNDEBUG $(INCS) -o $@ $^ $(LINKS)

debug_main: main.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c 
	gcc $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

test: debug_main
	./test.sh

.PHONY: clean
clean: 
	rm main