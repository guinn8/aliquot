SRCDIR=../src
INCS=-I../inc
LINKS= -lpthread -lm -lflint -lmpfr -lgmp -lpthread
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