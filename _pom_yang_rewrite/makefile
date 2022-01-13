SRCDIR=../src

INCS=-I../inc -I../PackedArray
LINKS= -lpthread -lm 
CFLAGS=-g -fopenmp -Wall 


main: main.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c ../PackedArray/PackedArray.c
	gcc $(CFLAGS) -O3 -DNDEBUG $(INCS) -o $@ $^ $(LINKS)

debug_main: main.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c ../PackedArray/PackedArray.c
	gcc $(CFLAGS) $(INCS) -o $@ $^ $(LINKS) -lflint -lmpfr -lgmp 

test: debug_main
	./test.sh

.PHONY: clean
clean: 
	rm main