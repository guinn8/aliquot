SRCDIR=../src
INCS=-I../inc -I../PackedArray
LINKS= -lpthread -lm 
CFLAGS=-g -fopenmp -Wall -Werror -Wextra -Wshadow -std=c99 -Wpedantic 

main: main.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c PackedArray.o
	gcc $(CFLAGS) -O3 $(INCS) -o $@ $^ $(LINKS)

debug_main: main.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c  PackedArray.o
	gcc $(CFLAGS) $(INCS) -DDEBUG_ASSERT_ON -o  $@ $^ $(LINKS)

PackedArray.o: ../PackedArray/PackedArray.c
	gcc -c  $^

.PHONY: clean
clean: 
	rm main debug_main