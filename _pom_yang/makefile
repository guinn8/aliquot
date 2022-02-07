SRCDIR=../src
INCS=-I../inc -I../PackedArray
LINKS= -lpthread -lm 
CFLAGS=-g -fopenmp -Wall -Werror -Wextra -Wshadow -std=c99 -Wpedantic -O3
CC=gcc

bin/pom_yang_cli: bin/pom_yang.o bin/PackedArray.o pom_yang_cli.c $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

# -DDEBUG_ASSERT_ON if you want to *slowly* double check the sieving results
bin/pom_yang.o: pom_yang.c 
	$(CC) -c $(CFLAGS) $^ -o $@ 

bin/PackedArray.o: ../PackedArray/PackedArray.c
	$(CC) -c -O3 $^ -o $@ 

test: bin/pom_yang_test
bin/pom_yang_test: pom_yang_test.c bin/pom_yang.o bin/PackedArray.o $(SRCDIR)/sieve_rewrite.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

.PHONY: clean
clean: 
	rm bin/*