# todo: compile sources all (non-main) sources into objects, look into grouped targets

INCS= -I../inc -I../PackedArray
LINKS= -lpthread -lm 
CFLAGS= -g -fopenmp -Wall -Werror -Wextra -Wshadow -std=c99 -Wpedantic -O3
CC= gcc

cli: src/pom_yang_cli_main.c src/sieve_rewrite.c src/properSumDiv.c obj/pom_yang.o obj/PackedArray.o
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

test: src/pom_yang_test_main.c src/sieve_rewrite.c src/properSumDiv.c src/brute_force_preimages.c obj/pom_yang.o obj/PackedArray.o
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

# -DDEBUG_ASSERT_ON if you want to *slowly* double check the sieving results and some misc. other asserts
obj/pom_yang.o: src/pom_yang.c 
	$(CC) -c $(CFLAGS) $^ -o $@ 

obj/PackedArray.o: ../PackedArray/PackedArray.c
	$(CC) -c -O3 $^ -o $@ 


.PHONY: clean
clean: 
	rm bin/*