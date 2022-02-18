# todo: compile sources all (non-main) sources into objects, look into grouped targets

INCS= -I../inc -I../PackedArray
LINKS= -lpthread -lm 
CFLAGS= -g -fopenmp -Wall -Wextra -Wshadow -std=c99 -Wpedantic -O3 -Wstrict-overflow -fno-strict-aliasing
CC= gcc

# -DDEBUG_ASSERT_ON if you want to *slowly* double check the sieving results and some misc. other asserts
# CFLAGS += -DDEBUG_ASSERT_ON

cli: src/pom_yang_cli_main.c obj/moews_moews_sieve.o obj/sumdiv.o obj/pom_yang.o obj/PackedArray.o
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

test: src/pom_yang_test_main.c obj/moews_moews_sieve.o obj/sumdiv.o obj/pom_yang.o obj/PackedArray.o obj/bruteforce_kparent.o
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS) 

profile: src/profile_pom_yang.c obj/moews_moews_sieve.o obj/sumdiv.o obj/pom_yang.o obj/PackedArray.o
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

obj/pom_yang.o: src/pom_yang.c
	$(CC) -c $(CFLAGS) $^ -o $@ 

obj/PackedArray.o: ../PackedArray/PackedArray.c
	$(CC) -c -O3 $^ -o $@ 

obj/sumdiv.o: src/sumdiv.c
	$(CC) -c $(CFLAGS) $^ -o $@

obj/moews_moews_sieve.o: src/moews_moews_sieve.c
	$(CC) -c $(CFLAGS) $^ -o $@

obj/bruteforce_kparent.o: src/bruteforce_kparent.c
	$(CC) -c $(CFLAGS) $^ -o $@

.PHONY: clean
clean: 
	rm obj/*