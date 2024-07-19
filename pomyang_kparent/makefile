INCS= -I../inc
LINKS= -lpthread -lm -lgmp
CFLAGS= -g -fopenmp -Wall -Wextra -Wshadow -std=c99 -Wpedantic -O3 -Wstrict-overflow -fno-strict-aliasing
CC= gcc

# -DDEBUG_ASSERT_ON if you want to *slowly* double check the sieving results and some misc. other asserts
# CFLAGS += -DDEBUG_ASSERT_ON
all: cli test profiler



cli: src/pomyang_cli_main.c obj/moewsmoews_sieve.o obj/factor.o obj/pomyang_kparent.o obj/PackedArray.o
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS)

test: src/pomyang_test_main.c obj/moewsmoews_sieve.o obj/factor.o obj/pomyang_kparent.o obj/PackedArray.o obj/bruteforce_kparent.o
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS) 

profiler: src/pomyang_profile_main.c obj/moewsmoews_sieve.o obj/factor.o obj/pomyang_kparent.o obj/PackedArray.o
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS)

geomean_sn: src/geomean_sn.c obj/moewsmoews_sieve.o obj/factor.o obj/pomyang_kparent.o obj/PackedArray.o
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS)

readfile: src/pomyang_readfile.c
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS)


obj/pomyang_kparent.o: src/pomyang_kparent.c
	mkdir -p obj
	$(CC) -c $(CFLAGS) $^ -o $@ 

obj/PackedArray.o: src/PackedArray.c
	mkdir -p obj
	$(CC) -c -O3 $^ -o $@ 

obj/factor.o: ../factor/factor.c
	mkdir -p obj
	$(CC) -c $(CFLAGS) $^ -o $@

obj/moewsmoews_sieve.o: src/moewsmoews_sieve.c
	mkdir -p obj
	$(CC) -c $(CFLAGS) $^ -o $@

obj/bruteforce_kparent.o: src/bruteforce_kparent.c
	mkdir -p obj
	$(CC) -c $(CFLAGS) $^ -o $@

.PHONY: clean
clean: 
	rm obj/*
	rm bin/*