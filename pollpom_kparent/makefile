SRCDIR=../src
INCS=-I../inc -I.
LINKS=-lm -lgmp
CFLAGS=-g -fopenmp -Wall -Wextra -Wshadow -std=c99 -Wpedantic -O3 -Wstrict-overflow -fno-strict-aliasing

all: bindir bin/pollpom_kparent_cli bin/pollpom_kparent_profile

bin/factor.o: ../factor/factor.c
	gcc -O3 -c $(INCS) -o $@ $^ $(LINKS)

bin/pollpom_kparent.o: pollpom_kparent.c
	gcc $(CFLAGS) -c $(INCS) -o $@ $^ $(LINKS)

bin/pollpom_kparent_cli: pollpom_kparent_cli.c bin/pollpom_kparent.o bin/factor.o 
	gcc $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

bin/pollpom_kparent_profile: pollpom_kparent_profile.c bin/pollpom_kparent.o bin/factor.o 
	gcc $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)

bindir:
	@mkdir -p bin

.PHONY: clean
clean: 
	rm -r bin

