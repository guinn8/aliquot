CFLAGS = -g -Wall -fopenmp -lm -O3 #-ftree-vectorize -fopt-info-vec-all 
SRCDIR=../src
INCS=-I../inc
LINKS=-lm -lpthread

all: sn_test test

sn_test: sn.c sn_test.c
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS)

# pomyang: pomyang.c sn.c $(SRCDIR)/sieve.c $(SRCDIR)/properSumDiv.c
	# $(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS) -DNDEBUG

# pomyang_debug: pomyang.c sn.c $(SRCDIR)/sieve.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o bin/$@ $^ $(LINKS)

test: sn_test test.sh
	chmod 777 test.sh 
	./test.sh

