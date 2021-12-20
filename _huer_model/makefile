SRCDIR=../src
INCS=-I../includes
LINKS=-lm -lpthread

huer_model: huer_model.c $(SRCDIR)/sieve.c  $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)