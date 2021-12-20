SRCDIR=../src
INCS=-I../includes
LINKS=-lm -lpthread

factor_ennum: factor_ennum.c $(SRCDIR)/properSumDiv.c
	$(CC) $(CFLAGS) $(INCS) -o $@ $^ $(LINKS)