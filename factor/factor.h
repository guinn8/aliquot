#include <stdint.h>

#include "gmp.h"

struct factors {
    mpz_t *p;
    unsigned long *e;
    long nfactors;
};

uint64_t factor_sigma(uint64_t n);
uint64_t factor_s(uint64_t n);
void factor(mpz_t, struct factors *);

void factor_clear(struct factors *factors);
