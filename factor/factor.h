/**
 * @file factor.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-03-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdint.h>
#include <stdlib.h>
#include "gmp.h"

/**
 * @brief Holds a prime factorization 
 * @struct factors
 * 
 * @var factors::p
 * List of input's prime factors
 * 
 * @var factors::e 
 * The exponent of each respective prime factor
 * 
 * @var factors::nfactors
 * Size of the factor lists
 */
struct factors {
    mpz_t *p;
    uint64_t *e;
    size_t nfactors;
};

uint64_t factor_sigma(uint64_t n);
uint64_t factor_s(uint64_t n);
void factor(mpz_t, struct factors *);

void factor_clear(struct factors *factors);
