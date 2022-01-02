/**
 * @file sieve_rewrite.c
 * @author Anton Mosunov and Gavin Guinn (gavinguinn1@gmail.com)
 * @brief functions to sieve blocks of sigma
 * @date Originally Apr 18, 2014, modified 2021-12-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef INC_SIEVE_REWRITE_H_
#define INC_SIEVE_REWRITE_H_

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

void init_sigma_sieve(const size_t bound);
void destroy_sigma_sieve(void);
void prime_sieve(const uint32_t max_prime, uint32_t *primes);
void sigma_sieve_odd(const uint64_t blocksize, const uint64_t l, uint64_t *sigma, const bool squared);

#endif  // INC_SIEVE_REWRITE_H_
