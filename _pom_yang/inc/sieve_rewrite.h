/**
 * @file sieve_rewrite.c
 * @author Anton Mosunov and Gavin Guinn (gavinguinn1@gmail.com)
 * @brief functions to sieve blocks of sigma
 * @date Originally Apr 18, 2014, modified 2021-12-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef _POM_YANG_INC_SIEVE_REWRITE_H_
#define _POM_YANG_INC_SIEVE_REWRITE_H_

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    size_t bound;
    size_t seg_len;
    size_t sigma_buf_len;
    size_t num_primes;
    uint32_t *primes;
} sieve_config_t;

typedef struct {
    sieve_config_t *cfg;
    bool squared;
    size_t seg_start;
    uint64_t *sigma_buf;
    uint64_t *numbers_buf;
    bool *is_prime;
} sieve_worker_t;

sieve_config_t *init_sigma_sieve(const size_t bound, const size_t seg_len);
sieve_worker_t *init_sieve_worker(sieve_config_t *cfg);
void destroy_sieve(sieve_config_t *cfg);
void prime_sieve(const uint32_t max_prime, uint32_t *primes);
void sigma_sieve_odd(sieve_worker_t *worker, const uint64_t seg_start, const bool squared);
uint64_t get_sigma_m(sieve_worker_t *worker, uint64_t m);
void get_primes(sieve_worker_t *worker);
bool is_prime(sieve_worker_t *worker, uint64_t m);
void destroy_worker(sieve_worker_t *worker);

#endif  // _POM_YANG_INC_SIEVE_REWRITE_H_
