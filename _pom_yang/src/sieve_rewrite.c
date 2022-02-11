/**
 * @file sieve_rewrite.c
 * @author Anton Mosunov and Gavin Guinn (gavinguinn1@gmail.com)
 * @brief functions to sieve blocks of sigma
 * @date Originally Apr 18, 2014, modified 2021-12-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "../inc/sieve_rewrite.h"

#include <assert.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EVEN(x) (0 == (x) % 2)

// #define DEBUG_SIEVE
#ifdef DEBUG_SIEVE
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif

static uint32_t *primes;
static size_t num_primes;

uint64_t init_sigma_sieve(const size_t bound) {
    const uint64_t max_prime = (uint64_t)sqrt(2 * bound);
    num_primes = (1.25506 * (max_prime + 1) / log(max_prime + 1)) + 1;
    primes = (uint32_t *)calloc(sizeof(uint32_t), num_primes);
    prime_sieve(max_prime, primes);
    return sizeof(uint32_t) * num_primes;
}

void destroy_sigma_sieve(void) {
    free(primes);
}

void prime_sieve(const uint32_t max_prime, uint32_t *returned_primes) {
    bool *is_prime = calloc(max_prime + 1, sizeof(bool));
    for (size_t i = 0; i < max_prime + 1; i++) {
        is_prime[i] = true;
    }
    returned_primes[0] = 0;

    for (size_t i = 2; i <= max_prime; i++) {
        if (is_prime[i]) {
            returned_primes[0]++;
            returned_primes[returned_primes[0]] = i;
            for (size_t j = 2 * i; j <= max_prime; j += i) {
                is_prime[j] = false;
            }
        }
    }
    free(is_prime);
}

void sigma_sieve_odd(const uint64_t seg_len, const uint64_t seg_start, uint64_t *sigma_buf, uint64_t *numbers, const bool squared) {
    assert(EVEN(seg_len));
    assert(EVEN(seg_start));

    const uint64_t max_prime = (uint64_t)sqrt(seg_start + seg_len);

    for (size_t i = 0; i < seg_len / 2; i++) {
        sigma_buf[i] = 1;
        numbers[i] = seg_start + 1 + (2 * i);
    }

    uint64_t offset[100];
    size_t prime_ind = 2;
    uint64_t p;
    while (max_prime >= (p = primes[prime_ind++])) {
        uint64_t k = 0;
        offset[k] = (p - (seg_start % p)) % p;

        if (EVEN(offset[k])) {
            offset[k] += p;
        }

        for (uint64_t p_pow = p; p_pow <= (seg_start + seg_len); p_pow *= p) {
            if (offset[k] > seg_len) {
                break;
            }

            uint64_t step = p_pow * p;
            offset[++k] = (step - (seg_start % step)) % step;
            if (EVEN(offset[k])) {
                offset[k] += step;
            }

            uint64_t h = offset[k - 1];
            for (uint64_t s = 0; s < p; s++) {
                if (h != offset[k]) {
                    for (uint64_t j = h; j < seg_len; j += 2 * step) {
                        numbers[(j - 1) / 2] /= p_pow;
                        if (squared) {
                            sigma_buf[(j - 1) / 2] *= ((step * p_pow - 1) / (p - 1));
                        } else {
                            sigma_buf[(j - 1) / 2] *= ((step - 1) / (p - 1));
                        }
                    }
                }
                h += (2 * p_pow);
            }
        }
    }

    for (size_t i = 0; i < seg_len / 2; i++) {
        if (numbers[i] > 1) {
            if (squared) {
                sigma_buf[i] *= (numbers[i] * (numbers[i] + 1L) + 1L);
            } else {
                sigma_buf[i] *= (numbers[i] + 1);
            }
        }
        DEBUG(printf("sigma(%lu)=%lu\n", seg_start + 1 + (2 * i), (uint64_t)sigma_buf[i]);)
    }
}
