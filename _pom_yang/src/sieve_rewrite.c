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

#include "../inc/properSumDiv.h"

#define EVEN(x) (0 == (x) % 2)
#define ODD(x) (0 == ((x) + 1) % 2)
#define SEG_OFFSET(seg_start, i) ((seg_start) + (2 * (i)) + 1)
#define SEG_DEOFFSET(seg_start, m) ((m - 1 - seg_start) / 2)
#define IS_M_PRIME(m, sigma_m) ((m) + 1 == (sigma_m))

#ifdef DEBUG_ASSERT_ON
#define DEBUG_ASSERT(x) x
#else
#define DEBUG_ASSERT(x)
#endif

void _sigma_sieve_odd(sieve_worker_t *worker, const uint64_t seg_start, const bool squared);

void sigma_sieve_odd(sieve_worker_t *worker, const uint64_t seg_start) {
    _sigma_sieve_odd(worker, seg_start, 0);
}

void sigma_sieve_odd_squared(sieve_worker_t *worker, const uint64_t seg_start) {
    _sigma_sieve_odd(worker, seg_start, 1);
}

size_t sieve_estimate_heap_usage(const size_t bound, const size_t seg_len, size_t num_workers) {
    size_t total = 0;

    const uint64_t max_prime = (uint64_t)sqrt(2 * bound);
    size_t num_primes = (1.25506 * (max_prime + 1) / log(max_prime + 1)) + 1;
    total += num_primes * sizeof(uint32_t);

    size_t sigma_buf_len = seg_len / 2;
    total += sigma_buf_len * sizeof(uint64_t) * num_workers;  // sigma_buf
    total += sigma_buf_len * sizeof(uint64_t) * num_workers;  // numbers_buf
    total += sigma_buf_len * sizeof(bool) * num_workers;      // is_prime

    return total;
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

sieve_config_t *init_sigma_sieve(const size_t bound, const size_t seg_len) {
    sieve_config_t *cfg = malloc(sizeof(sieve_config_t));

    const uint64_t max_prime = (uint64_t)sqrt(2 * bound);
    cfg->num_primes = (1.25506 * (max_prime + 1) / log(max_prime + 1)) + 1;
    cfg->primes = (uint32_t *)calloc(sizeof(uint32_t), cfg->num_primes);
    prime_sieve(max_prime, cfg->primes);

    cfg->bound = bound;
    cfg->seg_len = seg_len;
    cfg->sigma_buf_len = seg_len / 2;
    assert(EVEN(cfg->sigma_buf_len));

    return cfg;
}

void destroy_sieve(sieve_config_t *cfg) {
    free(cfg->primes);
    free(cfg);
}

void destroy_worker(sieve_worker_t *worker) {
    free(worker->sigma_buf);
    free(worker->numbers_buf);
    free(worker->is_prime);
    free(worker);
}

sieve_worker_t *init_sieve_worker(sieve_config_t *cfg) {
    sieve_worker_t *worker = malloc(sizeof(sieve_worker_t));
    worker->cfg = cfg;
    worker->sigma_buf = calloc(cfg->sigma_buf_len, sizeof(uint64_t));    // todo: move into sieve
    worker->numbers_buf = calloc(cfg->sigma_buf_len, sizeof(uint64_t));  // todo: move into sieve
    worker->is_prime = calloc(cfg->sigma_buf_len, sizeof(bool));         // todo: move into sieve
    return worker;
}

void _sigma_sieve_odd(sieve_worker_t *worker, const uint64_t seg_start, const bool squared) {
    assert(EVEN(worker->cfg->seg_len));
    assert(EVEN(seg_start));
    worker->seg_start = seg_start;
    worker->squared = squared;
    const uint64_t max_prime = (uint64_t)sqrt(seg_start + worker->cfg->seg_len);

    for (size_t i = 0; i < worker->cfg->seg_len / 2; i++) {
        worker->sigma_buf[i] = 1;
        worker->numbers_buf[i] = seg_start + 1 + (2 * i);
    }

    uint64_t offset[100];
    size_t prime_ind = 2;
    uint64_t p;
    while (max_prime >= (p = worker->cfg->primes[prime_ind++])) {
        uint64_t k = 0;
        offset[k] = (p - (seg_start % p)) % p;

        if (EVEN(offset[k])) {
            offset[k] += p;
        }

        for (uint64_t p_pow = p; p_pow <= (seg_start + worker->cfg->seg_len); p_pow *= p) {
            if (offset[k] > worker->cfg->seg_len) {
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
                    for (uint64_t j = h; j < worker->cfg->seg_len; j += 2 * step) {
                        worker->numbers_buf[(j - 1) / 2] /= p_pow;
                        if (squared) {
                            worker->sigma_buf[(j - 1) / 2] *= ((step * p_pow - 1) / (p - 1));
                        } else {
                            worker->sigma_buf[(j - 1) / 2] *= ((step - 1) / (p - 1));
                        }
                    }
                }
                h += (2 * p_pow);
            }
        }
    }

    for (size_t i = 0; i < worker->cfg->seg_len / 2; i++) {
        if (worker->numbers_buf[i] > 1) {
            if (squared) {
                worker->sigma_buf[i] *= (worker->numbers_buf[i] * (worker->numbers_buf[i] + 1L) + 1L);
            } else {
                worker->sigma_buf[i] *= (worker->numbers_buf[i] + 1);
            }
        }
    }
}

uint64_t get_sigma_m(sieve_worker_t *worker, uint64_t m) {
    DEBUG_ASSERT(assert(ODD(m)));
    size_t index = SEG_DEOFFSET(worker->seg_start, m);
    DEBUG_ASSERT(assert(index < worker->cfg->sigma_buf_len));

    uint64_t sigma_m = worker->sigma_buf[index];
    DEBUG_ASSERT(assert(sigma_m == wheelDivSigma((worker->squared) ? m * m : m)));

    return sigma_m;
}

void get_primes(sieve_worker_t *worker) {
    assert(worker->squared == false);
    for (size_t i = 0; i < worker->cfg->sigma_buf_len; i++) {
        uint64_t m = SEG_OFFSET(worker->seg_start, i);
        uint64_t sigma_m = worker->sigma_buf[i];
        worker->is_prime[i] = IS_M_PRIME(m, sigma_m);
    }
}

bool is_prime(sieve_worker_t *worker, uint64_t m) {
    size_t index = SEG_DEOFFSET(worker->seg_start, m);
    return worker->is_prime[index];
}
