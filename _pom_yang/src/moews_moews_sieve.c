/**
 * @file sieve_rewrite.c
 * @author Anton Mosunov and Gavin Guinn (gavinguinn1@gmail.com)
 * @brief functions to sieve blocks of sigma
 * @date Originally Apr 18, 2014, modified 2021-12-27
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 * PERFORMANCE SUMMARY:
 * The parameter seg_len has will increase sieve runtime by upto a order-of-magnitude if set to low/high.
 * If seg_len is set to low the sieve will execute far more instructions than needed to complete the same outcome,
 * as it's performance is improves with larger values of seg_len. I believe^ the amount of instructions increases...
 * exponentially as seg_len shrinks, judging by the cachegrind instruction lookup counts. However if seg_len is set...
 * to high the sieve will start thrashing at the cache missing reads/writes in the intermediate sieving steps...
 * thus degrading performance. I believe^^ this effect is more pronounced at bounds lower than ~10^9 as seg_len's that...
 * result in a reasonable amount of operations also happen to require sieving buffers...
 * that fit (very approximately) into my 32mb LL cache.
 *
 * TODO: test ^ and ^^
 *
 * PERFORMANCE:
 * The performance of this sieve is highly variable on the value of seg_len, unfortunately not in a clearly tractable way.
 * It is clear that if seg_len is set to low (whatever that means) the sieve will do a lot of additional operations for the same result.
 * This selected cachegrind/time output from 2 pom_yang run is illustrative:
 *
 * ./cli --bound=100000000 --seg_len=400000 --num_locks=1000000 --preimage_count_bits=1 --num_thread=12
 *  ==*== I   refs:      19,091,262,194
 *  ==*==
 *  ==*== D   refs:       26,107,875,360  (25,263,200,724 rd   + 844,674,636 wr)
 *  ==*== D1  misses:        268,703,487  (   255,239,558 rd   +  13,463,929 wr)
 *  ==*== LLd misses:            406,127  (       274,082 rd   +     132,045 wr)
 *  ==*== D1  miss rate:             1.0% (           1.0%     +         1.6%  )
 * 6.28s user 0.07s system 901% cpu 0.704 total
 *
 * ./cli --bound=100000000 --seg_len=12500 --num_locks=1000000 --preimage_count_bits=1 --num_thread=12
 *  ==*== I   refs:      190,995,976,640
 *  ==*==
 *  ==*== D   refs:       4,448,783,095  (3,658,655,908 rd   + 790,127,187 wr)
 *  ==*== D1  misses:       304,675,262  (  291,274,461 rd   +  13,400,801 wr)
 *  ==*== LLd misses:         8,850,772  (    5,407,998 rd   +   3,442,774 wr)
 *  ==*== D1  miss rate:            6.8% (          8.0%     +         1.7%  )
 * 28.27s user 0.04s system 1052% cpu 2.689 total
 *
 * Note that the seg_len=12500 is running an order of magnitude more instructions and far fewer references on the data cache,
 * these references are relatively inaccurate.
 *
 * The seg_len value can also degrade performace if it is set to high:
 *
 * ./cli --bound=100000000 --seg_len=4000000 --num_locks=$((10**6)) --preimage_count_bits=1 --num_thread=12
 *  ==*== I   refs:      13,701,814,468
 *  ==*==
 *  ==*== D   refs:       3,755,366,091  (2,955,265,779 rd   + 800,100,312 wr)
 *  ==*== D1  misses:       312,794,587  (  298,947,273 rd   +  13,847,314 wr)
 *  ==*== LLd misses:         6,388,217  (    2,678,555 rd   +   3,709,662 wr)
 *  ==*== D1  miss rate:            8.3% (         10.1%     +         1.7%  )
 * 12.94s user 0.36s system 744% cpu 1.786 total
 *
 * Comparing to the seg_len=400000 run we can observe that slightly less instructions are being executed but we are still in the...
 * same order of magnitude. More interestingly we can notice this run made an order of magnitude fewer references to the data cache and..
 * those references were much less accurate, likely accounting for the performace degradation. The vast majority of these misses...
 * originate from interactions with the sigma and numbers buffers.
 *
 * SEE: files in the cg directory for more details
 *
 * ALGORITHM:
 * This is a rewrite of Anton Mosunov's implementation of the [Moews and Moews] sieving algorithm, prepared to utilized...
 * in the tabulation of non-aliquots, sect. 3 of [Chum et al.]. This implementation includes undocumented optimization on...
 * the original method of [Moews and Moews] that are essential for performance, I can't make any sense of why it works however.
 * The _meows_naive_sieve implements the algorithm as described in the paper, it is much slower.
 *
 * CITATIONS:
 *  ->  [Moews and Moews] Moews, D. and Moews, P. C. (1991). A search for aliquot cycles below 10 10.
 *      Mathematics of Computation, 57(196):849.
 *  ->  [Chum et al.] Chum, K., Guy, R. K., Jacobson, J. M. J., and Mosunov, A. S. (2018).
 *      Numerical and statistical analysis of aliquot sequences. Experimental Mathematics, 29(4):414–425.
 */

#include "../inc/moews_moews_sieve.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "../inc/sumdiv.h"

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

static void _moews_sieve_odd(sieve_worker_t *worker, uint64_t seg_start, bool squared);
static void _prime_sieve(uint32_t max_prime, uint32_t *sieved_primes);
static void _meows_naive_sieve(uint64_t N, uint64_t M, uint64_t *sigma, const uint32_t *primes);

/**
 * @brief runs sieve for odd sigma(m) or odd sigma(m * m)
 *
 * @param worker see struct definiton
 * @param seg_start start of range to sieve
 * @param squared sigma(m) or odd sigma(m * m)
 */
static void _moews_sieve_odd(sieve_worker_t *worker, uint64_t seg_start, const bool squared) {
    assert(EVEN(worker->cfg->seg_len));
    assert(EVEN(seg_start));
    worker->seg_start = seg_start;
    worker->squared = squared;
    const uint64_t max_prime = (uint64_t)sqrt(seg_start + worker->cfg->seg_len);

    for (size_t i = 0; i < worker->cfg->buf_len; i++) {
        worker->sigma_buf[i] = 1;
        worker->q_buf[i] = seg_start + 1 + (2 * i);
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
                        worker->q_buf[(j - 1) / 2] /= p_pow;
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
        if (worker->q_buf[i] > 1) {
            if (squared) {
                worker->sigma_buf[i] *= (worker->q_buf[i] * (worker->q_buf[i] + 1L) + 1L);
            } else {
                worker->sigma_buf[i] *= (worker->q_buf[i] + 1);
            }
        }
    }
}

/* See header for documentation */
sieve_config_t *moews_init_sieve(size_t bound, size_t seg_len) {
    sieve_config_t *cfg = malloc(sizeof(sieve_config_t));

    const uint64_t max_prime = (uint64_t)sqrt(2 * bound);
    cfg->primes_len = (1.25506 * (max_prime + 1) / log(max_prime + 1)) + 1;
    cfg->primes = (uint32_t *)calloc(sizeof(uint32_t), cfg->primes_len);
    _prime_sieve(max_prime, cfg->primes);

    cfg->bound = bound;
    cfg->seg_len = seg_len;
    cfg->buf_len = seg_len / 2;
    assert(EVEN(cfg->buf_len));

    return cfg;
}

/* See header for documentation */
void moews_destroy_sieve(sieve_config_t *cfg) {
    free(cfg->primes);
    free(cfg);
}

/* See header for documentation */
sieve_worker_t *moews_init_worker(const sieve_config_t *cfg) {
    sieve_worker_t *worker = malloc(sizeof(sieve_worker_t));
    worker->cfg = cfg;
    worker->sigma_buf = calloc(cfg->buf_len, sizeof(uint64_t));
    worker->q_buf = calloc(cfg->buf_len, sizeof(uint64_t));
    worker->is_prime_buf = calloc(cfg->buf_len, sizeof(bool));
    worker->last_sieve_standard = -1;
    return worker;
}

/* See header for documentation */
void destroy_worker(sieve_worker_t *worker) {
    free(worker->sigma_buf);
    free(worker->q_buf);
    free(worker->is_prime_buf);
    free(worker);
}

/* See header for documentation */
void moews_sieve_odd_standard(sieve_worker_t *worker, const uint64_t seg_start) {
    _moews_sieve_odd(worker, seg_start, 0);
    worker->last_sieve_standard = seg_start;

    // use sigma data to store prime information
    for (size_t i = 0; i < worker->cfg->buf_len; i++) {
        uint64_t m = SEG_OFFSET(worker->seg_start, i);
        uint64_t sigma_m = worker->sigma_buf[i];
        worker->is_prime_buf[i] = IS_M_PRIME(m, sigma_m);
    }
}

/* See header for documentation */
void moews_sieve_odd_squared(sieve_worker_t *worker, const uint64_t seg_start) {
    _moews_sieve_odd(worker, seg_start, 1);
}

/* See header for documentation */
uint64_t moews_lookup_sigma_m(const sieve_worker_t *worker, uint64_t m) {
    DEBUG_ASSERT(assert(worker->squared == false));
    DEBUG_ASSERT(assert(ODD(m)));
    size_t index = SEG_DEOFFSET(worker->seg_start, m);
    DEBUG_ASSERT(assert(index < worker->cfg->buf_len));

    uint64_t sigma_m = worker->sigma_buf[index];
    DEBUG_ASSERT(assert(sigma_m == sumdiv_sigma(m)));
    return sigma_m;
}

/* See header for documentation */
uint64_t moews_lookup_sigma_m_squared(const sieve_worker_t *worker, uint64_t m) {
    DEBUG_ASSERT(assert(worker->squared == true));
    DEBUG_ASSERT(assert(ODD(m)));
    size_t index = SEG_DEOFFSET(worker->seg_start, m);
    DEBUG_ASSERT(assert(index < worker->cfg->buf_len));

    uint64_t sigma_m = worker->sigma_buf[index];
    DEBUG_ASSERT(assert(sigma_m == sumdiv_sigma(m * m)));
    return sigma_m;
}

/* See header for documentation */
bool moews_check_prime(const sieve_worker_t *worker, uint64_t m) {
    DEBUG_ASSERT(assert(-1 != worker->last_sieve_standard));
    DEBUG_ASSERT(assert(m > (uint64_t)worker->last_sieve_standard));
    DEBUG_ASSERT(assert(m < (uint64_t)worker->last_sieve_standard + worker->cfg->seg_len));

    size_t index = SEG_DEOFFSET(worker->seg_start, m);
    return worker->is_prime_buf[index];
}

/* See header for documentation */
size_t moews_estimate_heap_usage(size_t bound, size_t seg_len, size_t num_workers) {
    size_t total = 0;

    const uint64_t max_prime = (uint64_t)sqrt(2 * bound);
    size_t primes_len = (1.25506 * (max_prime + 1) / log(max_prime + 1)) + 1;
    total += primes_len * sizeof(uint32_t);

    size_t buf_len = seg_len / 2;
    total += buf_len * sizeof(uint64_t) * num_workers;  // sigma_buf
    total += buf_len * sizeof(uint64_t) * num_workers;  // q_buf
    printf("Sieving buffers take %.2fmb\n", (float)total / 1048576);
    total += buf_len * sizeof(bool) * num_workers;  // is_prime_buf

    return total;
}

/**
 * @brief prime sieve used to initialize the Moews and Moews sieve
 *
 * @param max_prime max prime to sieve
 * @param sieved_primes pre-alloc'ed buffer to fill with primes
 */
void _prime_sieve(uint32_t max_prime, uint32_t *sieved_primes) {
    bool *sieve_is_prime = calloc(max_prime + 1, sizeof(bool));
    for (size_t i = 0; i < max_prime + 1; i++) {
        sieve_is_prime[i] = true;
    }
    sieved_primes[0] = 0;

    for (size_t i = 2; i <= max_prime; i++) {
        if (sieve_is_prime[i]) {
            sieved_primes[0]++;
            sieved_primes[sieved_primes[0]] = i;
            for (size_t j = 2 * i; j <= max_prime; j += i) {
                sieve_is_prime[j] = false;
            }
        }
    }
    free(sieve_is_prime);
}

/**
 * Summary of sieving algorithm from [Meows and Meows]
 *
 * Suppose that we wish to compute sigma(M), ... , sigma(M + N - X).
 * The elements q[j] are initialized to M + j and r[j] to 1, for j = 0, ... , N - 1.
 *
 * Then for each prime power p^e with p <= sqrt(M + N - 1) and p^e <= M + N - X,
 * and for each j with M + j divisible by p^e and not by p^(e+1),
 * the r[j]-values are multiplied by sigma(p^e) and the q[j]-values are divided by p^e.
 *
 * At the end of this process the q[j]-values will be 1 unless
 * M + j had a prime factor p exceeding sqrt(M + N - 1), in which
 * case q[j] will be p. Hence, if we make a final pass multiplying each r[j] by
 * q[j] + 1 if q[j] != 1, we will have set each r[j] equal to sigma(M+j).
 *
 * @brief [Meows and Meows] sieving method to enumerate range of sigma
 *
 * @param N length of range to enumerate
 * @param M base number of range to enumerate
 * @param sigma buffer to be filled with enumeration
 * @param primes buffer containing all primes <= M + N
 *
 * NOTE:    this is included primarily for documentary purposes, it's much slower and untested compared to the...
 *          production sieve
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"  // included for documentary purposes
void _meows_naive_sieve(const uint64_t N, const uint64_t M, uint64_t *sigma, const uint32_t *primes) {
    uint64_t *q = (uint64_t *)calloc(N, sizeof(uint64_t));
    for (size_t j = 0; j < N; j++) {
        sigma[j] = 1;
        q[j] = M + j;
    }

    size_t p_ind = 1;
    uint64_t p;
    const uint64_t max_prime = sqrt(M + N);
    while ((p = primes[p_ind++]) <= max_prime) {  // for each prime power p^e with p <= sqrt(M + N - 1)
        uint64_t p_pow = p;
        while (p_pow <= (M + N)) {
            const uint64_t next_p_pow = p_pow * p;
            for (size_t j = 0; j < N; j++) {
                if ((0 == (M + j) % p_pow) && (0 != (M + j) % next_p_pow)) {
                    q[j] /= p_pow;
                    sigma[j] *= ((next_p_pow - 1) / (p - 1));
                }
            }
            p_pow = next_p_pow;
        }
    }

    for (size_t i = ((M == 0) ? 1 : 0); i < N; i++) {
        if (q[i] > 1) {
            sigma[i] *= (q[i] + 1);
        }
        // printf("sigma(%lu)=%lu\n", i + l , sigma[i]);
        printf("%lu\n", sigma[i]);
    }

    free(q);
}
#pragma GCC diagnostic pop
