/*
 * sieve.c
 *
 *  Created on: Apr 18, 2014
 *      Author: antonmosunov
 *  USED WITH PERMISSION
 */

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>

#include "sieve.h"

void prime_sieve(const uint max_prime, uint * primes)
{
    uint i, j;

    primes[0] = 0;

    char is_prime[max_prime + 1];
    memset(is_prime, 1, max_prime + 1);

    for (i = 2; i <= max_prime; i++)
    {
        if (is_prime[i])
        {
            primes[++primes[0]] = i;

            for (j = (i << 1); j <= max_prime; j += i)
            {
                is_prime[j] = 0;
            }
        }
    }

    primes[primes[0] + 1] = 0x7FFFFFFF;
}


void regular_sieve(const uint max_prime, const uint64_t blocksize, uint ** factors, const uint * primes, const int flags)
{
    uint64_t i, j;

    uint p, * f;

    factors[0][0] = 1;
    factors[0][1] = 0;
    factors[1][0] = 1;
    factors[1][1] = 1;

    for (i = 0; i < blocksize; i++)
    {
        factors[i][0] = 0;
    }

    for (i = 1, p = primes[i]; p < max_prime; p = primes[++i])
    {
        for (j = p; j < blocksize; j += p)
        {
            f = factors[j];
            f[++f[0]] = (flags & WITH_INDICES) ? i : p;
        }
    }
}


void segmented_sieve(uint max_prime, uint64_t blocksize, uint64_t l, uint ** factors, const uint * primes, const int flags)
{
    if (l == 0)
    {
        regular_sieve(max_prime, blocksize, factors, primes, flags);
        return;
    }

    uint i, j, k, p, offset, * f;

    for (k = 0; k < blocksize; k++)
    {
        factors[k][0] = 0;
    }

    for (i = 1, p = primes[i]; p < max_prime; p = primes[++i])
    {
        offset = (p - (l % p)) % p;

        for (j = offset; j < blocksize; j += p)
        {
            f = factors[j];
            f[++f[0]] = (flags & WITH_INDICES) ? i : p;
        }
    }
}

/**
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
 */

/**
 * @brief [Meows and Meows 91] sieving method to enumerate range of sigma 
 * 
 * @param N length of range to enumerate
 * @param M base number of range to enumerate
 * @param sigma buffer to be filled with enumeration
 * @param primes buffer containing all primes <= M + N
 */
void sum_of_divisors(const uint64_t N, const uint64_t M, uint64_t * sigma, const uint32_t * primes) {
    uint64_t *q = (uint64_t *)calloc(N, sizeof(uint64_t));

    // q[j] are init'ed to M + j and r[j] to 1
    for (size_t j = 0; j < N; j++) {
        sigma[j] = 1;
        q[j] = M + j;
    }

    size_t p_ind = 1;
    uint64_t p;
    // printf("********M = %ld********\n", M);

    const uint64_t max_prime = sqrt(M + N);
    while ((p = primes[p_ind++]) <= max_prime) {  // for each prime power p^e with p <= sqrt(M + N - 1)
        // printf("--------p = %ld--------\n", p);
        uint64_t p_pow = p;
        // uint64_t offset = (p_pow - (M % p_pow)) % p_pow;  // todo(Gavin): nicer as a do-while
        while (p_pow <= (M + N)) { // while (p_pow <= (M + N) && offset <= N) {
            // assert(0 == (M + offset) % p_pow);
            const uint64_t next_p_pow = p_pow * p;
            // printf("offset = %ld, p_pow=%ld, step=%ld\n", offset + M, p_pow, next_p_pow);


            #define NAIVE
            #ifdef NAIVE
            for (size_t j = 0; j < N; j++) {
                if ((0 == (M + j) % p_pow) && (0 != (M + j) % next_p_pow)) {
                    //printf("p_pow=%ld, next_p_pow=%ld, next_mult=%ld\n", p_pow, next_p_pow, M + j);
                    q[j] /= p_pow;
                    sigma[j] *= ((next_p_pow - 1) / (p - 1));
                }
            }
            #else
            uint64_t next_mult = offset + M;
            while (next_mult <= N + M) {
                if (0 != next_mult % next_p_pow) {
                    for (size_t j = next_mult; j < N + M; j += 2 * next_p_pow) {
                        if (false == hit[j - M]) {
                            assert(0 != j % next_p_pow);
                            assert(0 == (next_mult) % p_pow);
                            // printf("p_pow=%ld, next_p_pow=%ld, next_mult=%ld, j=%ld\n", p_pow, next_p_pow, next_mult, j);
                            q[j - M] /= p_pow;
                            sigma[j - M] *= ((next_p_pow - 1) / (p - 1));
                            hit[j - M] = true;
                        }
                    }
                }
                next_mult += p_pow;
            }
            #endif


            p_pow = next_p_pow;
            // offset = (p_pow - (M % p_pow)) % p_pow;
            // printf("\n"); p 
        }
    }

    for (size_t i = ((M == 0) ? 1 : 0) ; i < N; i++) {
        if (q[i] > 1) {
            sigma[i] *= (q[i] + 1);
        }
        // printf("sigma(%lu)=%lu\n", i + l , sigma[i]);
        printf("%lu\n", sigma[i]);
    }

    free(q);
}

#if BITSIZE == 32
void sum_of_divisors_odd(const uint64_t blocksize, const uint64_t l, uint * sigma, const uint * primes/*, uint64_t * max_val, uint64_t * max_sigma*/)
#else
void sum_of_divisors_odd(const uint64_t blocksize, const uint64_t l, uint64_t * sigma, const uint * primes/*, uint64_t * max_val, uint64_t * max_sigma*/)
#endif
{
    if (((blocksize & 1) == 1) || ((l & 1) == 1))
    {
        perror("The parameters blocksize and block_base must be even.\n");
        exit(1);
    }

    const uint64_t max_prime = (uint64_t) sqrt(l + blocksize);

    uint64_t h, i, j, k, p, s, p_pow, step;

    uint64_t * numbers = (uint64_t *) malloc(sizeof(uint64_t) * (blocksize >> 1));

    for (i = 0, j = l + 1; i < (blocksize >> 1); i++, j += 2)
    {
        sigma[i] = 1;
        numbers[i] = j;
    }

    if (l == 0)
    {
        sigma[0] = 1;
    }

    uint64_t offset[100];

    for (i = 2, p = primes[i]; p <= max_prime; p = primes[++i])
    {
        k = 0;
        offset[k] = (p - (l % p)) % p;

        if ((offset[k] & 1) == 0)
        {
            offset[k] += p;
        }

        for (p_pow = p, step = p * p; p_pow <= (l + blocksize); p_pow *= p, step *= p)
        {
            if (offset[k] > blocksize)
            {
                break;
            }

            // #ifdef DEBUG
            // printf("p_pow=%lu, offset=%lu, step=%lu\t", p_pow, offset[k], step);
            // #endif

            offset[++k] = (step - (l % step)) % step;
            assert(k < 100);
            if ((offset[k] & 1) == 0)
            {
                offset[k] += step;
            }

            for (s = 0, h = offset[k-1]; s < p; s++, h += (p_pow << 1))
            {
                // printf("offset[k-1] = %ld, offset[k] = %ld\n", offset[k-1], offset[k]);
                if (h != offset[k])
                {
                    for (j = h; j < blocksize; j += (step << 1))
                    {
                        #ifdef DEBUG
                        printf(", %lu (%lu)", l + j, (j-1)>>1);
                        #endif

                        numbers[(j-1) >> 1] /= p_pow;
                        sigma[(j-1) >> 1] *= ((step - 1) / (p - 1));
                    }
                }
            }

            #ifdef DEBUG
            printf("\n");
            #endif
        }
    }

    i = 0;
    j = l + 1;

    for (; i < (blocksize >> 1); i++, j += 2)
    {
        if (numbers[i] > 1)
        {
            sigma[i] *= (numbers[i] + 1);
        }
        
        #ifdef DEBUG
        printf("sigma(%lu)=%lu\n", j, (uint64_t) sigma[i]);
        #endif

        #ifdef WITH_PARI
        GEN g = stoi(j);
        GEN v = gsumdivk(g, 1);

        if (sigma[i] != gtolong(v))
        d ..cd
            char err[240];
            sprintf(err, "Error, sigma(%lu)=%lu, but we computed %lu.\n", j, gtolong(v), (uint64_t) sigma[i]);
            perror(err);
            cgiv(v);
            cgiv(g);
            pari_close();
            exit(1);
        }

        cgiv(v);
        cgiv(g);
        #endif
    }

    free(numbers);
}

#if BITSIZE == 32
void sum_of_divisors_odd2(const uint64_t blocksize, const uint64_t l, uint * sigma, const uint * primes)
#else
void sum_of_divisors_odd2(const uint64_t blocksize, const uint64_t l, uint64_t * sigma, const uint * primes)
#endif
{
    if (((blocksize & 1) == 1) || ((l & 1) == 1))
    {
        perror("The parameters blocksize and block_base must be even.\n");
        exit(1);
    }

    const uint64_t max_prime = (uint64_t) sqrt(l + blocksize);

    uint64_t h, i, j, k, p, s, p_pow, step;

    uint64_t * numbers = (uint64_t *) malloc(sizeof(uint64_t) * (blocksize >> 1));

    for (i = 0, j = l + 1; i < (blocksize >> 1); i++, j += 2)
    {
        sigma[i] = 1;
        numbers[i] = j;
    }

    if (l == 0)
    {
        sigma[0] = 1;
    }

    uint64_t offset[100];

    for (i = 2, p = primes[i]; p <= max_prime; p = primes[++i])
    {
        k = 0;
        offset[k] = (p - (l % p)) % p;

        if ((offset[k] & 1) == 0)
        {
            offset[k] += p;
        }

        for (p_pow = p, step = p * p; p_pow <= (l + blocksize); p_pow *= p, step *= p)
        {
            if (offset[k] > blocksize)
            {
                break;
            }

            #ifdef DEBUG
            printf("p_pow=%lu, offset=%lu, step=%lu\t", p_pow, offset[k], step);
            #endif

            offset[++k] = (step - (l % step)) % step;

            if ((offset[k] & 1) == 0)
            {
                offset[k] += step;
            }

            for (s = 0, h = offset[k-1]; s < p; s++, h += (p_pow << 1))
            {
                if (h != offset[k])
                {
                    for (j = h; j < blocksize; j += (step << 1))
                    {
                        #ifdef DEBUG
                        printf(", %lu (%lu)", l + j, (j-1)>>1);
                        #endif

                        numbers[(j-1) >> 1] /= p_pow;
                        sigma[(j-1) >> 1] *= ((step*p_pow - 1) / (p - 1));
                    }
                }
            }

            #ifdef DEBUG
            printf("\n");
            #endif
        }
    }

    uint64_t j_squared;

    for (i = 0, j = l + 1, j_squared = j * j; i < (blocksize >> 1); i++, j_squared += ((j + 1) << 2), j += 2)
    {
        if (numbers[i] > 1)
        {
            sigma[i] *= (numbers[i] * (numbers[i] + 1L) + 1L);
        }

        s = sigma[i] - j_squared;

        #ifdef DEBUG
        printf("sigma(%lu^2=%lu)=%lu, s(%lu)=%lu\n", j, j_squared, (uint64_t) sigma[i], j*j, s);
        #endif

        #ifdef WITH_PARI
        GEN g = stoi(j_squared);
        GEN v = gsumdivk(g, 1);

        if (sigma[i] != gtolong(v))
        {
            char err[240];
            sprintf(err, "Error, sigma(%lu^2=%lu)=%lu, but we computed %lu.\n", j, j_squared, gtolong(v), (uint64_t) sigma[i]);
            perror(err);
            cgiv(v);
            cgiv(g);
            pari_close();
            exit(1);
        }

        cgiv(v);
        cgiv(g);
        #endif
    }

    free(numbers);
}
