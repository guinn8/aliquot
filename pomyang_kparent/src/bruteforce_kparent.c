/**
 * @file bruteforce_kparent.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief   Uses a simple bruteforce method to compute number of preimages for even number upto bound.
 *          Useful as a backstop test to the pom_yang algorithm
 * @date 2022-02-11
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#include "../inc/bruteforce_kparent.h"

#include <assert.h>
#include <omp.h>

#include "../inc/sumdiv_util.h"

#define EVEN(x) (0 == (x) % 2)
#define F_OFFSET(x) ((x / 2) - 1)
#define F_DE_OFFSET(x) ((x + 1) * 2)

/* See header for documentation */
uint8_t *bf_kparent(size_t bound) {
    assert(EVEN(bound));

    const size_t f_len = bound / 2;
    uint8_t *f = calloc(f_len, sizeof(uint8_t));

#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 2; i <= 2 * bound; i += 2) {
            uint64_t s_n = sumdiv_s(i);
            if (s_n <= bound && EVEN(s_n)) {
#pragma omp atomic
                f[F_OFFSET(s_n)]++;
            }
        }

#pragma omp for
        for (size_t i = 1; i <= bound; i += 2) {
            uint64_t s_n = sumdiv_s(i * i);
            if (s_n <= bound && EVEN(s_n)) {
#pragma omp atomic
                f[F_OFFSET(s_n)]++;
            }
        }
    }
    return f;
}

/* See header for documentation */
uint64_t *bf_kparent_counts(size_t bound) {
    const size_t f_len = bound / 2;
    uint8_t *f = bf_kparent(bound);
    uint64_t *count = calloc(UINT8_MAX, sizeof(uint64_t));

    for (size_t i = 0; i < f_len; i++) {
        count[f[i]]++;
    }

    return count;
}