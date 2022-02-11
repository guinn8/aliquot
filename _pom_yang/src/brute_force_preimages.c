/**
 * @file brute_force_preimages.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief
 * @date 2022-02-11
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "../inc/brute_force_preimages.h"

#include <assert.h>
#include <omp.h>

#include "../inc/properSumDiv.h"

#define EVEN(x) (0 == (x) % 2)
#define F_OFFSET(x) ((x / 2) - 1)
#define F_DE_OFFSET(x) ((x + 1) * 2)

uint8_t *brute_force_preimages(size_t bound) {
    assert(EVEN(bound));

    const size_t f_len = bound / 2;
    uint8_t *f = calloc(f_len, sizeof(uint8_t));

#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 2; i <= 2 * bound; i += 2) {
            uint64_t s_n = wheelDivSum(i);
            if (s_n <= bound && EVEN(s_n)) {
#pragma omp atomic
                f[F_OFFSET(s_n)]++;  // TODO(gavin): check overflow
            }
        }

#pragma omp for
        for (size_t i = 1; i <= bound; i += 2) {
            uint64_t s_n = wheelDivSum(i * i);
            if (s_n <= bound && EVEN(s_n)) {
#pragma omp atomic
                f[F_OFFSET(s_n)]++;  // TODO(gavin): check overflow
            }
        }
    }
    return f;
}

uint64_t *brute_force_preimage_counts(size_t bound) {
    const size_t f_len = bound / 2;
    uint8_t *f = brute_force_preimages(bound);
    uint64_t *count = calloc(UINT8_MAX, sizeof(uint64_t));
    for (size_t i = 0; i < f_len; i++) {
        count[f[i]]++;
    }

    // for (size_t i = 0; i < 8; i++) {
    //     printf("%ld: %ld\n", i, count[i]);
    // }

    return count;
}