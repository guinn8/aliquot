/**
 * @file sieve_sn.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief
 * @date 2021-11-25
 *
 * @copyright Copyright (c) 2021
 *
 */
#include "../additive_sieve/sieve_sn.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
typedef struct enumerate_status_ {
    size_t max;
    size_t sqrt_max;
    size_t blk;
    uint64_t *iter;
    uint64_t *mult;
    enumerated_range_t range;
} enumerate_status_t;

enumerate_status_t *init_enumerate_sn(size_t max, size_t range_len) {
    assert(0 == max % 2);
    assert(0 == max % range_len);
    enumerate_status_t *sts = malloc(sizeof(enumerate_status_t));
    sts->max = max;
    sts->sqrt_max = sqrt(max);
    sts->range.len = range_len;
    sts->blk = range_len;
    sts->range.s = calloc(range_len, sizeof(uint64_t));  // normal array indexing

    // counts number of times multiple has been iterated
    sts->iter = calloc(sts->sqrt_max + 1, sizeof(uint64_t));  // indexing to bound
    for (size_t i = 1; i <= sts->sqrt_max; i++) {
        sts->iter[i] = 1;
    }

    // multiples of i stored at mult[i], not needed but provides constant speedboost
    sts->mult = calloc(sts->sqrt_max + 1, sizeof(uint64_t));  // indexing to bound
    for (size_t i = 1; i <= sts->sqrt_max; i++) {
        sts->mult[i] = i;
    }

    return sts;
}

enumerated_range_t *enumerate_sn(enumerate_handle_t sts) {
    enumerated_range_t *return_range = NULL;
    for (size_t i = 0; i < sts->range.len; i++) {
        sts->range.s[i] = 0;
    }

    if (sts->blk <= sts->max) {
        sts->range.base = (sts->blk - sts->range.len) + 1;  // offset from block index to integer
        uint32_t max_div = sts->sqrt_max;                   // max divisor for any item in block
        while (1 <= max_div) {
            for (size_t i = 1; i <= max_div; i++) {
                if (sts->mult[i] <= sts->blk) {
                    const size_t n = sts->mult[i] - sts->range.base;

                    sts->range.s[n] += i;
                    sts->range.s[n] += (sts->iter[i] > sts->sqrt_max) ? sts->iter[i] : 0;  // accumlate the divisor's pair

                    sts->iter[i]++;
                    sts->mult[i] += i;
                }
            }

            while (sts->mult[max_div] > sts->blk) {  // next max_div is largest multiple in block
                max_div--;
            }
        }

        sts->blk += sts->range.len;
        return_range = &sts->range;
    }

    return return_range;
}

void destroy_enumerate_status(enumerate_handle_t sts) {
    free(sts->mult);
    free(sts->iter);
    free(sts->range.s);
    free(sts);
}
