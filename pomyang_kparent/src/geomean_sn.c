/**
 * @file geomean_sn.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief Runs geomean calculation from [Chum et al.] Sect 2
 * @date 2022-02-18
 * @copyright Copyright (c) 2022
 *
 * CITATIONS
 * -----------
 *  -  [Chum et al.] Chum, K., Guy, R. K., Jacobson, J. M. J., and Mosunov, A. S. (2018).
 *      Numerical and statistical analysis of aliquot sequences. Experimental Mathematics, 29(4):414–425.
 */

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../inc/moewsmoews_sieve.h"
#include "../inc/pomyang_kparent.h"
#include "../inc/sumdiv_util.h"
#include "../inc/math_macros.h"

void weighted_geomean(size_t bound, size_t seg_len);
void geomean(size_t bound, size_t seg_len);

/** @brief CLI to geometric mean calculation */
int main(int argc, char const *argv[]) {
    assert(argc == 3);
    size_t bound = strtol(argv[1], NULL, 10);
    size_t seg_len = strtol(argv[2], NULL, 10);
    weighted_geomean(bound, seg_len);
}

/** @brief Computes the geomean of s(n) for even n numbers. */
void geomean(size_t bound, size_t seg_len) {
    assert(EVEN(bound));
    assert(bound > 0);
    assert(EVEN(seg_len));
    assert(seg_len > 0);
    assert(seg_len <= bound);
    assert(DIVIDES(seg_len, bound));

    double accumulation = 1;
    size_t total = 0;
    for (size_t seg_start = 0; seg_start < bound; seg_start += seg_len) {
        for (uint64_t m = seg_start + 2; m <= seg_start + seg_len; m += 2) {
            uint64_t s_m = sumdiv_s(m);
            double abundance = (double)s_m / (double)m;

            accumulation += log(abundance);
            total++;
        }
    }
    double geomean = accumulation * (2.0 / bound);

    printf("geomean= %f\n", geomean);
}

/** @brief Computes the geomean of s(n) for even n numbers weighted by the number of aliquot parents.*/
void weighted_geomean(size_t bound, size_t seg_len) {
    assert(EVEN(bound));
    assert(bound > 0);
    assert(EVEN(seg_len));
    assert(seg_len > 0);
    assert(seg_len <= bound);
    assert(DIVIDES(seg_len, bound));

    pomyang_config cfg = {
        .preimage_count_bits = 8,
        .bound = bound,
        .seg_len = bound / 100,
        .num_locks = bound / 10,
        .num_threads = 12,
        .est_heap = 0,
        .quiet = 1,
    };

    PackedArray *f = pomyang_algorithm(&cfg);

    double accumulation = 1;
    size_t total = 0;
    size_t total_preimages = 0;

    for (size_t seg_start = 0; seg_start < bound; seg_start += seg_len) {
        for (uint64_t m = seg_start + 2; m <= seg_start + seg_len; m += 2) {
            uint64_t s_m = sumdiv_s(m);
            double abundance = (double)s_m / (double)m;
            uint8_t num_preimages = PackedArray_get(f, total);
            // num_preimages = (num_preimages == 0) ? 1 : num_preimages;
            total_preimages += num_preimages;
            accumulation += num_preimages * log(abundance);
            total++;
        }
    }

    double geomean = 1 + (accumulation * (1.0 / (double)total_preimages));
    printf("geomean= %f\n", geomean);
}
