/**
 * @file main.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief Implementation of algorithm to enumerate the number preimages under s(n) for all
 *        odd numbers under the user supplied bound
 * @date 2021-12-21
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "pom_yang.h"

#include <omp.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "../inc/properSumDiv.h"
#include "../inc/sieve_rewrite.h"

#define EVEN(x) (0 == (x) % 2)
#define SQUARE(x) ((x) * (x))
#define IS_M_PRIME(m, sigma_m) ((m) + 1 == (sigma_m))
#define SEG_OFFSET(seg_start, i) ((seg_start) + (2 * (i)) + 1)
#define TWO_THIRDS .6666666666666666666666
#define BYTES_TO_GB 0.000000001

#ifdef DEBUG_ASSERT_ON
#define DEBUG_ASSERT(x) x
#else
#define DEBUG_ASSERT(x)
#endif

static inline void set_sigma(uint64_t *m, uint64_t *sigma_m, const uint64_t set_m, const uint64_t set_sigma_m);
static inline void buffered_record_image(uint64_t x, PackedArray *f, uint64_t *writebuf, size_t writebuf_len, size_t *bufind);
static void flush_buf(PackedArray *f, uint64_t *writebuf, size_t *bufind);
static uint64_t *tabulate_aliquot_parents(PackedArray *f);
// static void print_config(const PomYang_config *cfg, double odd_comp_bound, size_t odd_comp_bound_seg);

PackedArray *Pomerance_Yang_aliquot(const PomYang_config *cfg) {
    assert(EVEN(cfg->bound));
    assert(EVEN(cfg->seg_len));
    assert(EVEN(cfg->seg_len / 2));
    assert(cfg->seg_len <= cfg->bound);
    assert(0 == cfg->bound % cfg->seg_len);
    assert(cfg->preimage_count_bits >= 1 && cfg->preimage_count_bits <= 8);

    const double odd_comp_bound = pow(cfg->bound, TWO_THIRDS);  // we need to compute values for all odd && composite numbers <= x^(2/3)
    const size_t odd_comp_bound_seg = ceil((float)odd_comp_bound / (float)cfg->seg_len) * cfg->seg_len;  // largest segment to sieve to reach max bound
    const size_t sigma_len = cfg->seg_len / 2;

    init_sigma_sieve(cfg->bound);
    PackedArray *f = PackedArray_create(cfg->preimage_count_bits, cfg->bound / 2);
#pragma omp parallel shared(f)
    {
        double thread_start = omp_get_wtime();
        uint64_t *sigma = (uint64_t *)calloc(sigma_len, sizeof(uint64_t));
        uint64_t *sigma_scratch = (uint64_t *)calloc(sigma_len, sizeof(uint64_t));
        uint64_t *writebuf = (uint64_t *)calloc(cfg->writebuf_len, sizeof(uint64_t));

        uint64_t m, sigma_m;
        size_t write_buf_ind = 0;

#pragma omp for schedule(guided)
        for (size_t seg_start = 0; seg_start < cfg->bound; seg_start += cfg->seg_len) {
            sigma_sieve_odd(cfg->seg_len, seg_start, sigma, sigma_scratch, 0);

            for (size_t i = 0; i < sigma_len; i++) {
                set_sigma(&m, &sigma_m, SEG_OFFSET(seg_start, i), sigma[i]);

                if (EVEN(sigma_m)) {
                    uint64_t t = (3 * sigma_m) - (2 * m);
                    while (t <= cfg->bound) {
                        buffered_record_image(t, f, writebuf, cfg->writebuf_len, &write_buf_ind);
                        t = (2 * t) + sigma_m;
                    }
                }

                if (IS_M_PRIME(m, sigma_m)) {
                    buffered_record_image(m + 1, f, writebuf, cfg->writebuf_len, &write_buf_ind);
                }
            }

            if (seg_start < odd_comp_bound_seg) {
                sigma_sieve_odd(cfg->seg_len, seg_start, sigma, sigma_scratch, 1);

                for (size_t i = 0; i < sigma_len; i++) {
                    set_sigma(&m, &sigma_m, SQUARE(SEG_OFFSET(seg_start, i)), sigma[i]);

                    if ((float)SEG_OFFSET(seg_start, i) > odd_comp_bound) {  // this test would be better pulled out of the loop (probably)
                        printf("\n %ld > %.2f (odd_comp_bound)reached by thread %d, breaking...\n\n",
                               SEG_OFFSET(seg_start, i), odd_comp_bound, omp_get_thread_num());
                        break;
                    }

                    uint64_t s_m = sigma_m - m;
                    if (!IS_M_PRIME(m, sigma_m) && (m > 1) && (s_m <= cfg->bound)) {
                        buffered_record_image(s_m, f, writebuf, cfg->writebuf_len, &write_buf_ind);
                    }
                }
            }
        }

        flush_buf(f, writebuf, &write_buf_ind);
        free(sigma);
        free(sigma_scratch);

        printf("thread %d completed standard in %.2fs\n", omp_get_thread_num(), omp_get_wtime() - thread_start);
    }

    destroy_sigma_sieve();
    return f;
}

uint64_t *count_kparent_aliquot(const PomYang_config *cfg) {
    PackedArray *f = Pomerance_Yang_aliquot(cfg);
    return tabulate_aliquot_parents(f);
}

inline void set_sigma(uint64_t *m, uint64_t *sigma_m, uint64_t set_m, uint64_t set_sigma_m) {
    *m = set_m;
    *sigma_m = set_sigma_m;
    DEBUG_ASSERT(assert(wheelDivSigma(*m) == *sigma_m));
}

void flush_buf(PackedArray *f, uint64_t *writebuf, size_t *bufind) {
    // double s = omp_get_wtime();
#pragma omp critical
    {
        for (size_t i = 0; i < *bufind; i++) {
            size_t offset = (writebuf[i] / 2) - 1;
            uint8_t count = PackedArray_get(f, offset);
            if (count + 1 < (1 << f->bitsPerItem)) {  // 2^f->bitsPerItem, prevents overflow
                PackedArray_set(f, offset, count + 1);
            }
        }
    }
    // printf("Thread %d drained buffer in %.2f, %.0f per second \n", omp_get_thread_num(), omp_get_wtime() - s, (double)*bufind / (omp_get_wtime() - s));

    *bufind = 0;
}

static inline void buffered_record_image(uint64_t x, PackedArray *f, uint64_t *writebuf, size_t writebuf_len, size_t *bufind) {
    DEBUG_ASSERT(assert(x > 0));
    DEBUG_ASSERT(assert(EVEN(x)));
    DEBUG_ASSERT(assert(*bufind + 1 < writebuf_len));

    writebuf[*bufind] = x;
    *bufind = *bufind + 1;

    if (*bufind + 1 == writebuf_len) {
        flush_buf(f, writebuf, bufind);
        *bufind = 0;
    }
}

uint64_t *tabulate_aliquot_parents(PackedArray *f) {
    double time_tabulate = omp_get_wtime();
    uint64_t *count = calloc(sizeof(uint64_t), UINT8_MAX + 1);  // we can have num_preimages == UINT8_MAX
    for (size_t i = 0; i < f->count; i++) {
        const uint8_t num_preimages = PackedArray_get(f, i);
        count[num_preimages]++;
    }

    printf("\nTabulation completed in %.2fs\n", omp_get_wtime() - time_tabulate);
    printf("\nCount of odd k-parent numbers under %d\n", f->count * 2);
    for (size_t i = 0; i < 8; i++) {
        printf("%ld: %ld\n", i, count[i]);
    }

    return count;
}

// void print_config(const PomYang_config *cfg, double odd_comp_bound, size_t odd_comp_bound_seg) {
//     printf("\nPomerance-Yang Algorithm Config Information\n");
//     printf("-> Using %ld bits per number, count upto 0-%d preimages\n", cfg->preimage_count_bits, (1 << cfg->preimage_count_bits) - 1);
//     printf("-> Bound = %ld\n", cfg->bound);
//     printf("-> Segment Length = %ld\n", cfg->seg_len);
//     printf("-> Number of segments = %ld\n", cfg->bound / cfg->seg_len);
//     printf("-> Max number of threads = %d\n", omp_get_max_threads());
//     printf("-> Bound^(2/3) = %.2f\n", odd_comp_bound);
//     printf("-> Bound^(2/3) Max Segment = %ld\n", odd_comp_bound_seg);

//     if (true == cfg->just_config) {
//         exit(EXIT_SUCCESS);
//     }
// }