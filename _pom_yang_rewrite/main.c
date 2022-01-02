/**
 * @file main.cpp
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2021-12-21
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <assert.h>
#include <flint/fmpz.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>

#include "../inc/sieve_rewrite.h"

#define EVEN(x) (0 == (x) % 2)
#define SQUARE(x) ((x) * (x))
#define IS_M_PRIME(m, sigma_m) ((m) + 1 == (sigma_m))
#define SEG_OFFSET(seg_start, i) ((seg_start) + (2 * (i)) + 1)
#define TWO_THIRDS .6666666666666666666666

#ifndef NDEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif

uint64_t sigma(uint64_t n);
void set_sigma(uint64_t *m, uint64_t *sigma_m, const uint64_t set_m, const uint64_t set_sigma_m);
void record_image(uint64_t x, uint8_t *f);

int main(int argc, const char **argv) {
    if (argc != 3) {
        printf("USAGE: <bound><segment_len>\n");
        exit(EXIT_FAILURE);
    }

    const size_t bound = strtol(argv[1], NULL, 10);
    const size_t seg_len = strtol(argv[2], NULL, 10);
    assert(EVEN(bound));
    assert(EVEN(seg_len));
    assert(EVEN(seg_len / 2));
    assert(seg_len <= bound);
    assert(0 == bound % seg_len);

    const double odd_comp_bound_float = pow(bound, TWO_THIRDS);
    const size_t odd_comp_bound = round(odd_comp_bound_float);  // ! unsure of best way to convert to int

    uint8_t *f = (uint8_t *)calloc(bound / 2, sizeof(uint8_t));  // counts preimages for odd numbers
    init_sigma_sieve(bound);

    printf("Pomerance-Yang Algorithm\n-> Bound = %ld\n-> Segment Length = %ld\n-> Number of segments = %ld\n"
           "-> Bound^(2/3) = %f\n-> Integer(Bound^(2/3)) = %ld\n",
            bound, seg_len, bound/seg_len, odd_comp_bound_float, odd_comp_bound);

    #pragma omp parallel shared(f)
    {
        uint64_t *sigma_buf = (uint64_t *)calloc(seg_len / 2, sizeof(uint64_t));
        uint64_t m, sigma_m;

        #pragma omp for schedule(dynamic)
        for (size_t seg_start = 0; seg_start < bound; seg_start += seg_len) {
            sigma_sieve_odd(seg_len, seg_start, sigma_buf, 0);

            for (size_t i = 0; i < seg_len / 2; i++) {
                set_sigma(&m, &sigma_m, SEG_OFFSET(seg_start, i), sigma_buf[i]);
                if (EVEN(sigma_m)) {
                    uint64_t t = (3 * sigma_m) - (2 * m);
                    while (t <= bound) {
                        record_image(t, f);
                        t = (2 * t) + sigma_m;
                    }
                }

                if (IS_M_PRIME(m, sigma_m)) {
                    record_image(m + 1, f);
                }
            }
        }

        // #pragma omp for nowait schedule(dynamic)
        // for (size_t seg_start = 0; seg_start < odd_comp_bound; seg_start += seg_len) {
        //     sigma_sieve_odd(seg_len, seg_start, sigma_buf, 1);

        //     for (size_t i = 0; i < seg_len / 2; i++) {
        //         set_sigma(&m, &sigma_m, SQUARE(SEG_OFFSET(seg_start, i)), sigma_buf[i]);
        //         if (!IS_M_PRIME(m, sigma_m) && (m > 1) && (sigma_m - m <= bound)) {
        //             record_image(sigma_m - m, f);
        //         }
        //     }
        // }
        #pragma omp for schedule(dynamic)
        for (size_t m = 1; m < odd_comp_bound; m += 2) {
            set_sigma(&m, &sigma_m, SQUARE(m), sigma(SQUARE(m)));
            uint64_t s_m = sigma_m - m;
            if (!IS_M_PRIME(m, sigma_m) && (m > 1) && (s_m <= bound)) {
                record_image(s_m, f);
            }
        }
        free(sigma_buf);
    }

    size_t count[UINT8_MAX] = {0};
    for (size_t i = 0; i < bound / 2; i++) {
        count[f[i]]++;
    }

    printf("\nCount of odd k-parent numbers under %ld\n", bound);
    for (size_t i = 0; i < 8; i++) {
        printf("%ld: %ld\n", i, count[i]);
    }

    free(f);
    exit(EXIT_SUCCESS);
}

uint64_t sigma(uint64_t n) {
    fmpz_t tmp = {0};
    fmpz_init_set_ui(tmp, n);
    fmpz_t sigma_tmp = {0};
    fmpz_divisor_sigma(sigma_tmp, tmp, 1);
    assert(fmpz_abs_fits_ui(sigma_tmp));
    return fmpz_get_ui(sigma_tmp);
}

inline void set_sigma(uint64_t *m, uint64_t *sigma_m, const uint64_t set_m, const uint64_t set_sigma_m) {
    *m = set_m;
    *sigma_m = set_sigma_m;
    DEBUG(printf("sieved_sigma(%ld) = %ld, sigma(%ld) = %ld\n", *m, *sigma_m, *m, sigma(*m));)
    assert(sigma(*m) == *sigma_m);
}

inline void record_image(uint64_t x, uint8_t *f) {
    assert(x > 0);
    assert(EVEN(x));
    size_t offset = (x / 2) - 1;
    #pragma omp atomic
        f[offset]++;
}
