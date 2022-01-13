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
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifndef NDEBUG
#include <flint/fmpz.h>
#endif  // NDEBUG

#include "../PackedArray/PackedArray.h"
#include "../inc/properSumDiv.h"
#include "../inc/sieve_rewrite.h"

#define EVEN(x) (0 == (x) % 2)
#define SQUARE(x) ((x) * (x))
#define IS_M_PRIME(m, sigma_m) ((m) + 1 == (sigma_m))
#define SEG_OFFSET(seg_start, i) ((seg_start) + (2 * (i)) + 1)
#define TWO_THIRDS .6666666666666666666666
#define BYTES_TO_GB 0.000000001
#define OUTPUT_FILE "counts.csv"

#ifndef NDEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif

uint64_t sigma(uint64_t n);
void set_sigma(uint64_t *m, uint64_t *sigma_m, const uint64_t set_m, const uint64_t set_sigma_m);
void record_image(uint64_t x, PackedArray *f);
void buffered_record_image(uint64_t x, PackedArray *f, uint64_t *buf, size_t buflen, size_t *bufind);
void flush_buf(PackedArray *f, uint64_t *buf, size_t *bufind);

int main(int argc, const char **argv) {
    if (argc != 3 && argc != 4 && argc != 5) {
        printf("USAGE: <bound><segment_len><(optional)num_threads><(optional)just_print_config>\n");
        exit(EXIT_FAILURE);
    }

    double start_time = omp_get_wtime();

    const size_t bound = strtol(argv[1], NULL, 10);
    const size_t seg_len = strtol(argv[2], NULL, 10);

    int just_config = 0;
    if (argc == 5) {
        just_config = strtol(argv[3], NULL, 10);
        omp_set_num_threads(strtol(argv[4], NULL, 10));
    }

    if (!EVEN(bound) ||
        !EVEN(seg_len) ||
        !EVEN(seg_len / 2) ||
        seg_len > bound ||
        0 != bound % seg_len) {
        printf("\nInvalid input!\n\n");
        exit(EXIT_FAILURE);
    }

    const double odd_comp_bound_float = pow(bound, TWO_THIRDS);
    const size_t odd_comp_bound = round(odd_comp_bound_float);                                 // ! unsure of best way to convert to int
    const size_t odd_comp_bound_seg = ceil((float)odd_comp_bound / (float)seg_len) * seg_len;  // largest segment to sieve to reach max bound

    const size_t f_len = bound / 2;
    const uint64_t f_bytes = f_len * sizeof(uint8_t);

    const size_t sigma_buf_len = seg_len / 2;
    const uint64_t sigma_buf_bytes = sigma_buf_len * sizeof(uint64_t);

    printf(
        "\nPomerance-Yang Algorithm\n-> Bound = %ld\n-> Segment Length = %ld\n-> Number of segments = %ld\n\n"
        "-> Bound^(2/3) = %.2f\n-> Integer(Bound^(2/3)) = %ld\n-> Bound^(2/3) Max Segment = %ld\n\n-> Max number of threads = %d\n",
        bound, seg_len, bound / seg_len, odd_comp_bound_float, odd_comp_bound, odd_comp_bound_seg, omp_get_max_threads());

    init_sigma_sieve(bound);

    uint64_t min_heap_alloc = f_bytes + 2 * sigma_buf_bytes * omp_get_max_threads();
    printf("\nHeap usage\n-> f: %.2fgb\n-> thread local buffers: %.2fgb\n-> Minimum heap allocation: %.2fgb\n",
           BYTES_TO_GB * f_bytes, BYTES_TO_GB * 2 * sigma_buf_bytes * omp_get_max_threads(), BYTES_TO_GB * min_heap_alloc);
    printf("Note: Flint allocates memory for itself which is excluded from this estimation.\n");

    if (just_config == 1) {
        exit(EXIT_SUCCESS);
    }

    // uint8_t *f = (uint8_t *)malloc(f_bytes);  // counts preimages for odd numbers
    PackedArray *f = PackedArray_create(1, f_len);
    // memset(f, 0, f_bytes);
    size_t write_buf_len = 5000000;


#pragma omp parallel shared(f)
    {
        double thread_start = omp_get_wtime();
        uint64_t *sigma_buf = (uint64_t *)malloc(sigma_buf_bytes);
        uint64_t *sieve_buf = (uint64_t *)malloc(sigma_buf_bytes);
        uint64_t m, sigma_m; 

        uint64_t *write_buf = (uint64_t *)calloc(write_buf_len, sizeof(uint64_t));
        size_t write_buf_ind = 0;

#pragma omp for schedule(dynamic)
        for (size_t seg_start = 0; seg_start < bound; seg_start += seg_len) {
            sigma_sieve_odd(seg_len, seg_start, sigma_buf, sieve_buf, 0);

            for (size_t i = 0; i < sigma_buf_len; i++) {
                set_sigma(&m, &sigma_m, SEG_OFFSET(seg_start, i), sigma_buf[i]);

                if (EVEN(sigma_m)) {
                    uint64_t t = (3 * sigma_m) - (2 * m);
                    while (t <= bound) {
                        buffered_record_image(t, f, write_buf, write_buf_len, &write_buf_ind);
                        t = (2 * t) + sigma_m;
                    }
                }

                if (IS_M_PRIME(m, sigma_m)) {
                    // record_image(m + 1, f);
                    buffered_record_image(m + 1, f, write_buf, write_buf_len, &write_buf_ind);

                }
            }
        }

        printf("thread %d completed standard in %.2f\n", omp_get_thread_num(), omp_get_wtime() - thread_start);
        thread_start = omp_get_wtime();

        /**
         * odd_comp_bound is not necessarily divisable by seg_len, to enumerate to the bound
         * we sieve to the next segment but halt reading the array once the bound has been reached 
         */
#pragma omp for schedule(dynamic)
        for (size_t seg_start = 0; seg_start < odd_comp_bound_seg; seg_start += seg_len) {
            sigma_sieve_odd(seg_len, seg_start, sigma_buf, sieve_buf, 1);

            for (size_t i = 0; i < sigma_buf_len; i++) {
                set_sigma(&m, &sigma_m, SQUARE(SEG_OFFSET(seg_start, i)), sigma_buf[i]);

                if (SEG_OFFSET(seg_start, i) > odd_comp_bound) {  // this test would be better pulled out of the loop (probably)
                    printf("\n %ld > %ld (odd_comp_bound)reached by thread %d, breaking...\n\n",
                           SEG_OFFSET(seg_start, i), odd_comp_bound, omp_get_thread_num());
                    break;
                }

                uint64_t s_m = sigma_m - m;
                if (!IS_M_PRIME(m, sigma_m) && (m > 1) && (s_m <= bound)) {
                    buffered_record_image(s_m, f, write_buf, write_buf_len, &write_buf_ind);

                    // record_image(s_m, f);
                }
            }
        }
        printf("thread %d completed squares in %.2f\n", omp_get_thread_num(), omp_get_wtime() - thread_start);

        flush_buf(f, write_buf, &write_buf_ind);
        free(sigma_buf);
        free(sieve_buf);
    }

    double time_tabulate = omp_get_wtime();
    size_t count[UINT8_MAX + 1] = {0};
    for (size_t i = 0; i < f_len; i++) {
        const uint8_t num_preimages = PackedArray_get(f, i);
        count[num_preimages]++;
    }
    printf("\nTabulation completed in %.2f\n", omp_get_wtime() - time_tabulate);

    printf("\nCount of odd k-parent numbers under %ld\n", bound);
    for (size_t i = 0; i < 8; i++) {
        printf("%ld: %ld\n", i, count[i]);
    }

    destroy_sigma_sieve();
    free(f);

    double end_time = omp_get_wtime();
    printf("\nCompleted in %.2f seconds\n\n", end_time - start_time);

    const size_t max_line_len = 4162;
    char *header_line = calloc(max_line_len, sizeof(char));
    if (0 != access(OUTPUT_FILE, F_OK)) {
        snprintf(header_line, max_line_len, "timestamp, bound, segment_length, timing, num_threads");
        for (size_t i = 0; i <= UINT8_MAX; i++) {
            snprintf(header_line + strlen(header_line), max_line_len - strlen(header_line), ", %ld", i);
        }
    }

    FILE *fp = fopen("counts.csv", "a");
    fprintf(fp, "%s\n", header_line);

    fprintf(fp, "%ld, ", time(NULL));              // timestamp
    fprintf(fp, "%ld, ", bound);                   // bound
    fprintf(fp, "%ld, ", seg_len);                 // segment length
    fprintf(fp, "%.2f, ", end_time - start_time);  // timing
    fprintf(fp, "%d", omp_get_max_threads());      // threads

    for (size_t i = 0; i <= UINT8_MAX; i++) {
        fprintf(fp, ", %ld", count[i]);
    }
    fclose(fp);

    exit(EXIT_SUCCESS);
}

#ifndef NDEBUG
uint64_t sigma(uint64_t n) {
    fmpz_t tmp = {0};
    fmpz_init_set_ui(tmp, n);
    fmpz_t sigma_tmp = {0};
    fmpz_divisor_sigma(sigma_tmp, tmp, 1);
    assert(fmpz_abs_fits_ui(sigma_tmp));
    return fmpz_get_ui(sigma_tmp);
}
#endif

inline void set_sigma(uint64_t *m, uint64_t *sigma_m, const uint64_t set_m, const uint64_t set_sigma_m) {
    *m = set_m;
    *sigma_m = set_sigma_m;
    // DEBUG(printf("sieved_sigma(%ld) = %ld, sigma(%ld) = %ld\n", *m, *sigma_m, *m, sigma(*m));)
    assert(sigma(*m) == *sigma_m);
}

inline void record_image(uint64_t x, PackedArray *f) {
    assert(x > 0);
    assert(EVEN(x));
    size_t offset = (x / 2) - 1;
#pragma omp critical
    PackedArray_set(f, offset, 1);
}

inline void flush_buf(PackedArray *f, uint64_t *buf, size_t *bufind) {
#pragma omp critical
    {
        for (size_t i = 0; i < *bufind; i++) {
            // printf("buf[%ld] = %ld\n", i, buf[i]);
            size_t offset = (buf[i] / 2) - 1;
            PackedArray_set(f, offset, 1);
        }
    }
    *bufind = 0;
}

inline void buffered_record_image(uint64_t x, PackedArray *f, uint64_t *buf, size_t buflen, size_t *bufind) {
    assert(x > 0);
    assert(EVEN(x));
    assert(*bufind + 1 < buflen);

    buf[*bufind] = x;
    *bufind = *bufind + 1;

    if (*bufind + 1 == buflen) {
        flush_buf(f, buf, bufind);
        *bufind = 0;
    }
}

