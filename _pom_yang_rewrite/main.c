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

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

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

#ifdef DEBUG_ASSERT_ON
#define DEBUG_ASSERT(x) x
#else
#define DEBUG_ASSERT(x)
#endif

typedef struct {
    char name[32];
    size_t len;
    size_t item_size_bits;
    size_t num_instances;
} buf_info_t;

typedef struct {
    const buf_info_t *info;
    uint64_t *buf;
} u64_buf_t;

static struct option long_options[] = {
    {"bound", required_argument, NULL, 'b'},
    {"seg_len", required_argument, NULL, 's'},
    {"writebuf_len", required_argument, NULL, 'w'},
    {"just_config", no_argument, NULL, 'j'},
    {"num_threads", required_argument, NULL, 't'},
    {"preimage_count_bits", required_argument, NULL, 'p'},
    {0, 0, 0, 0}};

void set_sigma(uint64_t *m, uint64_t *sigma_m, const uint64_t set_m, const uint64_t set_sigma_m);
void record_image(uint64_t x, PackedArray *f);
void buffered_record_image(uint64_t x, PackedArray *f, u64_buf_t *writebuf, size_t *bufind);
void flush_buf(PackedArray *f, u64_buf_t *writebuf, size_t *bufind);
size_t print_array_info(buf_info_t x);

void usage(void) {
    printf("\nFast and Memory effiecent implementation of the Pomerance-Yang algorithm enumerating preimages under s()\n");
    printf("Usage: [--bound=][--seg_len=][--writebuf_len=][(OPTIONAL)--just_config][(OPTIONAL)--num_threads=][(OPTIONAL)--preimage_count_bits]\n\n");
}

int main(int argc, char **argv) {
    size_t preimage_count_bits = 8;
    size_t bound = 0;
    size_t seg_len = 0;
    size_t writebuf_len = 0;
    bool just_config = 0;

    int opt_code, option_index;
    while (-1 != (opt_code = getopt_long(argc, argv, "", long_options, &option_index))) {
        switch (opt_code) {
            case 'b':
                bound = strtol(optarg, NULL, 10);
                break;
            case 's':
                seg_len = strtol(optarg, NULL, 10);
                break;
            case 'w':
                writebuf_len = strtol(optarg, NULL, 10);
                break;
            case 'j':
                just_config = true;
                break;
            case 't':
                omp_set_num_threads(strtol(optarg, NULL, 10));
                break;
            case 'p':
                preimage_count_bits = strtol(optarg, NULL, 10);
                break;
            default:
                usage();
                return EXIT_FAILURE;
        }
    }

    if (0 == bound ||
        0 == seg_len ||
        0 == writebuf_len) {  // required args
        usage();
        return EXIT_FAILURE;
    }

    assert(EVEN(bound));
    assert(EVEN(seg_len));
    assert(EVEN(seg_len / 2));
    assert(seg_len <= bound);
    assert(0 == bound % seg_len);
    assert(preimage_count_bits >= 1 && preimage_count_bits <= 8);

    // this is only for the heap estimation, use the PackedArray helpers in practice
    const buf_info_t f_info = {
        .name = "f",
        .num_instances = 1,
        .item_size_bits = preimage_count_bits,
        .len = bound / 2,
    };
    const buf_info_t sigma_info = {
        .name = "sigma",
        .num_instances = omp_get_max_threads(),
        .item_size_bits = sizeof(uint64_t) * 8,  // ! measuring by bits
        .len = seg_len / 2,
    };
    const buf_info_t sieve_info = {
        .name = "sieve",
        .num_instances = omp_get_max_threads(),
        .item_size_bits = sizeof(uint64_t) * 8,  // ! measuring by bits
        .len = seg_len / 2,
    };
    const buf_info_t writebuf_info = {
        .name = "writebuf",
        .num_instances = omp_get_max_threads(),
        .item_size_bits = sizeof(uint64_t) * 8,  // ! measuring by bits
        .len = writebuf_len,
    };

    const double odd_comp_bound_float = pow(bound, TWO_THIRDS);  // we need to compute values for all odd && composite numbers <= x^(2/3)
    const size_t odd_comp_bound_seg = ceil((float)odd_comp_bound_float / (float)seg_len) * seg_len;  // largest segment to sieve to reach max bound

    printf("\nPomerance-Yang Algorithm Config Information\n");
    printf("-> Using %ld bits per number, count upto 0-%d preimages\n", preimage_count_bits, (1 << preimage_count_bits) - 1);
    printf("-> Bound = %ld\n", bound);
    printf("-> Segment Length = %ld\n", seg_len);
    printf("-> Number of segments = %ld\n", bound / seg_len);
    printf("-> Bound^(2/3) = %.2f\n", odd_comp_bound_float);
    printf("-> Bound^(2/3) Max Segment = %ld\n", odd_comp_bound_seg);
    printf("-> Max number of threads = %d\n", omp_get_max_threads());

    size_t min_heap_alloc = 0;
    printf("\nEstimated Minimum Heap Allocation\n");
    min_heap_alloc += print_array_info(f_info);
    min_heap_alloc += print_array_info(sigma_info);
    min_heap_alloc += print_array_info(sieve_info);
    min_heap_alloc += print_array_info(writebuf_info);
    printf("This configuration requires a minimum of %0.2fgb\n\n", BYTES_TO_GB * min_heap_alloc);

    if (true == just_config) {
        exit(EXIT_SUCCESS);
    }

    double start_time = omp_get_wtime();
    init_sigma_sieve(bound);
    PackedArray *f = PackedArray_create(f_info.item_size_bits, f_info.len);

#pragma omp parallel shared(f)
    {
        double thread_start = omp_get_wtime();
        u64_buf_t sigma = {
            .info = &sigma_info,
            .buf = (uint64_t *)calloc(sigma_info.len, sigma_info.item_size_bits / 8), 
        };
        u64_buf_t sieve = {
            .info = &sieve_info,
            .buf = (uint64_t *)calloc(sieve_info.len, sieve_info.item_size_bits / 8),
        };
        u64_buf_t writebuf = {
            .info = &writebuf_info,
            .buf = (uint64_t *)calloc(writebuf_info.len, writebuf_info.item_size_bits / 8),
        };

        uint64_t m, sigma_m;
        size_t write_buf_ind = 0;

#pragma omp for schedule(dynamic)
        for (size_t seg_start = 0; seg_start < bound; seg_start += seg_len) {
            sigma_sieve_odd(seg_len, seg_start, sigma.buf, sieve.buf, 0);

            for (size_t i = 0; i < sigma_info.len; i++) {
                set_sigma(&m, &sigma_m, SEG_OFFSET(seg_start, i), sigma.buf[i]);

                if (EVEN(sigma_m)) {
                    uint64_t t = (3 * sigma_m) - (2 * m);
                    while (t <= bound) {
                        buffered_record_image(t, f, &writebuf, &write_buf_ind);
                        t = (2 * t) + sigma_m;
                    }
                }

                if (IS_M_PRIME(m, sigma_m)) {
                    buffered_record_image(m + 1, f, &writebuf, &write_buf_ind);
                }
            }
        }

        printf("thread %d completed standard in %.2fs\n", omp_get_thread_num(), omp_get_wtime() - thread_start);
        thread_start = omp_get_wtime();

        // odd_comp_bound is not necessarily divisable by seg_len, to enumerate to the bound
        // we sieve to the next segment but halt reading the array once the bound has been reached
#pragma omp for schedule(dynamic)
        for (size_t seg_start = 0; seg_start < odd_comp_bound_seg; seg_start += seg_len) {
            sigma_sieve_odd(seg_len, seg_start, sigma.buf, sieve.buf, 1);

            for (size_t i = 0; i < sigma_info.len; i++) {
                set_sigma(&m, &sigma_m, SQUARE(SEG_OFFSET(seg_start, i)), sigma.buf[i]);

                if ((float)SEG_OFFSET(seg_start, i) > odd_comp_bound_float) {  // this test would be better pulled out of the loop (probably)
                    printf("\n %ld > %.2f (odd_comp_bound_float)reached by thread %d, breaking...\n\n",
                           SEG_OFFSET(seg_start, i), odd_comp_bound_float, omp_get_thread_num());
                    break;
                }

                uint64_t s_m = sigma_m - m;
                if (!IS_M_PRIME(m, sigma_m) && (m > 1) && (s_m <= bound)) {
                    buffered_record_image(s_m, f, &writebuf, &write_buf_ind);
                }
            }
        }
        printf("thread %d completed squares in %.2f\n", omp_get_thread_num(), omp_get_wtime() - thread_start);

        flush_buf(f, &writebuf, &write_buf_ind);
        free(sigma.buf);
        free(sieve.buf);
    }

    double time_tabulate = omp_get_wtime();
    size_t count[UINT8_MAX + 1] = {0};
    for (size_t i = 0; i < f_info.len; i++) {
        const uint8_t num_preimages = PackedArray_get(f, i);
        count[num_preimages]++;
    }
    printf("\nTabulation completed in %.2fs\n", omp_get_wtime() - time_tabulate);

    printf("\nCount of odd k-parent numbers under %ld\n", bound);
    for (size_t i = 0; i < 8; i++) {
        printf("%ld: %ld\n", i, count[i]);
    }

    destroy_sigma_sieve();
    free(f);

    double end_time = omp_get_wtime();
    printf("\nCompleted in %.2fs\n\n", end_time - start_time);

    const size_t max_line_len = 65536;
    char *header_line = calloc(max_line_len, sizeof(char));
    if (0 != access(OUTPUT_FILE, F_OK)) {
        snprintf(header_line, max_line_len, "timestamp, bound, segment_length, timing, num_threads");
        for (size_t i = 0; i <= UINT8_MAX; i++) {
            snprintf(header_line + strlen(header_line), max_line_len - strlen(header_line), ", %ld", i);
        }
    }

    FILE *fp = fopen("counts.csv", "a");
    fprintf(fp, "%s\n", header_line);

    fprintf(fp, "%ld, ", time(NULL));
    fprintf(fp, "%ld, ", bound);
    fprintf(fp, "%ld, ", seg_len);
    fprintf(fp, "%.2f, ", end_time - start_time);
    fprintf(fp, "%d", omp_get_max_threads());

    for (size_t i = 0; i <= UINT8_MAX; i++) {
        fprintf(fp, ", %ld", count[i]);
    }
    fclose(fp);

    exit(EXIT_SUCCESS);
}

inline void set_sigma(uint64_t *m, uint64_t *sigma_m, uint64_t set_m, uint64_t set_sigma_m) {
    *m = set_m;
    *sigma_m = set_sigma_m;
    DEBUG_ASSERT(assert(wheelDivSigma(*m) == *sigma_m));
}

inline void flush_buf(PackedArray *f, u64_buf_t *writebuf, size_t *bufind) {
#pragma omp critical
    {
        for (size_t i = 0; i < *bufind; i++) {
            size_t offset = (writebuf->buf[i] / 2) - 1;
            uint8_t count = PackedArray_get(f, offset);
            if (count + 1 < (1 << f->bitsPerItem)) {  // 2^f->bitsPerItem, prevents overflow
                PackedArray_set(f, offset, count + 1);
            }
        }
    }

    *bufind = 0;
}

inline void buffered_record_image(uint64_t x, PackedArray *f, u64_buf_t *writebuf, size_t *bufind) {
    DEBUG_ASSERT(assert(x > 0));
    DEBUG_ASSERT(assert(EVEN(x)));
    DEBUG_ASSERT(assert(*bufind + 1 < writebuf->info->len));

    writebuf->buf[*bufind] = x;
    *bufind = *bufind + 1;

    if (*bufind + 1 == writebuf->info->len) {
        flush_buf(f, writebuf, bufind);
        *bufind = 0;
    }
}

size_t print_array_info(const buf_info_t x) {
    size_t size_bytes = (x.item_size_bits * x.len) / 8;
    size_t total_bytes = size_bytes * x.num_instances;
    printf("%s: num_instances=%ld len=%ld size_bytes=%ld total_gb=%0.2f\n",
           x.name, x.num_instances, x.len, size_bytes, BYTES_TO_GB * total_bytes);
    return total_bytes;
}
