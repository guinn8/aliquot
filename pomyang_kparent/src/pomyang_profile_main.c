/**
 * @cond Doxygen_Suppress
 * @file pomyang_profile_main.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief a scratch file to programatically run ranges of the pomerance-yang algorithm and dump to file
 * @date 2022-02-17
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../inc/pomyang_kparent.h"

void run_seglen_perms(size_t sub_bound, const char *filename);
void run_ordermag_intervals(size_t bound, const char *filename);

int main(int argc, char const *argv[]) {
    assert(argc == 4);
    size_t minbound = strtol(argv[1], NULL, 10);
    printf("minbound: %ld\n", minbound);
    size_t maxbound = strtol(argv[2], NULL, 10);
    const char *filename = argv[3];
    size_t bound = minbound;
    while (bound <= maxbound) {
        // pomyang_config cfg = {
        //     .preimage_count_bits = 8,
        //     .bound = bound,
        //     .seg_len = bound / 2,
        //     .num_locks = bound / 2,
        //     .num_threads = 12,
        //     .est_heap = 0,
        //     .quiet = 1,
        // };

        // clock_t start = clock();
        // uint64_t *count = pomyang_count_kparent(&cfg);
        // print_to_file(&cfg, filename, count, (clock() - start));
        // bound *= 2;
        run_ordermag_intervals(bound, filename);

        bound *= 10;
    }
    return 0;
}

void run_ordermag_intervals(size_t bound, const char *filename) {
    for (size_t i = 1; i < 10; i++) {
        size_t sub_bound = i * bound;
        pomyang_config cfg = {
            .preimage_count_bits = 8,
            .bound = sub_bound,
            .seg_len = sub_bound,
            .num_locks = sub_bound / 10,
            .num_threads = 12,
            .est_heap = 0,
            .quiet = 1,
        };

        clock_t start = clock();
        uint64_t *count = pomyang_count_kparent(&cfg);
        print_to_file(&cfg, filename, count, (clock() - start));
    }
    bound = bound * 10;
}

void run_seglen_perms(size_t sub_bound, const char *filename) {
    for (size_t seg_len = 1; seg_len <= sub_bound; seg_len++) {
        size_t num_segs = sub_bound / seg_len;
        if (0 == sub_bound % seg_len &&
            0 == seg_len % 2 &&
            0 == (seg_len / 2) % 2 &&
            num_segs >= 12 &&
            num_segs <= (sub_bound / 12)) {
            pomyang_config cfg = {
                .preimage_count_bits = 1,
                .bound = sub_bound,
                .seg_len = seg_len,
                .num_locks = sub_bound / 100,
                .num_threads = 12,
                .est_heap = 0,
                .quiet = 1,
            };

            clock_t start = clock();
            uint64_t *count = pomyang_count_kparent(&cfg);
            print_to_file(&cfg, filename, count, (clock() - start));
        }
    }
}
/** @endcond */
