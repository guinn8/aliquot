/**
 * @file profile_pom_yang.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief a scratch file to programatically run ranges of the pomerance-yang algorithm and dump to file
 * @date 2022-02-17
 * 
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 * 
 */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../inc/pom_yang.h"

#define OUTPUT_FILE "profile.csv"

int main(int argc, char const *argv[]) {
    (void)argc;
    (void)argv;
    const size_t min = 10000;
    const size_t max = 10000000000;
    size_t bound = min;
    while (bound <= max) {
        for (size_t i = 0; i < 10; i++) {
            size_t sub_bound = i * bound;
            for (size_t seg_len = 1; seg_len <= sub_bound; seg_len++) {
                size_t num_segs = sub_bound / seg_len;
                if (0 == sub_bound % seg_len &&
                    0 == seg_len % 2 &&
                    0 == (seg_len / 2) % 2 &&
                    num_segs >= 12 &&
                    num_segs <= (bound / 12)) {
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
                    print_to_file(&cfg, OUTPUT_FILE, count, (clock() - start));
                }
            }
        }

        bound = bound * 10;
    }
    return 0;
}
