/**
 * @file pom_yang_test.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief
 * @date 2022-02-06
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <assert.h>
#include <stdio.h>

#include "../inc/brute_force_preimages.h"
#include "../inc/pom_yang.h"

// I had to manually type these from Numerical and Statistical Analysis of Aliquot Sequences (Chum et al.)
// these are non aliquot counts
#define _10_TO_THE_4 1212
#define _10_TO_THE_5 13863
#define _10_TO_THE_6 150232
#define _10_TO_THE_7 1574973
#define _10_TO_THE_8 16246940
#define _10_TO_THE_9 165826606
#define _10_TO_THE_10 1681871718
#define _10_TO_THE_11 16988116409
#define _10_TO_THE_12 171128671374
#define _2_TO_THE_40 188206399403

void simple_10_to_the_4(void);
void simple_10_to_the_5(void);
void simple_10_to_the_6(void);
void simple_10_to_the_7(void);
void simple_10_to_the_8(void);
void simple_10_to_the_9(void);
void simple_10_to_the_10(void);

void seg_len_0_edge_10_to_the_6(void);
void seg_len_1_edge_10_to_the_6(void);
void seg_len_2_edge_10_to_the_6(void);

void writebuf_0_edge_10_to_the_6(void);
void writebuf_1_edge_10_to_the_6(void);

void numbits_1_10_to_the_5(void);
void numbits_2_10_to_the_5(void);
void numbits_3_10_to_the_5(void);
void numbits_4_10_to_the_5(void);
void numbits_5_10_to_the_5(void);
void numbits_6_10_to_the_5(void);
void numbits_7_10_to_the_5(void);
void numbits_8_10_to_the_5(void);

void bruteforce_kparent_aliquot_10_to_the_5(void);
void bruteforce_kparent_aliquot_2_to_the_20(void);
void bruteforce_kparent_aliquot_2_to_the_3(void);
void bruteforce_kparent_aliquot_24(void);
void bruteforce_kparent_aliquot_10024(void);

int main(int argc, char **argv) {
    (void)argc, (void)argv;

    printf("Pomerance-Yang Test Suite, all tests have passed if the program does not assert out :)\n");
    bruteforce_kparent_aliquot_10_to_the_5();
    bruteforce_kparent_aliquot_2_to_the_20();
    bruteforce_kparent_aliquot_2_to_the_3();
    bruteforce_kparent_aliquot_24();
    bruteforce_kparent_aliquot_10024();

    numbits_1_10_to_the_5();
    numbits_2_10_to_the_5();
    numbits_3_10_to_the_5();
    numbits_4_10_to_the_5();
    numbits_5_10_to_the_5();
    numbits_6_10_to_the_5();
    numbits_7_10_to_the_5();
    numbits_8_10_to_the_5();

    writebuf_0_edge_10_to_the_6();
    writebuf_1_edge_10_to_the_6();

    seg_len_0_edge_10_to_the_6();
    seg_len_1_edge_10_to_the_6();
    seg_len_2_edge_10_to_the_6();

    simple_10_to_the_4();
    simple_10_to_the_5();
    simple_10_to_the_6();
    simple_10_to_the_7();
    simple_10_to_the_8();
    simple_10_to_the_9();
    simple_10_to_the_10();  // takes a long time but is needed
}

void bruteforce_kparent_aliquot_10024(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 10024,
        .seg_len = 56,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    uint64_t *brute_force_count = brute_force_preimage_counts(cfg.bound);

    for (size_t i = 0; i < UINT8_MAX; i++) {
        assert(count[i] == brute_force_count[i]);
    }

}

void bruteforce_kparent_aliquot_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    uint64_t *brute_force_count = brute_force_preimage_counts(cfg.bound);

    for (size_t i = 0; i < UINT8_MAX; i++) {
        assert(count[i] == brute_force_count[i]);
    }

    assert(count[0] + 1 == _10_TO_THE_5);
}

void bruteforce_kparent_aliquot_2_to_the_20(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 262144,
        .seg_len = 65536,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    uint64_t *brute_force_count = brute_force_preimage_counts(cfg.bound);

    for (size_t i = 0; i < UINT8_MAX; i++) {
        assert(count[i] == brute_force_count[i]);
    }
}

void bruteforce_kparent_aliquot_2_to_the_3(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 8,
        .seg_len = 4,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    uint64_t *brute_force_count = brute_force_preimage_counts(cfg.bound);

    for (size_t i = 0; i < UINT8_MAX; i++) {
        assert(count[i] == brute_force_count[i]);
    }
}

void bruteforce_kparent_aliquot_24(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 24,
        .seg_len = 8,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    uint64_t *brute_force_count = brute_force_preimage_counts(cfg.bound);

    for (size_t i = 0; i < UINT8_MAX; i++) {
        assert(count[i] == brute_force_count[i]);
    }
}

void numbits_1_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_2_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 2,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_3_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 3,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_4_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 4,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_5_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 5,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_6_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 6,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_7_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 7,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void numbits_8_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 100000,
        .seg_len = 10000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void writebuf_0_edge_10_to_the_6(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000,
        .seg_len = 10000,
        .writebuf_len = 1,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_6);
    free(count);
}

void writebuf_1_edge_10_to_the_6(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000,
        .seg_len = 10000,
        .writebuf_len = 9999999,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_6);
    free(count);
}

void seg_len_0_edge_10_to_the_6(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000,
        .seg_len = 4,
        .writebuf_len = 100,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_6);
    free(count);
}

void seg_len_1_edge_10_to_the_6(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000,
        .seg_len = 1000000,
        .writebuf_len = 100,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_6);
    free(count);
}

void seg_len_2_edge_10_to_the_6(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000,
        .seg_len = 62500,
        .writebuf_len = 100,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_6);
    free(count);
}

void simple_10_to_the_4(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 10000,
        .seg_len = 100,
        .writebuf_len = 100,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_4);
    free(count);
}

void simple_10_to_the_5(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 100000,
        .seg_len = 1000,
        .writebuf_len = 1000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
}

void simple_10_to_the_6(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000,
        .seg_len = 10000,
        .writebuf_len = 10000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_6);
    free(count);
}

void simple_10_to_the_7(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 10000000,
        .seg_len = 100000,
        .writebuf_len = 100000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_7);
    free(count);
}

void simple_10_to_the_8(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 100000000,
        .seg_len = 1000000,
        .writebuf_len = 100000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_8);
    free(count);
}

void simple_10_to_the_9(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 1000000000,
        .seg_len = 5000000,
        .writebuf_len = 1000000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_9);
    free(count);
}

// this one will take a long time, but is required.
// it is the first simple test greater than 2^32
void simple_10_to_the_10(void) {
    PomYang_config cfg = {
        .preimage_count_bits = 1,
        .bound = 10000000000,
        .seg_len = 50000000,
        .writebuf_len = 1000000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_10);
    free(count);
}
