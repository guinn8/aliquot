/**
 * @file pom_yang_test.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief
 * @date 2022-02-06
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "../_pom_yang/pom_yang.h"

#include <assert.h>
// #include <omp.h>
#include <stdio.h>

#include "../inc/properSumDiv.h"

#define EVEN(x) (0 == (x) % 2)
#define F_OFFSET(x) ((x / 2) - 1)
#define F_DE_OFFSET(x) ((x + 1) * 2)

// I had to manually type these from Numerical and Statistical Analysis of Aliquot Sequences (Chum et al.)
#define _10_TO_THE_4     1212
#define _10_TO_THE_5     13863
#define _10_TO_THE_6     150232
#define _10_TO_THE_7     1574973
#define _10_TO_THE_8     16246940
#define _10_TO_THE_9     165826606
#define _10_TO_THE_10    1681871718
#define _10_TO_THE_11    16988116409
#define _10_TO_THE_12    171128671374
#define _2_TO_THE_40     188206399403

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

uint64_t *brute_force_preimage_counts(size_t bound);
uint8_t *brute_force_preimages(size_t bound);

int main(int argc, char **argv) {
    (void)argc;
    (void)argv;

    printf("Pomerance-Yang Test Suite, all tests have passed if the program does not assert out :)\n");
    // // bruteforce_kparent_aliquot_10_to_the_5();
    uint8_t *f_bruteforce = brute_force_preimages(100000);


    PomYang_config cfg = {
        .preimage_count_bits = 8,
        .bound = 100000,
        .seg_len = 100,
        .writebuf_len = 10,
        .num_threads = 1,
        .just_config = false,
        .quiet = true,
    };
    PackedArray *f = Pomerance_Yang_aliquot(&cfg);

    for (size_t i = 0; i < 100000 / 2; i++){
        if(f_bruteforce[i] < PackedArray_get(f, i)) {
            printf("f_bruteforce[%ld] = %d, PackedArray_get(f, %ld) = %d\n", F_DE_OFFSET(i), f_bruteforce[i], F_DE_OFFSET(i),  PackedArray_get(f, i));
        }
    }
    

    // numbits_1_10_to_the_5();
    // numbits_2_10_to_the_5();
    // numbits_3_10_to_the_5();
    // numbits_4_10_to_the_5();
    // numbits_5_10_to_the_5();
    // numbits_6_10_to_the_5();
    // numbits_7_10_to_the_5();
    // numbits_8_10_to_the_5();

    // writebuf_0_edge_10_to_the_6();
    // writebuf_1_edge_10_to_the_6();

    // seg_len_0_edge_10_to_the_6();
    // seg_len_1_edge_10_to_the_6();
    // seg_len_2_edge_10_to_the_6();

    simple_10_to_the_4();
    simple_10_to_the_5();
    simple_10_to_the_6();
    simple_10_to_the_7();
    simple_10_to_the_8();
    simple_10_to_the_9();
    simple_10_to_the_10();  // takes a long time but is needed
}

uint8_t *brute_force_preimages(size_t bound) {
    assert(EVEN(bound));

    const size_t f_len = bound / 2;
    uint8_t *f = calloc(f_len, sizeof(uint8_t));

#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 2; i <= 2 * bound; i += 2) {
            uint64_t s_n = wheelDivSum(i);
            if(s_n == 4) printf("s(%ld) = %ld\n", i, s_n);
            if (s_n <= bound && EVEN(s_n)) {
#pragma omp atomic
                f[F_OFFSET(s_n)]++;  // TODO(gavin): check overflow
            }
        }

#pragma omp for
        for (size_t i = 1; i <= bound; i += 2) {
            uint64_t s_n = wheelDivSum(i * i);
            if(s_n == 4) printf("2s(%ld) = %ld\n", i * i, s_n);
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

    for (size_t i = 0; i < 8; i++) {
        printf("%ld: %ld\n", i, count[i]);
    }

    return count;
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
    uint64_t *brute_force_count = brute_force_preimage_counts(100000);

    for (size_t i = 0; i < 8; i++) {
        printf("count[%ld] = %ld, brute_force_count[%ld] = %ld\n", i, count[i], i, brute_force_count[i]);
        // assert(count[i] == brute_force_count[i]);
    }

    assert(count[0] + 1 == _10_TO_THE_5);
    free(count);
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
        .bound = 1000000000,
        .seg_len = 50000000,
        .writebuf_len = 1000000,
        .num_threads = 8,
        .just_config = false,
        .quiet = true,
    };

    uint64_t *count = count_kparent_aliquot(&cfg);
    assert(count[0] + 1 == _10_TO_THE_9);
    free(count);
}
