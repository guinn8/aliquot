/**
 * @file pom_yang.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief api for pom_yang algorithm operations
 * @date 2021-2-17
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#ifndef _POM_YANG_INC_POM_YANG_H_
#define _POM_YANG_INC_POM_YANG_H_

#include <stdbool.h>
#include <stdlib.h>

#include "../../PackedArray/PackedArray.h"

// I had to manually type these non-aliquot counts from Numerical and Statistical Analysis of Aliquot Sequences (Chum et al.)
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"  // falsely reports unused
static char *preimage_count_bits_desc =
    "preimage_count_bits: input to the PackedArray data struct, determines the actual number of bits used...\n"
    "to count the number of preimages some number has. Can be 1, 2, 4, or 8.";
static char *bound_desc = "bound: determine counts of all even k-parent numbers under this bound. Must be even.";
static char *seg_len_desc = "seg_len: segment length for sieving blocks of sigma(n), the algorithms input.";
static char *num_locks_desc = "num_locks: the size of the mutex buffer used to protect the shared counter buffer during multi-threaded access.";
static char *num_threads_desc = "num_threads: how many threads to use.";
static char *est_heap_desc = "est_heap: estimates heap usage, useful for large bound runs.";
static char *quiet_desc = "quiet: quiets some logging.";
#pragma GCC diagnostic pop

typedef struct {
    size_t preimage_count_bits;
    size_t bound;
    size_t seg_len;
    size_t num_locks;
    size_t num_threads;
    bool est_heap;
    bool quiet;
} pomyang_config;

/**
 * @brief runs the Pomerance-Yang algorithm.
 *
 * @param cfg See struct defintion
 * @return PackedArray* containing the number of preimages for even number upto bound. Free'd by caller.
 */
PackedArray *pomyang_algorithm(const pomyang_config *cfg);

/**
 * @brief runs the Pomerance-Yang algorithm and counts occurrence of kparent numbers.
 *
 * @param cfg See struct definition
 * @return uint64_t* buffer of length (UINT8_MAX + 1) with occurrence counts, must be free'd by caller
 */
uint64_t *pomyang_count_kparent(const pomyang_config *cfg);


void print_to_file(pomyang_config *cfg, char *filename, uint64_t *count, float runtime);


#endif  // _POM_YANG_INC_POM_YANG_H_
