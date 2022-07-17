/**
 * @file pomyang_kparent.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief api for pom_yang algorithm operations
 * @date 2021-2-17
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#ifndef POMYANG_KPARENT_INC_POMYANG_KPARENT_H_
#define POMYANG_KPARENT_INC_POMYANG_KPARENT_H_

#include <stdbool.h>
#include <stdlib.h>

#include "../inc/PackedArray.h"
#include "../inc/math_macros.h"

/**
 * @brief Configure the Pomerance-Yang algorithm.
 * @struct pomyang_config_t
 * 
 * @var pomyang_config_t::preimage_count_bits
 * Input to the PackedArray data struct, determines the actual number of bits used 
 * to count the number of preimages some number has. Can be 1, 2, 4, or 8.
 * 
 * @var pomyang_config_t::bound
 * Determine counts of all even k-parent numbers under this bound. Must be even.
 * 
 * @var pomyang_config_t::seg_len
 * Segment length for sieving blocks of sigma(n), the algorithms input.
 * 
 * @var pomyang_config_t::num_locks
 * The size of the mutex buffer used to protect the shared counter buffer during multi-threaded access.
 * 
 * @var pomyang_config_t::num_threads
 * How many threads to use.
 * 
 * @var pomyang_config_t::est_heap
 * Estimates heap usage, useful for large bound runs.
 * 
 * @var pomyang_config_t::quiet
 * Quiets some logging.
 */
typedef struct {
    size_t preimage_count_bits;
    size_t bound;
    size_t seg_len;
    size_t num_locks;
    size_t num_threads;
    bool est_heap;
    bool quiet;
    char filename[32];
} pomyang_config_t;

/**
 * @brief Runs the Pomerance-Yang algorithm.
 *
 * @param cfg See struct defintion
 * @return PackedArray* containing the number of preimages for even number upto bound. Free'd by caller.
 */
PackedArray *pomyang_algorithm(const pomyang_config_t *cfg);

/**
 * @brief runs the Pomerance-Yang algorithm and counts occurrence of kparent numbers.
 *
 * @param cfg See struct definition
 * @return uint64_t* buffer of length (UINT8_MAX + 1) with occurrence counts, must be free'd by caller
 */
uint64_t *pomyang_count_kparent(const pomyang_config_t *cfg);

/**
 * @brief Prints Pomerance-Yang algorithm configuration.
 * 
 * @param cfg Pointer to config struct.
 * @param filename Name of file to write.
 * @param count Array of k-parent counts.
 * @param runtime CPU seconds used.
 */
void print_to_file(pomyang_config_t *cfg, const char *filename, uint64_t *count, float runtime);
#endif  // POMYANG_KPARENT_INC_POMYANG_KPARENT_H_
