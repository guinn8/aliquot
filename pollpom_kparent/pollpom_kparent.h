/**
 * @file pollpom_kparent.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-03-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdint.h>
#include <stdlib.h>

/**
 * @brief Config struct for the kparent generalization of the Pollack-Pomerance heuristic model for non-aliquots
 * @struct pollpom_config_t
 * 
 * @var pollpom_config_t::accum
 * Array where index k holds the summation for kparent numbers.
 * 
 * @var pollpom_config_t::bound
 * Compute the limit upto this bound.
 * 
 * @var pollpom_config_t::num_chunks
 * Computed value, the number of chunks to sieve. No sieve currently implemented.
 * 
 * @var pollpom_config_t::chunk_len
 * Size of chuck to sieve at a time. No sieve currently implemented. 
 * 
 * @var pollpom_config_t::buffer_len
 * Size of local sigma buffer.
 * 
 * @var pollpom_config_t::filename
 * Name of the file to write the output to
 */
typedef struct {
    double *accum;
    uint64_t bound;
    size_t num_chunks;
    size_t chunk_len;
    size_t buffer_len;
    char filename[128];
} pollpom_config_t;

/**
 * @brief 
 * 
 * @param cfg 
 */
void pollpom_kparent(pollpom_config_t *cfg);

