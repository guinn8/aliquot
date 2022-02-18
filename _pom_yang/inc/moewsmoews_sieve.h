/**
 * @file sieve_rewrite.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief api to the Moews and Moews sieving algorithm
 * @date Originally Apr 18, 2014, modified 2021-12-27
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#ifndef _POM_YANG_INC_MOEWS_MOEWS_SIEVE_H_
#define _POM_YANG_INC_MOEWS_MOEWS_SIEVE_H_

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

/**
 * @brief configuration and computed values needed to operate sieve, only 1 needed for multiple threads
 *
 * @param bound sieve all odd sigma upto to bound
 * @param seg_len size of range to sieve in a step
 * @param buf_len half seg_len as we are only sieving evens
 * @param primes_len length of primes buffer
 * @param primes buffer of primes needed to operate sieve
 */
typedef struct {
    size_t bound;
    size_t seg_len;
    size_t buf_len;
    size_t primes_len;
    uint32_t *primes;
} sieve_config_t;

/**
 * @brief worker which holds state about a single thread running the sieve
 *
 * @param cfg initalized sieve_config_t
 * @param squared was the last sieve block sigma(m * m)
 * @param seg_start start of the segment to be sieved
 * @param sigma_buf buffer to be filled will odd sigma in range (seg_start, seg_start + seg_len)
 * @param q_buf buffer used interally by the sieve
 * @param is_prime_buf <after running moews_sieve_odd_standard> this buffer store whether odd numbers...
 *                     in range (seg_start, seg_start + seg_len) are prime
 */
typedef struct {
    const sieve_config_t *cfg;
    size_t seg_start;
    uint64_t *sigma_buf;
    uint64_t *q_buf;
    ssize_t last_sieve_standard;
    bool squared;
    bool *is_prime_buf;
} sieve_worker_t;

/**
 * @brief runs sieve for odd sigma(m)
 *
 * @param worker see struct definiton
 * @param seg_start start of range to sieve
 */
void moews_sieve_odd_standard(sieve_worker_t *worker, uint64_t seg_start);

/**
 * @brief takes odd m in range and runs sieve for sigma(m * m)
 *
 * @param worker see struct definiton
 * @param seg_start start of range to sieve
 */
void moews_sieve_odd_squared(sieve_worker_t *worker, uint64_t seg_start);

/**
 * @brief lookup sigma(m) from m in a sieved block
 *
 * @param worker see definition, has run moews_sieve_odd_standard last
 * @param m number to lookup
 * @return uint64_t sigma(m)
 */
uint64_t moews_lookup_sigma_m(const sieve_worker_t *worker, uint64_t m);

/**
 * @brief lookup sigma(m * m) from m in a sieved block
 *
 * @param worker see definition, has run moews_sieve_odd_squared last
 * @param m number to lookup
 * @return uint64_t sigma(m * m)
 */
uint64_t moews_lookup_sigma_m_squared(const sieve_worker_t *worker, uint64_t m);

/**
 * @brief creates a sieve_config_t which is shared by a pool of workers
 *
 * @param bound sieving all odd sigma upto to bound
 * @param seg_len size of range to sieve in a step
 * @return sieve_config_t* config struct which is free'd by calling moews_destroy_sieve
 */
sieve_config_t *moews_init_sieve(size_t bound, size_t seg_len);

/**
 * @brief free's memory associated with sieve_config_t
 *
 * @param cfg to be free'd
 */
void moews_destroy_sieve(sieve_config_t *cfg);

/**
 * @brief creates a worker object which can be used to run the sieve
 *
 * @param cfg initalized sieve configuration struct
 * @return sieve_worker_t*
 */
sieve_worker_t *moews_init_worker(const sieve_config_t *cfg);

/**
 * @brief free's memory associated with sieve_worker_t
 *
 * @param worker to be free'd
 */
void destroy_worker(sieve_worker_t *worker);

/**
 * @brief checks if odd numbers in range (seg_start, seg_start + seg_len) is prime
 *
 * @param worker that has previously run moews_sieve_odd_standard
 * @param m check if this odd number in range (seg_start, seg_start + seg_len) is prime
 * @return bool is prime
 */
bool moews_check_prime(const sieve_worker_t *worker, uint64_t m);

/**
 * @brief estimates how much heap the sieving buffers will require
 *
 * @param bound sieve all odd sigma upto to bound
 * @param seg_len size of range to sieve in a step
 * @param num_workers number of threads working
 * @return size_t approx bytes to be consumed
 */
size_t moews_estimate_heap_usage(size_t bound, size_t seg_len, size_t num_workers);

#endif  // _POM_YANG_INC_MOEWS_MOEWS_SIEVE_H_
