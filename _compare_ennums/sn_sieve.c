/**
 * @file sn_sieve.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2021-11-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "sieve.h"

uint32_t chooseNumChunks(uint64_t max_bound);
// #define ODD
#ifdef ODD
int main(int argc, char **argv) {
    if (2 != argc) {
        printf("./sn_sieve <bound>");
        assert(0);
    }

    size_t bound = strtol(argv[1], NULL, 10);
    assert(0 == bound % 2);

    const uint64_t max_prime = sqrt(2 * bound);
    const uint64_t prime_bound = ceil(1.25506 * (max_prime + 1) / log(max_prime + 1));
    uint32_t * primes = calloc(prime_bound + 1, sizeof(uint32_t));
    prime_sieve(max_prime, primes);

    uint32_t numChunks = 100;//chooseNumChunks(bound);
    uint64_t chunk_size = bound / numChunks;
    uint64_t * sigma = calloc(chunk_size / 2, sizeof(uint64_t));

    
     for (uint32_t i = 0; i < numChunks; i++) {
        uint64_t m = (i * chunk_size);
        sum_of_divisors_odd(chunk_size, m, sigma, primes);
        for (size_t j = 0; j < chunk_size / 2; j++) {
            printf("%ld\n", sigma[j]);
        }
    }
}
#else
int main(int argc, char **argv) {
    if (2 != argc) {
        printf("./sn_sieve <bound>");
        assert(0);
    }

    size_t bound = strtol(argv[1], NULL, 10);
    assert(0 == bound % 2);

    const uint64_t max_prime = sqrt(2 * bound);
    const uint64_t prime_bound = ceil(1.25506 * (max_prime + 1) / log(max_prime + 1));
    uint32_t * primes = calloc(prime_bound + 1, sizeof(uint32_t));
    prime_sieve(max_prime, primes);

    uint32_t numChunks = 100;  // chooseNumChunks(bound);
    uint64_t chunk_size = bound / numChunks;
    uint64_t * sigma = calloc(chunk_size, sizeof(uint64_t));

    for (uint32_t i = 0; i < numChunks; i++) {
        uint64_t m = (i * chunk_size);
        sum_of_divisors(chunk_size, m, sigma, primes);
        // for (size_t j = 0; j < chunk_size; j++) {
        //     printf("%ld\n", sigma[j]);
        // }
    }
}
#endif // ALL


/**
 * @brief  When threading in this manner it is needed to choose a number of chunks to
 * split the problem into such that each thread can work through the problem a chunk at a time
 * chooses the first divisor greater than the sqrt of max bound to be the number of chunks
 * 
 * @param max_bound bound to be split into chunks
 * @return uint32_t amount of chucks
 */
uint32_t chooseNumChunks(uint64_t max_bound) {
    uint32_t chunks;
    uint32_t root = sqrt(max_bound);
    for (uint32_t i=1; i <= max_bound; i++) {
        if (0 == max_bound % i) {
            if (i >= root && (max_bound / i) % 2 == 0) {
                chunks = i;
                return chunks;
            }
        }
    }

    printf("Could not find a suitable divisor i of maxbound st maxbound/i is even");
    exit(EXIT_FAILURE);
}
