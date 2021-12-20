/**
 * @file pomyang.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2021-11-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#define FILENAME "counts.csv"

#include "sieve.h"
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include "properSumDiv.h"
#include "sn.h"

/**
 * uints64 are used ensure consistant behaviour about over unsigned wrap
 * We are safe as long as all values are below 18,446,744,073,709,551,615
 * The largest number that can possibly be produced is: max_bound ^ 4/3
 * This means that approx 281474976700000 (2.81x10^14) is the largest safe value for max_bound
 */

void writeBuffer(uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc);
void tabStats(uint8_t * characArr);
uint32_t writeImage(uint64_t image, uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc);
uint32_t chooseNumChunks(uint64_t max_bound);

uint32_t buffer_size;  // size of buffer alloced to each chunk
uint64_t max_bound;


int main(int argc, char **argv) {
    if (argc == 0) {
        printf("What tf are you trying to input\n");
        printf("./[max_bound] OR ./[max_bound][buffer_size]");
        exit(0);
    } else if (argc == 2) {
        max_bound = atol(argv[1]);
        buffer_size = 10000;
    } else if (argc == 3) {
        max_bound = atol(argv[1]);
        buffer_size = atol(argv[2]);
    } else {
        printf("What tf are you trying to input\n");
        printf("./[max_bound] OR ./[max_bound][buffer_size]");
        exit(0);
    }

    assert(max_bound % 2 == 0);
    assert(max_bound < 281474976700000);  // max safe value for max_bound
    assert(buffer_size > 0);

    // chooseNumChunks(max_bound);
    uint32_t numChunks = 1000;  // number of pieces we break the problem into for threading
    uint64_t chunk_size = max_bound/numChunks;  // amount of numbers over which one thread will run the PY algorithm

    // byte array that accumulates parent information
    // characFunc[i] holds the number of parents for (i + 1) * 2
    uint8_t *characFunc = malloc(max_bound/2);
    memset(characFunc, 0, max_bound/2);

    const uint64_t max_prime = (uint64_t) sqrt(2 * max_bound);
    const uint64_t prime_bound = (uint64_t) (1.25506 * (max_prime + 1) / log(max_prime + 1));
    uint32_t *primes = calloc(prime_bound + 1, sizeof(uint32_t));
    prime_sieve(max_prime, primes);

    // These vars handle the odd composite squares
    // these numbers are pushed into a buffer to be processed in the second part of the algor
    const uint64_t compositeBound = (uint64_t)pow(max_bound, .6666666666666666666);
    uint64_t *compSquares = calloc(compositeBound, sizeof(uint64_t));
    uint64_t compSquareCounter = 0;

    uint64_t *sigma = calloc (chunk_size / 2, sizeof(uint64_t));

    // this buffer holds preimages before they are written to disk in a batch
    uint64_t *imageBuffer = calloc(buffer_size, sizeof(uint64_t));
    uint32_t bufferInd;  // counter for the imageChunk buffer

    // enumerate_handle_t ennum_handle = init_enumerate_sn(max_bound, chunk_size);

    for (uint32_t i = 0; i < numChunks; i++) {
        bufferInd = 0;

        uint64_t m = (i * chunk_size);  // the number currently being processed by algo

        sum_of_divisors_odd(chunk_size, m, sigma, primes);  // Antons again, very fast sum of divisors
        // enumerated_range_t *range = enumerate_sn(ennum_handle);
        // for (size_t k = 0; k < range->len / 2; k++) {
        //     sigma[k] = range->s[2 *  k];
        // }

        for (uint64_t j = 0; j < chunk_size / 2; j++) {
            m = (i * chunk_size + 1) + (2 * j);  // sigma(m) = sigma[j]
            // assert(chunk_size == range->len);
            // printf("sigma[%ld] = %ld, range->s[%ld] = %ld\n", j, sigma[j], 2 * j, range->s[2 *  j]);
            // assert(sigma[j] == range->s[2 * j]);

            // catches evens values of sigma[j] and runs them through recurrance
            if (!(sigma[j] & 1)) {
                uint64_t  t = 3 * sigma[j] - (2 * m);
                while (t <= max_bound) {
                    bufferInd = writeImage(t, imageBuffer, bufferInd, characFunc);
                    t = (2 * t) + sigma[j];
                }
            }

            // catches primes are records them appropriately
            if (sigma[j] == m+1) {
                bufferInd = writeImage(m+1, imageBuffer, bufferInd, characFunc);
            } else if (m <= compositeBound && m != 1) {
                // We are counting the number of preimages for even numbers [1,max_bound] since s(1) = 0 it falls outside that range and should be excluded
                // catches odd composite values of m which are pushed to buffer for later processing
                compSquares[compSquareCounter] = m;
                compSquareCounter++;
            }
        }

        if (bufferInd > 0) {
            writeBuffer(imageBuffer, bufferInd, characFunc);
            bufferInd = 0;
        }
    }

    free(sigma);

    for (uint64_t  i = 0; i < compSquareCounter; i++) {
        uint64_t s_mSq = wheelDivSum(compSquares[i]*compSquares[i]);
        assert(compSquares[i] % 2 != 0);

        if (s_mSq <= max_bound) {
            bufferInd = writeImage(s_mSq, imageBuffer, bufferInd, characFunc);
        }
    }

    if (bufferInd > 0) {
        writeBuffer(imageBuffer, bufferInd, characFunc);
        bufferInd = 0;
    }

    tabStats(characFunc);
}

/**
 * @brief writes a image of s(n) to the chunks buffer and if the buffer is full it is written to the charac function
 * 
 * @param image s(n) = image found by the PY algorithm
 * @param imageBuffer buffer of images being held for a later write to the charac func
 * @param bufferInd index currently being filled by preimages
 * @param characFunc charac function where statistics for each number is accumulated
 * @return uint32_t number of images written
 */
uint32_t writeImage(uint64_t image, uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc) {
    imageBuffer[bufferInd] = image;
    bufferInd++;
    if (bufferInd == buffer_size) {
        writeBuffer(imageBuffer, bufferInd, characFunc);
        bufferInd = 0;
    }
    return bufferInd;
}

/**
 * @brief This function writes a buffer of preimages to the characFile
 * 
 * @param imageBuffer buffer full of images to be acc'ed in characFunc
 * @param bufferInd index currently being filled by preimages
 * @param characFunc charac function where statistics for each number is accumulated
 */
void writeBuffer(uint64_t * imageBuffer, uint32_t bufferInd, uint8_t * characFunc) {
    for (uint32_t i = 0; i < bufferInd; i++) {
        characFunc[(imageBuffer[i]/2)-1]++;
    }
}

/**
 * @brief This functions processes the characArray and tabulates statisics about k-parent numbers characArr must be max_bound/2 in size
 * 
 * @param characArr char array where the number of preimages for each evenb number is written
 */
void tabStats(uint8_t * characArr) {
    FILE * fp = fopen(FILENAME, "a");
    uint64_t accPreimages[256] = {0};

    for (uint64_t i = 0;  i < max_bound/2; i++) {
        accPreimages[characArr[i]]++;
    }

    fprintf(fp, "%lu", max_bound);
    for (uint32_t j = 0;  j < 255; j++) {
       fprintf(fp, ",%lu",  accPreimages[j]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

/**
 * @brief   When threading in this manner it is necessary to choose a number of chunks to 
 *          split the problem into such that each thread can work through the problem a chunk at a time
 *          chooses the first divisor greater than the sqrt of max bound to be the number of chunks
 * @param max_bound bound
 * @return uint32_t number of chunks to split problem into
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
