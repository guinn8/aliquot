/**
 * @file main.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief use https://www2.math.upenn.edu/~deturck/m170/wk3/lecture/sumdiv.html in a sieve
 * @date 2022-03-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
// https://oeis.org/A000040/list
static const uint32_t primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
                                  61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
                                  131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
                                  193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
                                  263, 269, 271};

static const uint32_t sigma[] = {1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24, 31, 18, 39,
                                 20, 42, 32, 36, 24, 60, 31, 42, 40, 56, 30, 72, 32, 63, 48, 54,
                                 48, 91, 38, 60, 56, 90, 42, 96, 44, 84, 78, 72, 48, 124, 57, 93,
                                 72, 98, 54, 120, 72, 120, 80, 90, 60, 168, 62, 96, 104, 127, 84,
                                 144, 68, 126, 96, 144};

int main() {
    uint64_t bound = 1000000;
    uint32_t s = sqrt(bound);
    for (size_t i = 0; i < 58; i++) {
        uint64_t p = primes[i];
        for (uint64_t j = p; j <= sqrt(bound) / p; j *= p) {
            assert(j > 0);
            usleep(10000);
            uint64_t x = 0;
            while (x <= s - j) {
                x += j;
                printf("%lu\t", x);
            }
            printf("\n");
        }
        printf("\n\n");
    }

    // while ((pow * p) < UINT64_MAX) {
    // 	pow *= p;
    //     printf("%d\n", pow);
    // }

    return (0);
}
