/**
 * @file sumdiv_util.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief   Uses factoring wheel (https://en.wikipedia.org/wiki/Wheel_factorization) to compute the sum of divisors.
 *          Not well tested but mostly used a backstop test for sigma sieves.
 * @date 2022-02-17
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#include "../inc/sumdiv_util.h"

/** @brief Number wheel used to check divisors.*/
const int initWheel[3] = {2, 3, 5};

/** @brief ?*/
const int inc[8] = {4, 2, 4, 2, 4, 6, 2, 6};

/**
 * @brief compute sigma of n (sum-of-divisors)
 * 
 * @param n input
 * @return uint64_t sigma(n)
 */
uint64_t sumdiv_sigma(uint64_t n) {
    if (n == 0) {
        return 0;
    }

    uint64_t curr_sum = 1;
    uint64_t curr_term = 1;
    uint64_t res = 1;
    for (int i = 0; i < 3; i++) {
        if (n % initWheel[i] == 0) {
            do {
                n = n / initWheel[i];
                curr_term *= initWheel[i];
                curr_sum += curr_term;
            } while (n % initWheel[i] == 0);

            res *= curr_sum;
            curr_sum = curr_term = 1;
        }
    }

    uint64_t k = 7;
    int i = 0;
    while (k <= n) {
        if (n % k == 0) {
            do {
                n = n / k;
                curr_term *= k;
                curr_sum += curr_term;
            } while (n % k == 0);
            res *= curr_sum;
            curr_sum = curr_term = 1;
        } else {
            k += inc[i];
            if (i < 7)
                i++;
            else
                i = 0;
        }
    }
    res *= curr_sum;
    return res;
}

/**
 * @brief compute s of n (sum-of-proper-divisors)
 * 
 * @param n input
 * @return uint64_t s(n)
 */
uint64_t sumdiv_s(uint64_t n) {
    return sumdiv_sigma(n) - n;
}
