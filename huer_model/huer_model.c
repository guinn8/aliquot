/**
 * @file huer_model.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-02-24
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define UPPERPARENTS 16

int numChunks = 10;
int chunk_size;
int buffer_size;

uint64_t s(uint64_t n);
uint64_t factorial(uint64_t n);
double accumulator(uint64_t start_a, char k, uint64_t *propSumDiv, long double *acc);

// This program implements an generaliztion of Conj. 1.4 of Pollack/Pomerance "Some problems of Erdos on the Sum of Divisors Function"
// Instead of estimating the natural density of only aliqout orphans this program also estimates the density of k-parent aliquot numbers
// n is a k-parent aliqout number iff there are k distinct natural numbers m st s(m) = n
// An aliquot orphan is a 0-parent aliquot number
// let delta-k be the estimated density of k-parent aliquot numbers and s(n) be the sum-of-proper-divisors function
// delta-k = 1/log(max_bound) * sum(forall a <= max_bound)( (a^(k-1) * e^(-a/s(a)) / k! * s(a)^k) )
int main(int argc, char *argv[]) {
    uint64_t max_bound;

    long double acc[UPPERPARENTS] = {0};

    if (argc < 1) {
        printf("./[max_bound]");
        exit(0);
    }

    max_bound = atol(argv[1]);
    chunk_size = max_bound / numChunks;
    buffer_size = chunk_size / 2;

    // loop through the chunks for multi-threading
    for (int i = 0; i < numChunks; i++) {
        uint64_t m = (chunk_size * i) + 2;  // first even number in chunk
        uint64_t *sigma = malloc(buffer_size * sizeof(uint64_t));

        for (int j = 0; j < buffer_size; j++) {
            sigma[j] = sumdiv_s(m + 2 * j);
        }

        for (int g = 0; g < UPPERPARENTS; g++) {
            accumulator(m, g, sigma, acc);
        }
    }

    for (int i = 0; i < UPPERPARENTS; i++) {
        acc[i] *= 1 / (log(max_bound) * factorial(i));
        printf("delta %d = %Lf\n", i, acc[i]);
    }
}

uint64_t s(uint64_t n) {
    return 0;  // fmpz_get_ui(0);
}

uint64_t factorial(uint64_t n) {
    if (n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}

double accumulator(uint64_t start_a, char k, uint64_t *propSumDiv, long double *acc) {
    double numer;
    double denom;
    double a;
    uint64_t s_a;

    for (uint64_t i = 0; i < buffer_size; i++) {
        a = (2 * i) + start_a;
        s_a = propSumDiv[i];

        // This needs to be treated casewise because computers dont understand neg. exp. correctly
        if (k == 0) {
            numer = exp(-a / s_a);
            denom = a;
        } else {
            numer = pow(a, k - 1) * exp(-a / s_a);
            denom = pow(s_a, k);
        }

        acc[k] += numer / denom;
    }
}
