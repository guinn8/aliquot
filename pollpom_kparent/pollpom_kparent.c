/**
 * @file pollpom_kparent.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief
 * This program implements an generaliztion of Conj. 1.4 of Pollack/Pomerance "Some problems of Erdos on the Sum of Divisors Function"
 * Instead of estimating the natural density of only aliqout orphans this program also estimates the density of k-parent aliquot numbers
 * n is a k-parent aliqout number iff there are k distinct natural numbers m st s(m) = n
 * An aliquot orphan is a 0-parent aliquot number
 * let delta-k be the estimated density of k-parent aliquot numbers and s(n) be the sum-of-proper-divisors function
 * delta-k = 1/log(bound) * sum(forall a <= bound)( (a^(k-1) * e^(-a/s(a)) / k! * s(a)^k) )
 * @date 2022-02-24
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "../pollpom_kparent/pollpom_kparent.h"

#include <alloca.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../factor/factor.h"
#include "../pomyang_kparent/inc/math_macros.h"

#define MAX_PARENT 16

static uint64_t factorial(uint64_t n);
static void accumulator(pollpom_config_t *cfg, uint64_t start_a, uint32_t k, uint64_t *propSumDiv, double *acc);
static void print_to_file(pollpom_config_t *cfg, float runtime);


/* See header for documentation */
void pollpom_kparent(pollpom_config_t *cfg) {
    assert(DIVIDES(cfg->chunk_len, cfg->bound));
    cfg->num_chunks = cfg->bound / cfg->chunk_len;
    cfg->accum = calloc(MAX_PARENT, sizeof(double));
    cfg->buffer_len = cfg->chunk_len / 2;

    clock_t start = clock();
#pragma omp parallel
    {
        uint64_t *sigma = malloc(cfg->buffer_len * sizeof(uint64_t));

#pragma omp for
        for (size_t i = 0; i < cfg->num_chunks; i++) {
            uint64_t m = (cfg->chunk_len * i) + 2;  // first even number in chunk

            for (size_t j = 0; j < cfg->buffer_len; j++) {
                sigma[j] = factor_s(m + 2 * j);
            }

            for (size_t g = 0; g < 16; g++) {
                accumulator(cfg, m, g, sigma, cfg->accum);
            }
        }
        free(sigma);
    }

    for (int i = 0; i < 16; i++) {
        cfg->accum[i] *= 1 / (double)(log(cfg->bound) * factorial(i));
        printf("delta%d\t= %f %% \n", i, 100 * cfg->accum[i]);
    }
    printf("\n");

    print_to_file(cfg, clock() - start);
    free(cfg->accum);

}

static void accumulator(pollpom_config_t *cfg, uint64_t start_a, uint32_t k, uint64_t *propSumDiv, double *acc) {
    double numer;
    double denom;
    double a;
    uint64_t s_a;

    for (uint64_t i = 0; i < cfg->buffer_len; i++) {
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

#pragma omp atomic
        acc[k] += numer / denom;
    }
}

static void print_to_file(pollpom_config_t *cfg, float runtime) {
    const size_t max_line_len = 10000;
    char *header_line = alloca(max_line_len * sizeof(char));
    memset(header_line, 0, max_line_len * sizeof(char));
    if (0 != access(cfg->filename, F_OK)) {
        snprintf(header_line, max_line_len, "timestamp, bound, chunk_size, timing");
        for (size_t i = 0; i < MAX_PARENT; i++) {
            snprintf(header_line + strlen(header_line), max_line_len - strlen(header_line), ", %ld", i);
        }
    }

    FILE *fp = fopen(cfg->filename, "a");
    fprintf(fp, "%s\n", header_line);
    fprintf(fp, "%ld, ", time(NULL));
    fprintf(fp, "%ld, ", cfg->bound);
    fprintf(fp, "%ld, ", cfg->buffer_len);
    fprintf(fp, "%.2f, ", runtime);

    for (size_t i = 0; i < MAX_PARENT; i++) {
        fprintf(fp, "%0.16f, ", cfg->accum[i]);
    }
    fclose(fp);
}

/**
 * @brief Computes factorial
 * 
 * @param n input number
 * @return uint64_t n!
 */
static uint64_t factorial(uint64_t n) {
    if (n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}
