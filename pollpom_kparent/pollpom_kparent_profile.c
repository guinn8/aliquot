/**
 * @file pollpom_kparent_profile.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief Run a number of inputs into the Pollack and Pomerance model generalization
 * @date 2022-03-07
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <stdio.h>
#include <string.h>
#include "../pollpom_kparent/pollpom_kparent.h"

#define FILENAME "profile.csv"

int main(int argc, char const *argv[]) {
    (void)argc;
    (void)argv;

    size_t start = 16384;
    size_t end = 1099511627776;

    size_t current = start;
    while (current <= end) {
        pollpom_config_t cfg = {
            .bound = current,
            .chunk_len = start / 2,
            .filename = FILENAME,
        };
        pollpom_kparent(&cfg);
        current *= 2;
    }


    return 0;
}
