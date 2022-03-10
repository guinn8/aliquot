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

    size_t start = 100;
    size_t end = 1000000000;

    while (start <= end) {
        for (size_t i = 1; i <= 9; i++) {
            pollpom_config_t cfg = {
                .bound = i * start,
                .chunk_len = start / 10,
                .filename = FILENAME,
            };
            pollpom_kparent(&cfg);
        }
        start *= 10;
    }

    pollpom_config_t cfg = {
        .bound = end,
        .chunk_len = 100000,
    };
    pollpom_kparent(&cfg);


    return 0;
}
