/**
 * @file pollpom_kparent_cli.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-03-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <string.h>

#include "../pollpom_kparent/pollpom_kparent.h"

/** @brief Name of output file */
#define FILENAME "densities.csv"

/** @brief Interface to Pollack and Pomerance kParent generalization */
int main(int argc, char **argv) {
    if (argc != 3) {
        printf("USAGE: ./[bound][chunk_len]\n");
        exit(0);
    }

    pollpom_config_t cfg = {
        .bound = atol(argv[1]),
        .chunk_len = atol(argv[2]),
    };

    strncpy(cfg.filename, FILENAME, sizeof(cfg.filename));

    pollpom_kparent(&cfg);
}
