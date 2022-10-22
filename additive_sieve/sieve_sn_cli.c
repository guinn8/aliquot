/**
 * @file sieve_sn_cli.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief
 * @date 2021-11-28
 *
 * @copyright Copyright (c) 2021
 *
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "../additive_sieve/sieve_sn.h"

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("./sn <max> <blk_len>\n");
        assert(0);
    }

    size_t max = strtol(argv[1], NULL, 10);
    size_t blk_len = strtol(argv[2], NULL, 10);

    enumerate_handle_t sts = init_enumerate_sn(max, blk_len);
    enumerated_range_t *range;
    while (NULL != (range = enumerate_sn(sts))) {
        for (size_t i = 0; i < range->len; i++) {
            printf("%ld\n", range->s[i] - (i + range->base));
        }
    }

    destroy_enumerate_status(sts);
    return 0;
}
