/**
 * @file sn.c
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2021-11-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

void accum_array_multiples(uint64_t *sn_buf, uint64_t *multiple_buf, size_t sqrt_max, size_t current_block, size_t block_size, size_t max_divisor);
void iterate_multiples(uint64_t *multiple_buf, size_t sqrt_max, size_t current_block, size_t max_divisor);
void compute_sn_buffer_vec_add_sqrt(size_t max, size_t block_size);
size_t find_max_divisor(size_t previous_max_divisor, uint64_t *multiple_buf, size_t current_block);

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("./sn <max> <block_size>\n");
        assert(0);
    }

    size_t max = strtol(argv[1], NULL, 10);
    size_t block_size = strtol(argv[2], NULL, 10);
    compute_sn_buffer_vec_add_sqrt(max, block_size);

    return 0;
}

void compute_sn_buffer_vec_add_sqrt(size_t max, size_t block_size) {
    assert(0 == max % 2);
    assert(0 == max % block_size);

    const size_t sqrt_max = sqrt(max);

    uint64_t *sn_buf = calloc(block_size, sizeof(uint64_t));
    for (size_t i = 1; i < block_size; i++) {
        sn_buf[i] = 1;
    }
    sn_buf[0] = 0;

    uint64_t *multiple_buf = calloc(sqrt_max + 1, sizeof(uint64_t));
    for (size_t i = 0; i <= sqrt_max; i++) {
        multiple_buf[i] = 2 * i;
    }

    for (size_t current_block = block_size; current_block <= max; current_block += block_size) {
        // printf("\nSTARTED BLOCK %ld\n", current_block);
        size_t max_divisor = sqrt_max;
        accum_array_multiples(sn_buf, multiple_buf, sqrt_max, current_block, block_size, max_divisor);

        while (1 < (max_divisor = find_max_divisor(max_divisor, multiple_buf, current_block))) {
            iterate_multiples(multiple_buf, sqrt_max, current_block, max_divisor);
            accum_array_multiples(sn_buf, multiple_buf, sqrt_max, current_block, block_size, max_divisor);
        }

        for (size_t i = 0; i < block_size; i++) {
            // printf("s(%ld) = %ld\n", i  + (current_block - block_size) + 1, sn_buf[i]);
            size_t n = i + (current_block - block_size) + 1;
            #define ODD
            #ifdef ODD
            if (0 != n % 2)
            #endif  // ODD
            {
                printf("%ld\n", sn_buf[i] + n);
            }
            sn_buf[i] = 1;
        }
        // printf("\nFINISHED BLOCK %ld\n", current_block);
    }

    free(sn_buf);
    free(multiple_buf);
}

size_t find_max_divisor(size_t previous_max_divisor, uint64_t *multiple_buf, size_t current_block) {
    for (size_t l = previous_max_divisor; l >= 0; l--) {
        if (multiple_buf[l] <= current_block) {
            // printf("previous_max_divisor = %ld\tnew_max_divisor = %ld\n", previous_max_divisor, l);
            return l;
        }
    }
}

void accum_array_multiples(uint64_t *sn_buf, uint64_t *multiple_buf, size_t sqrt_max, size_t current_block, size_t block_size, size_t max_divisor) {
    const size_t array_offset = (current_block - block_size) + 1;
    for (size_t i = 2; i <= max_divisor; i++) {
        if (multiple_buf[i] <= current_block) {
            size_t ele_offset = multiple_buf[i] - array_offset;
            sn_buf[ele_offset] += i;
            // printf("%ld divides %ld", i, multiple_buf[i]);
            if (multiple_buf[i]/i > sqrt_max) {
                // printf(", %ld divides %ld", multiple_buf[i]/i, multiple_buf[i]);
                sn_buf[ele_offset] += multiple_buf[i]/i;
            }
            // printf("\n");
        }
    }
    // printf("\n");
}

void iterate_multiples(uint64_t *multiple_buf, size_t sqrt_max, size_t current_block, size_t max_divisor) {
    size_t bound = (current_block < max_divisor) ? current_block : max_divisor;
    for (size_t i = 2; i <= bound; i++) {
        if (multiple_buf[i] <= current_block) {
            multiple_buf[i] += i;
            // printf("%ld divides %ld\n", i, multiple_buf[i]);
        }
    }
    // printf("\n");
}
