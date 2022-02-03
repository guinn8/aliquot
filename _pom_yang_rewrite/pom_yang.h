/**
 * @file pom_yang.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-02-02
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdlib.h>
#include <stdbool.h>

#include "../PackedArray/PackedArray.h"


typedef struct {
    size_t preimage_count_bits;
    size_t bound;
    size_t seg_len;
    size_t writebuf_len;
    size_t num_threads;
    bool just_config;
} PomYang_config;

PackedArray *Pomerance_Yang_aliquot(const PomYang_config *cfg);
uint64_t *count_kparent_aliquot(const PomYang_config *cfg);
