/**
 * @file sieve_sn.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2021-11-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef FAST_ENNUM_SN_SN_H_
#define FAST_ENNUM_SN_SN_H_

#include <stdint.h>

typedef struct {
    uint64_t *s;
    uint64_t len;
    uint64_t base;
} enumerated_range_t;

struct enumerate_status_;
typedef struct enumerate_status_* enumerate_handle_t;

enumerate_handle_t init_enumerate_sn(size_t max, size_t blk_len);
enumerated_range_t *enumerate_sn(enumerate_handle_t sts);
void destroy_enumerate_status(enumerate_handle_t sts);

#endif  // FAST_ENNUM_SN_SN_H_
