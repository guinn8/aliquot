/**
 * @file brute_force_preimages.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-02-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _POM_YANG_INC_BRUTE_FORCE_PREIMAGES_H_
#define _POM_YANG_INC_BRUTE_FORCE_PREIMAGES_H_

#include <stdlib.h>
#include <stdint.h>

uint8_t *brute_force_preimages(size_t bound);
uint64_t *brute_force_preimage_counts(size_t bound);

#endif  // _POM_YANG_INC_BRUTE_FORCE_PREIMAGES_H_
