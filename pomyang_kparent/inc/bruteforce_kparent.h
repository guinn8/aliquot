/**
 * @file bruteforce_kparent.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief api to brute force number of preimages
 * @date 2022-02-11
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#ifndef _POM_YANG_INC_BRUTEFORCE_KPARENT_H_
#define _POM_YANG_INC_BRUTEFORCE_KPARENT_H_

#include <stdint.h>
#include <stdlib.h>

/**
 * @brief compute number of preimages for even number upto bound
 *
 * @param bound
 * @return uint8_t* buffer sized bound/2 counting counts of preimages, free'ed by caller.
 */
uint8_t *bf_kparent(size_t bound);

/**
 * @brief count number of even kparent numbers upto bound
 *
 * @param bound
 * @return uint64_t* buffer of size UINT8_MAX containing preimage counts, free'ed by caller.
 */
uint64_t *bf_kparent_counts(size_t bound);

#endif  // _POM_YANG_INC_BRUTEFORCE_KPARENT_H_
