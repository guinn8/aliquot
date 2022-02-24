/**
 * @file sumdiv_util.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief api for sumdiv functions
 * @date 2022-02-17
 *
 * @copyright Public Domain (Please credit me; if you find this code useful I would love to hear about your work!)
 *
 */

#ifndef _POM_YANG_INC_SUMDIV_H_
#define _POM_YANG_INC_SUMDIV_H_

#include <stdint.h>

uint64_t sumdiv_s(uint64_t n);
uint64_t sumdiv_sigma(uint64_t n);

#endif  // _POM_YANG_INC_SUMDIV_H_
