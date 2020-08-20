/*
 * sieve.h
 *
 *  Created on: Apr 18, 2014
 *      Author: antonmosunov
 */

#ifndef SIEVE_H_
#define SIEVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 2^20=1048576
#define FAC_TOTAL 1048576
#define WITH_INDICES 1

#ifndef ulong
#define ulong unsigned long 
#endif

#ifndef uint
#define uint unsigned int
#endif

void prime_sieve(const uint max_prime, uint * primes);

void regular_sieve(const uint max_prime, const ulong blocksize, uint ** factors, const uint * primes, const int flags);

void segmented_sieve(const uint max_prime, const ulong blocksize, const ulong l, uint ** factors, const uint * primes, const int flags);

#if BITSIZE == 32
void sum_of_divisors(const ulong blocksize, const ulong l, uint * sigma, const uint * primes, ulong * max_val, ulong * max_sigma);

void sum_of_divisors_odd(const ulong blocksize, const ulong l, uint * sigma, const uint * primes/*, ulong * max_val, ulong * max_sigma*/);

// used to tabulate s(m^2) in Pomerance-Yang algorithm
void sum_of_divisors_odd2(const ulong blocksize, const ulong l, uint * sigma, const uint * primes);
#else
void sum_of_divisors(const ulong blocksize, const ulong l, ulong * sigma, const uint * primes);

void sum_of_divisors_odd(const ulong blocksize, const ulong l, ulong * sigma, const uint * primes/*, ulong * max_val, ulong * max_sigma*/);

// used to tabulate s(m^2) in Pomerance-Yang algorithm
void sum_of_divisors_odd2(const ulong blocksize, const ulong l, ulong * sigma, const uint * primes);
#endif



#endif /* SIEVE_H_ */
