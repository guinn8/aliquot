/*
 * sieve.c
 *
 *  Created on: Apr 18, 2014
 *      Author: antonmosunov
 *  USED WITH PERMISSION
 */

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

#include "sieve.h"

void prime_sieve(const uint max_prime, uint * primes)
{
	uint i, j;

	primes[0] = 0;

	char is_prime[max_prime + 1];
	memset(is_prime, 1, max_prime + 1);

	for (i = 2; i <= max_prime; i++)
	{
		if (is_prime[i])
		{
			primes[++primes[0]] = i;

			for (j = (i << 1); j <= max_prime; j += i)
			{
				is_prime[j] = 0;
			}
		}
	}

	primes[primes[0] + 1] = 0x7FFFFFFF;
}


void regular_sieve(const uint max_prime, const ulong blocksize, uint ** factors, const uint * primes, const int flags)
{
	ulong i, j;

	uint p, * f;

	factors[0][0] = 1;
	factors[0][1] = 0;
	factors[1][0] = 1;
	factors[1][1] = 1;

	for (i = 0; i < blocksize; i++)
	{
		factors[i][0] = 0;
	}

	for (i = 1, p = primes[i]; p < max_prime; p = primes[++i])
	{
		for (j = p; j < blocksize; j += p)
		{
			f = factors[j];
			f[++f[0]] = (flags & WITH_INDICES) ? i : p;
		}
	}
}


void segmented_sieve(uint max_prime, ulong blocksize, ulong l, uint ** factors, const uint * primes, const int flags)
{
	if (l == 0)
	{
		regular_sieve(max_prime, blocksize, factors, primes, flags);
		return;
	}

	uint i, j, k, p, offset, * f;

	for (k = 0; k < blocksize; k++)
	{
		factors[k][0] = 0;
	}

	for (i = 1, p = primes[i]; p < max_prime; p = primes[++i])
	{
		offset = (p - (l % p)) % p;

		for (j = offset; j < blocksize; j += p)
		{
			f = factors[j];
			f[++f[0]] = (flags & WITH_INDICES) ? i : p;
		}
	}
}

//Not sure if this is currently functional
#if BITSIZE == 32
void sum_of_divisors(const ulong blocksize, const ulong l, uint * sigma, const uint * primes, ulong * max_val, ulong * max_sigma)
#else
void sum_of_divisors(const ulong blocksize, const ulong l, ulong * sigma, const uint * primes)
#endif
{
	const ulong max_prime = (ulong) sqrt(l + blocksize);

	ulong h, i, j, k, p, p_pow, step, offset;

	ulong * numbers = (ulong *) malloc(sizeof(ulong) * blocksize);

	for (i = 0, j = l; i < blocksize; i++, j++)
	{
		sigma[i] = 1;
		numbers[i] = j;
	}

	if (l == 0)
	{
		sigma[0] = 0;
	}

	for (i = 1, p = primes[i]; p <= max_prime; p = primes[++i])
	{
		for (p_pow = p, step = p * p; p_pow <= (l + blocksize); p_pow *= p, step *= p)
		{
			offset = (p_pow - (l % p_pow)) % p_pow;

			if (offset > blocksize)
			{
				break;
			}

			//printf("p_pow=%lu, offset=%lu, step=%lu\t", p_pow, offset, step);
			for (k = 0, h = offset; k < p; k++, h += p_pow)
			{
				if ((l + h) % step != 0)
				{
					for (j = h; j < blocksize; j += (step << 1))
					{
						//printf(", %lu", l + j);
						numbers[j] /= p_pow;
						sigma[j] *= ((step - 1) / (p - 1));
					}
				}
			}
			//printf("\n");
		}
	}

	if (l == 0)
	{
		i = j = 1;
	}
	else
	{
		i = 0;
		j = l;
	}

	for (; i < blocksize; i++, j++)
	{
		if (numbers[i] > 1)
		{
			sigma[i] *= (numbers[i] + 1);
		}

		

		#if BITSIZE == 32
		printf("sigma(%lu)=%u\n", j, sigma[i]);
		#else
		printf("sigma(%lu)=%lu\n", j, sigma[i]);
		#endif
	}

	free(numbers);
}

#if BITSIZE == 32
void sum_of_divisors_odd(const ulong blocksize, const ulong l, uint * sigma, const uint * primes/*, ulong * max_val, ulong * max_sigma*/)
#else
void sum_of_divisors_odd(const ulong blocksize, const ulong l, ulong * sigma, const uint * primes/*, ulong * max_val, ulong * max_sigma*/)
#endif
{
	if (((blocksize & 1) == 1) || ((l & 1) == 1))
	{
		perror("The parameters blocksize and l must be even.\n");
		exit(1);
	}

	const ulong max_prime = (ulong) sqrt(l + blocksize);

	ulong h, i, j, k, p, s, p_pow, step;

	ulong * numbers = (ulong *) malloc(sizeof(ulong) * (blocksize >> 1));

	for (i = 0, j = l + 1; i < (blocksize >> 1); i++, j += 2)
	{
		sigma[i] = 1;
		numbers[i] = j;
	}

	if (l == 0)
	{
		sigma[0] = 1;
	}

	ulong offset[100];

	for (i = 2, p = primes[i]; p <= max_prime; p = primes[++i])
	{
		k = 0;
		offset[k] = (p - (l % p)) % p;

		if ((offset[k] & 1) == 0)
		{
			offset[k] += p;
		}

		for (p_pow = p, step = p * p; p_pow <= (l + blocksize); p_pow *= p, step *= p)
		{
			if (offset[k] > blocksize)
			{
				break;
			}

			#ifdef DEBUG
			printf("p_pow=%lu, offset=%lu, step=%lu\t", p_pow, offset[k], step);
			#endif

			offset[++k] = (step - (l % step)) % step;

			if ((offset[k] & 1) == 0)
			{
				offset[k] += step;
			}

			for (s = 0, h = offset[k-1]; s < p; s++, h += (p_pow << 1))
			{
				if (h != offset[k])
				{
					for (j = h; j < blocksize; j += (step << 1))
					{
						#ifdef DEBUG
						printf(", %lu (%lu)", l + j, (j-1)>>1);
						#endif

						numbers[(j-1) >> 1] /= p_pow;
						sigma[(j-1) >> 1] *= ((step - 1) / (p - 1));
					}
				}
			}

			#ifdef DEBUG
			printf("\n");
			#endif
		}
	}

	i = 0;
	j = l + 1;

	for (; i < (blocksize >> 1); i++, j += 2)
	{
		if (numbers[i] > 1)
		{
			sigma[i] *= (numbers[i] + 1);
		}
		
		#ifdef DEBUG
		printf("sigma(%lu)=%lu\n", j, (ulong) sigma[i]);
		#endif

		#ifdef WITH_PARI
		GEN g = stoi(j);
		GEN v = gsumdivk(g, 1);

		if (sigma[i] != gtolong(v))
		d ..cd
			char err[240];
			sprintf(err, "Error, sigma(%lu)=%lu, but we computed %lu.\n", j, gtolong(v), (ulong) sigma[i]);
			perror(err);
			cgiv(v);
			cgiv(g);
			pari_close();
			exit(1);
		}

		cgiv(v);
		cgiv(g);
		#endif
	}

	free(numbers);
}

#if BITSIZE == 32
void sum_of_divisors_odd2(const ulong blocksize, const ulong l, uint * sigma, const uint * primes)
#else
void sum_of_divisors_odd2(const ulong blocksize, const ulong l, ulong * sigma, const uint * primes)
#endif
{
	if (((blocksize & 1) == 1) || ((l & 1) == 1))
	{
		perror("The parameters blocksize and l must be even.\n");
		exit(1);
	}

	const ulong max_prime = (ulong) sqrt(l + blocksize);

	ulong h, i, j, k, p, s, p_pow, step;

	ulong * numbers = (ulong *) malloc(sizeof(ulong) * (blocksize >> 1));

	for (i = 0, j = l + 1; i < (blocksize >> 1); i++, j += 2)
	{
		sigma[i] = 1;
		numbers[i] = j;
	}

	if (l == 0)
	{
		sigma[0] = 1;
	}

	ulong offset[100];

	for (i = 2, p = primes[i]; p <= max_prime; p = primes[++i])
	{
		k = 0;
		offset[k] = (p - (l % p)) % p;

		if ((offset[k] & 1) == 0)
		{
			offset[k] += p;
		}

		for (p_pow = p, step = p * p; p_pow <= (l + blocksize); p_pow *= p, step *= p)
		{
			if (offset[k] > blocksize)
			{
				break;
			}

			#ifdef DEBUG
			printf("p_pow=%lu, offset=%lu, step=%lu\t", p_pow, offset[k], step);
			#endif

			offset[++k] = (step - (l % step)) % step;

			if ((offset[k] & 1) == 0)
			{
				offset[k] += step;
			}

			for (s = 0, h = offset[k-1]; s < p; s++, h += (p_pow << 1))
			{
				if (h != offset[k])
				{
					for (j = h; j < blocksize; j += (step << 1))
					{
						#ifdef DEBUG
						printf(", %lu (%lu)", l + j, (j-1)>>1);
						#endif

						numbers[(j-1) >> 1] /= p_pow;
						sigma[(j-1) >> 1] *= ((step*p_pow - 1) / (p - 1));
					}
				}
			}

			#ifdef DEBUG
			printf("\n");
			#endif
		}
	}

	ulong j_squared;

	for (i = 0, j = l + 1, j_squared = j * j; i < (blocksize >> 1); i++, j_squared += ((j + 1) << 2), j += 2)
	{
		if (numbers[i] > 1)
		{
			sigma[i] *= (numbers[i] * (numbers[i] + 1L) + 1L);
		}

		s = sigma[i] - j_squared;

		#ifdef DEBUG
		printf("sigma(%lu^2=%lu)=%lu, s(%lu)=%lu\n", j, j_squared, (ulong) sigma[i], j*j, s);
		#endif

		#ifdef WITH_PARI
		GEN g = stoi(j_squared);
		GEN v = gsumdivk(g, 1);

		if (sigma[i] != gtolong(v))
		{
			char err[240];
			sprintf(err, "Error, sigma(%lu^2=%lu)=%lu, but we computed %lu.\n", j, j_squared, gtolong(v), (ulong) sigma[i]);
			perror(err);
			cgiv(v);
			cgiv(g);
			pari_close();
			exit(1);
		}

		cgiv(v);
		cgiv(g);
		#endif
	}

	free(numbers);
}
