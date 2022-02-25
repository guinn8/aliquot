/**
 * @file math_macros.h
 * @author Gavin Guinn (gavinguinn1@gmail.com)
 * @brief 
 * @date 2022-02-24
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef POMYANG_KPARENT_INC_MATH_MACROS_H_
#define POMYANG_KPARENT_INC_MATH_MACROS_H_

/** @brief Convert bytes to gigabytes. */
#define BYTES_TO_GB 0.000000001

/** @brief Constant for 2/3. */
#define TWO_THIRDS .6666666666666666666666

/** @brief 1 if a number is even. */
#define EVEN(x) (0 == (x) % 2)

/** @brief Square a number. */
#define SQUARE(x) ((x) * (x))

/** @brief Return the smaller number. */
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/** @brief Check if x divides y evenly. */
#define DIVIDES(x, y) (0 == y % x)

/** @brief Offset an even number into an array. */
#define F_OFFSET(x) ((x / 2) - 1)

/** @brief Get the even number represented by this array index. */
#define F_DE_OFFSET(x) ((x + 1) * 2)

/** @brief quietable logging macro. */
#define LOG(shush, fmt, ...)                          \
    do {                                              \
        if (!shush) {                                 \
            fprintf(stderr, "[%s:%d] " fmt, __FILE__, \
                    __LINE__, __VA_ARGS__);           \
        }                                             \
    } while (0)


#ifdef DEBUG_ASSERT_ON
/** @brief Turn off asserts for production runs.*/
#define DEBUG_ASSERT(x) x
#else
/** @brief Turn off asserts for production runs.*/
#define DEBUG_ASSERT(x)
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"  // uint64_t used in enum

/**
 * @brief Counts of nonaliquots typed out from (Chum et al.)
 * 
 */
typedef enum {
    _10_TO_THE_4 = 1212,
    _10_TO_THE_5 = 13863,
    _10_TO_THE_6 = 150232,
    _10_TO_THE_7 = 1574973,
    _10_TO_THE_8 = 16246940,
    _10_TO_THE_9 = 165826606,
    _10_TO_THE_10 = 1681871718,
    _10_TO_THE_11 = 16988116409,
    _10_TO_THE_12 = 171128671374,
    _2_TO_THE_40 = 188206399403,
} nonaliquot_counts_t;
#pragma GCC diagnostic pop

#endif  // POMYANG_KPARENT_INC_MATH_MACROS_H_
