/** @copyright see source */
#ifndef POMYANG_KPARENT_INC_PACKEDARRAY_H_
#define POMYANG_KPARENT_INC_PACKEDARRAY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <omp.h>
#include <stdint.h>

/*

PackedArray principle:
  . compact storage of <= 32 bits items
  . items are tightly packed into a buffer of uint32_t integers

PackedArray requirements:
  . you must know in advance how many bits are needed to hold a single item
  . you must know in advance how many items you want to store
  . when packing, behavior is undefined if items have more than bitsPerItem bits

PackedArray general in memory representation:
  |-------------------------------------------------- - - -
  |       b0       |       b1       |       b2       |
  |-------------------------------------------------- - - -
  | i0 | i1 | i2 | i3 | i4 | i5 | i6 | i7 | i8 | i9 |
  |-------------------------------------------------- - - -

  . items are tightly packed together
  . several items end up inside the same buffer cell, e.g. i0, i1, i2
  . some items span two buffer cells, e.g. i3, i6

*/

/**
 * @brief Protects buffer underlying PackedArray with an array of mutexs.
 * @struct lock_info_t
 *
 * @var lock_info_t::num_locks
 * Number of mutexes protecting buffer undlying the PackedArray from concorent access.
 *
 * @var lock_info_t::lock_interval
 * Computed value used internally.
 *
 * @var lock_info_t::locks
 * Calloc'ed buffer of locks.
 */
typedef struct {
    uint32_t num_locks;
    uint32_t lock_interval;
    omp_lock_t* locks;
} lock_info_t;

/**
 * @brief PackedArray handle.
 * @struct PackedArray
 *
 * @var PackedArray::bitsPerItem
 * Bits to be used per number stored in the PackedArray.
 *
 * @var PackedArray::count
 * Length of PackedArray.
 *
 * @var PackedArray::lock_info
 * Mutex array handle.
 *
 * @var PackedArray::padding
 * ?
 *
 * @var PackedArray::buffer
 * Buffer underlying PackedArray.
 */
typedef struct {
    uint32_t bitsPerItem;
    uint64_t count;
    lock_info_t lock_info;
    uint32_t padding[2];
    uint32_t buffer[];
} PackedArray;

// creation / destruction
PackedArray* PackedArray_create(uint32_t bitsPerItem, uint64_t count, uint32_t num_locks);
void PackedArray_destroy(PackedArray* a);
uint64_t PackedArray_estimate_heap(uint32_t bitsPerItem, uint64_t count, uint32_t num_locks);

// packing / unpacking
// offset is expressed in number of elements
void PackedArray_pack(PackedArray* a, const uint64_t offset, const uint32_t* in, uint64_t count);
void PackedArray_unpack(const PackedArray* a, const uint64_t offset, uint32_t* out, uint64_t count);

// single item access
void PackedArray_set(PackedArray* a, const uint64_t offset, const uint32_t in);
uint32_t PackedArray_get(const PackedArray* a, const uint64_t offset);

// helpers
uint32_t PackedArray_bufferSize(const PackedArray* a);
uint32_t PackedArray_computeBitsPerItem(const uint32_t* in, uint64_t count);
void PackedArray_unlock_offset(PackedArray* a, const uint64_t offset);
void PackedArray_lock_offset(PackedArray* a, const uint64_t offset);
#ifdef __cplusplus
}
#endif

#endif  // POMYANG_KPARENT_INC_PACKEDARRAY_H_
