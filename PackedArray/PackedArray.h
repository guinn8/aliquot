#ifndef PACKEDARRAY_H
#define PACKEDARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/*

PackedArray principle:
  . compact storage of <= 32 bits items
  . items are tightly packed into a buffer of uint64_t integers

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

struct _PackedArray
{
  uint64_t bitsPerItem;
  uint64_t count;

  uint64_t padding[2];
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4200)
#endif // #ifdef _MSC_VER
  uint64_t buffer[];
#ifdef _MSC_VER
#pragma warning(pop)
#endif // #ifdef _MSC_VER
};
typedef struct _PackedArray PackedArray;

// creation / destruction
PackedArray* PackedArray_create(uint64_t bitsPerItem, uint64_t count);
void PackedArray_destroy(PackedArray* a);

// packing / unpacking
// offset is expressed in number of elements
void PackedArray_pack(PackedArray* a, const uint64_t offset, const uint64_t* in, uint64_t count);
void PackedArray_unpack(const PackedArray* a, const uint64_t offset, uint64_t* out, uint64_t count);

// single item access
void PackedArray_set(PackedArray* a, const uint64_t offset, const uint64_t in);
uint64_t PackedArray_get(const PackedArray* a, const uint64_t offset);

// helpers
uint64_t PackedArray_bufferSize(const PackedArray* a);
uint64_t PackedArray_computeBitsPerItem(const uint64_t* in, uint64_t count);

#ifdef __cplusplus
}
#endif

#endif // #ifndef PACKEDARRAY_H
