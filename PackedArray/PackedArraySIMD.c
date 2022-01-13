// see README.md for usage instructions.
// (‑●‑●)> released under the WTFPL v2 license, by Gregory Pakosz (@gpakosz)

#ifndef PACKEDARRAY_SELF
#define PACKEDARRAY_SELF "PackedArraySIMD.c"
#endif

#ifdef PACKEDARRAY_IMPL

#ifndef PACKEDARRAY_JOIN
#define PACKEDARRAY_JOIN(lhs, rhs)    PACKEDARRAY_JOIN_(lhs, rhs)
#define PACKEDARRAY_JOIN_(lhs, rhs)   PACKEDARRAY_JOIN__(lhs, rhs)
#define PACKEDARRAY_JOIN__(lhs, rhs)  lhs##rhs
#endif // #ifndef PACKEDARRAY_JOIN

#ifndef PACKEDARRAY_IMPL_BITS_PER_ITEM
#error PACKEDARRAY_IMPL_BITS_PER_ITEM undefined
#endif // #ifndef PACKEDARRAY_IMPL_BITS_PER_ITEM

#if defined(PACKEDARRAY_IMPL_PACK_CASES) || defined(PACKEDARRAY_IMPL_UNPACK_CASES)

#ifndef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 0
#elif PACKEDARRAY_IMPL_CASE_I == 0
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 1
#elif PACKEDARRAY_IMPL_CASE_I == 1
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 2
#elif PACKEDARRAY_IMPL_CASE_I == 2
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 3
#elif PACKEDARRAY_IMPL_CASE_I == 3
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 4
#elif PACKEDARRAY_IMPL_CASE_I == 4
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 5
#elif PACKEDARRAY_IMPL_CASE_I == 5
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 6
#elif PACKEDARRAY_IMPL_CASE_I == 6
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 7
#elif PACKEDARRAY_IMPL_CASE_I == 7
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 8
#elif PACKEDARRAY_IMPL_CASE_I == 8
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 9
#elif PACKEDARRAY_IMPL_CASE_I == 9
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 10
#elif PACKEDARRAY_IMPL_CASE_I == 10
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 11
#elif PACKEDARRAY_IMPL_CASE_I == 11
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 12
#elif PACKEDARRAY_IMPL_CASE_I == 12
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 13
#elif PACKEDARRAY_IMPL_CASE_I == 13
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 14
#elif PACKEDARRAY_IMPL_CASE_I == 14
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 15
#elif PACKEDARRAY_IMPL_CASE_I == 15
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 16
#elif PACKEDARRAY_IMPL_CASE_I == 16
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 17
#elif PACKEDARRAY_IMPL_CASE_I == 17
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 18
#elif PACKEDARRAY_IMPL_CASE_I == 18
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 19
#elif PACKEDARRAY_IMPL_CASE_I == 19
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 20
#elif PACKEDARRAY_IMPL_CASE_I == 20
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 21
#elif PACKEDARRAY_IMPL_CASE_I == 21
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 22
#elif PACKEDARRAY_IMPL_CASE_I == 22
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 23
#elif PACKEDARRAY_IMPL_CASE_I == 23
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 24
#elif PACKEDARRAY_IMPL_CASE_I == 24
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 25
#elif PACKEDARRAY_IMPL_CASE_I == 25
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 26
#elif PACKEDARRAY_IMPL_CASE_I == 26
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 27
#elif PACKEDARRAY_IMPL_CASE_I == 27
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 28
#elif PACKEDARRAY_IMPL_CASE_I == 28
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 29
#elif PACKEDARRAY_IMPL_CASE_I == 29
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 30
#elif PACKEDARRAY_IMPL_CASE_I == 30
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 31
#elif PACKEDARRAY_IMPL_CASE_I == 31
#undef PACKEDARRAY_IMPL_CASE_I
#define PACKEDARRAY_IMPL_CASE_I 32
#endif // #ifndef PACKEDARRAY_IMPL_CASE_I

#ifndef PACKEDARRAY_IMPL_BITS_AVAILABLE
#define PACKEDARRAY_IMPL_BITS_AVAILABLE (32 - ((PACKEDARRAY_IMPL_CASE_I * PACKEDARRAY_IMPL_BITS_PER_ITEM) % 32))
#endif
#ifndef PACKEDARRAY_IMPL_START_BIT
#define PACKEDARRAY_IMPL_START_BIT ((PACKEDARRAY_IMPL_CASE_I * PACKEDARRAY_IMPL_BITS_PER_ITEM) % 32)
#endif

#if defined(PACKEDARRAY_IMPL_PACK_CASES)

#ifndef PACKEDARRAY_IMPL_PACK_CASE_BREAK
#define PACKEDARRAY_IMPL_PACK_CASE_BREAK
#endif

      case PACKEDARRAY_IMPL_CASE_I:
#if (PACKEDARRAY_IMPL_BITS_PER_ITEM <= PACKEDARRAY_IMPL_BITS_AVAILABLE)
        in_4 = PackedArray_loadu_uint32x4(in);
        packed = PackedArray_vsli0_uint32x4(packed, in_4, PACKEDARRAY_IMPL_START_BIT);
        in += 4;
#if (PACKEDARRAY_IMPL_BITS_PER_ITEM == PACKEDARRAY_IMPL_BITS_AVAILABLE)
        PackedArray_store_uint32x4(out, packed);
        out += 4;
        packed = PackedArray_uint32x4_zero;
#endif
#else
        in_4 = PackedArray_loadu_uint32x4(in);
        packed = PackedArray_vsli0_uint32x4(packed, in_4, PACKEDARRAY_IMPL_START_BIT);
        PackedArray_store_uint32x4(out, packed);
        out += 4;
        packed = PackedArray_shr_uint32x4(in_4, PACKEDARRAY_IMPL_BITS_AVAILABLE);
        in += 4;
#endif
        PACKEDARRAY_IMPL_PACK_CASE_BREAK

#if PACKEDARRAY_IMPL_CASE_I < 31
#include PACKEDARRAY_SELF
#else
#undef PACKEDARRAY_IMPL_CASE_I
#undef PACKEDARRAY_IMPL_PACK_CASE_BREAK
#undef PACKEDARRAY_IMPL_PACK_CASES
#endif

#elif defined(PACKEDARRAY_IMPL_UNPACK_CASES) // #if defined(PACKEDARRAY_IMPL_PACK_CASES)

#ifndef PACKEDARRAY_IMPL_UNPACK_CASE_BREAK
#define PACKEDARRAY_IMPL_UNPACK_CASE_BREAK
#endif

      case PACKEDARRAY_IMPL_CASE_I:
#if (PACKEDARRAY_IMPL_BITS_PER_ITEM <= PACKEDARRAY_IMPL_BITS_AVAILABLE)
        out_4 = PackedArray_and_uint32x4(PackedArray_shr_uint32x4(packed, PACKEDARRAY_IMPL_START_BIT), PackedArray_set_uint32x4(PACKEDARRAY_IMPL_MASK));
        PackedArray_storeu_uint32x4(out, out_4);
        out += 4;
        PACKEDARRAY_IMPL_UNPACK_CASE_BREAK
#if (PACKEDARRAY_IMPL_CASE_I < 31) && (PACKEDARRAY_IMPL_BITS_PER_ITEM == PACKEDARRAY_IMPL_BITS_AVAILABLE)
        in += 4;
        packed = PackedArray_load_uint32x4(in);
#endif
#else
        out_4 = PackedArray_shr_uint32x4(packed, PACKEDARRAY_IMPL_START_BIT);
        in += 4;
        packed = PackedArray_load_uint32x4(in);
        out_4 = PackedArray_vsli0_uint32x4(out_4, packed, PACKEDARRAY_IMPL_BITS_AVAILABLE);
        out_4 = PackedArray_and_uint32x4(out_4, PackedArray_set_uint32x4(PACKEDARRAY_IMPL_MASK));
        PackedArray_storeu_uint32x4(out, out_4);
        out += 4;
        PACKEDARRAY_IMPL_UNPACK_CASE_BREAK
#endif

#if PACKEDARRAY_IMPL_CASE_I < 31
#include PACKEDARRAY_SELF
#else
#undef PACKEDARRAY_IMPL_CASE_I
#undef PACKEDARRAY_IMPL_UNPACK_CASE_BREAK
#undef PACKEDARRAY_IMPL_UNPACK_CASES
#endif

#endif // #elif defined(PACKEDARRAY_IMPL_UNPACK_CASES)

#else // #if defined(PACKEDARRAY_IMPL_PACK_CASES) || defined(PACKEDARRAY_IMPL_UNPACK_CASES)

#ifndef PACKEDARRAY_IMPL_MASK
#define PACKEDARRAY_IMPL_MASK (uint32_t)((1ULL << PACKEDARRAY_IMPL_BITS_PER_ITEM) - 1)
#endif

void PACKEDARRAY_JOIN(__PackedArray_pack_, PACKEDARRAY_IMPL_BITS_PER_ITEM)(uint32_t* __restrict buffer, uint32_t offset, const uint32_t* __restrict in, uint32_t count)
{
  uint32_t pre, post;
  uint32_t* __restrict out;
  const uint32_t* __restrict end;
  uint32_t startBit;
  PackedArray_uint32x4_t packed, in_4, mask;
  uint32_t offset_4;

  pre = (offset + 3) / 4 * 4 - offset;
  pre = pre > count ? count : pre;

  if (pre > 0)
  {
    __PackedArray_pack_scalar(buffer, PACKEDARRAY_IMPL_BITS_PER_ITEM, PACKEDARRAY_IMPL_MASK, offset, in, pre);
    offset += pre;
    in += pre;
    count -= pre;
  }

  post = count % 4;
  count -= post;

  if (count > 0)
  {
    out = &buffer[(offset / 4 * PACKEDARRAY_IMPL_BITS_PER_ITEM) / 32 * 4];
    startBit = (offset / 4 * PACKEDARRAY_IMPL_BITS_PER_ITEM) % 32;
    packed = PackedArray_load_uint32x4(out);
    mask = PackedArray_sub_uint32x4(PackedArray_shl_uint32x4(PackedArray_set_uint32x4(1), startBit), PackedArray_set_uint32x4(1));
    packed = PackedArray_and_uint32x4(packed, mask);

    offset_4 = offset % 128;
    offset += count;

    if (count >= 128 - offset_4)
    {
      int32_t n;

      n = (count + offset_4) / 128;
      count -= 128 * n - offset_4;
      switch (offset_4 / 4)
      {
        do
        {
#define PACKEDARRAY_IMPL_PACK_CASES
#include PACKEDARRAY_SELF
        } while (--n > 0);
      }

      if (count == 0)
        goto PACKEDARRAY_JOIN(PACKEDARRAY_JOIN(__PackedArray_pack_, PACKEDARRAY_IMPL_BITS_PER_ITEM), _post);

      offset_4 = 0;
      startBit = 0;
    }

    end = in + count;
    switch (offset_4 / 4)
    {
#define PACKEDARRAY_IMPL_PACK_CASES
#define PACKEDARRAY_IMPL_PACK_CASE_BREAK \
      if (in == end)\
        break;
#include PACKEDARRAY_SELF
    }
    PACKEDARRAY_ASSERT(in == end);
    if ((count / 4 * PACKEDARRAY_IMPL_BITS_PER_ITEM + startBit) % 32)
    {
      in_4 = PackedArray_loadu_uint32x4(out);
      mask = PackedArray_sub_uint32x4(PackedArray_shl_uint32x4(PackedArray_set_uint32x4(1), ((count / 4 * PACKEDARRAY_IMPL_BITS_PER_ITEM + startBit - 1) % 32) + 1), PackedArray_set_uint32x4(1));
      in_4 = PackedArray_andnot_uint32x4(in_4, mask);
      packed = PackedArray_or_uint32x4(packed, in_4);
      PackedArray_store_uint32x4(out, packed);
    }
  }

  PACKEDARRAY_JOIN(PACKEDARRAY_JOIN(__PackedArray_pack_, PACKEDARRAY_IMPL_BITS_PER_ITEM), _post):
  if (post > 0)
    __PackedArray_pack_scalar(buffer, PACKEDARRAY_IMPL_BITS_PER_ITEM, PACKEDARRAY_IMPL_MASK, offset, in, post);
}

void PACKEDARRAY_JOIN(__PackedArray_unpack_, PACKEDARRAY_IMPL_BITS_PER_ITEM)(const uint32_t* __restrict buffer, uint32_t offset, uint32_t* __restrict out, uint32_t count)
{
  uint32_t pre, post;
  const uint32_t* __restrict in;
  const uint32_t* __restrict end;
  PackedArray_uint32x4_t packed, out_4;
  uint32_t offset_4;

  pre = (offset + 3) / 4 * 4 - offset;
  pre = pre > count ? count : pre;

  if (pre > 0)
  {
    __PackedArray_unpack_scalar(buffer, PACKEDARRAY_IMPL_BITS_PER_ITEM, PACKEDARRAY_IMPL_MASK, offset, out, pre);
    offset += pre;
    out += pre;
    count -= pre;
  }

  post = count % 4;
  count -= post;

  if (count > 0)
  {
    in = &buffer[(offset / 4 * PACKEDARRAY_IMPL_BITS_PER_ITEM) / 32 * 4];
    packed = PackedArray_load_uint32x4(in);

    offset_4 = offset % 128;
    offset += count;

    if (count >= 128 - offset_4)
    {
      int32_t n;

      n = (count + offset_4) / 128;
      count -= 128 * n - offset_4;
      switch (offset_4 / 4)
      {
        do
        {
          in += 4;
          packed = PackedArray_load_uint32x4(in);
#define PACKEDARRAY_IMPL_UNPACK_CASES
#include PACKEDARRAY_SELF
        } while (--n > 0);
      }

      if (count == 0)
        goto PACKEDARRAY_JOIN(PACKEDARRAY_JOIN(__PackedArray_unpack_, PACKEDARRAY_IMPL_BITS_PER_ITEM), _post);

      in += 4;
      packed = PackedArray_loadu_uint32x4(in);
      offset_4 = 0;
    }

    end = out + count;
    switch (offset_4 / 4)
    {
#define PACKEDARRAY_IMPL_UNPACK_CASES
#define PACKEDARRAY_IMPL_UNPACK_CASE_BREAK \
      if (out == end)\
        break;
#include PACKEDARRAY_SELF
    }
    PACKEDARRAY_ASSERT(out == end);
  }

  PACKEDARRAY_JOIN(PACKEDARRAY_JOIN(__PackedArray_unpack_, PACKEDARRAY_IMPL_BITS_PER_ITEM), _post):
  if (post > 0)
    __PackedArray_unpack_scalar(buffer, PACKEDARRAY_IMPL_BITS_PER_ITEM, PACKEDARRAY_IMPL_MASK, offset, out, post);
}

#undef PACKEDARRAY_IMPL_BITS_PER_ITEM
#undef PACKEDARRAY_IMPL_BITS_AVAILABLE
#undef PACKEDARRAY_IMPL_START_BIT
#undef PACKEDARRAY_IMPL_START_MASK

#endif // #if defined(PACKEDARRAY_IMPL_PACK_CASES) || defined(PACKEDARRAY_IMPL_UNPACK_CASES)

#else

#include "PackedArray.h"

#if !defined(PACKEDARRAY_ASSERT)
#include <assert.h>
#define PACKEDARRAY_ASSERT(expression) assert(expression)
#endif

#include <stddef.h>

static void __PackedArray_pack_scalar(uint32_t* buffer, const uint32_t bitsPerItem, const uint32_t mask, uint32_t offset, const uint32_t* in, uint32_t count)
{
  uint32_t* __restrict out;
  uint32_t startBit;
  uint32_t bitsAvailable;

  PACKEDARRAY_ASSERT(buffer != NULL);
  PACKEDARRAY_ASSERT(in != NULL);
  PACKEDARRAY_ASSERT(count != 0);

  while (count--)
  {
    uint32_t value = *in++;

    PACKEDARRAY_ASSERT(0 == (~mask & value));

    out = &buffer[((uint64_t)offset / 4 * (uint64_t)bitsPerItem) / 32 * 4 + offset % 4];
    startBit = ((uint64_t)offset / 4 * (uint64_t)bitsPerItem) % 32;
    ++offset;

    bitsAvailable = 32 - startBit;

    if (bitsPerItem <= bitsAvailable)
    {
      out[0] = (out[0] & ~(mask << startBit)) | (value << startBit);
    }
    else
    {
      // value spans 2 buffer cells
      uint32_t low, high;

      low = value << startBit;
      high = value >> bitsAvailable;

      out[0] = (out[0] & ~(mask << startBit)) | low;

      out[4] = (out[4] & ~(mask >> bitsAvailable)) | high;
    }
  }
}

static void __PackedArray_unpack_scalar(const uint32_t* buffer, const uint32_t bitsPerItem, const uint32_t mask, uint32_t offset, uint32_t* out, uint32_t count)
{
  const uint32_t* __restrict in;
  uint32_t startBit;
  uint32_t bitsAvailable;

  PACKEDARRAY_ASSERT(buffer != NULL);
  PACKEDARRAY_ASSERT(out != NULL);
  PACKEDARRAY_ASSERT(count != 0);

  while (count--)
  {
    uint32_t value;

    in = &buffer[((uint64_t)offset / 4 * (uint64_t)bitsPerItem) / 32 * 4 + offset % 4];
    startBit = ((uint64_t)offset / 4 * (uint64_t)bitsPerItem) % 32;
    ++offset;

    bitsAvailable = 32 - startBit;

    if (bitsPerItem <= bitsAvailable)
    {
      value = (in[0] >> startBit) & mask;
    }
    else
    {
      // value spans 2 buffer cells
      uint32_t low, high;

      low = in[0] >> startBit;
      high = in[4] << bitsAvailable;

      value = (low | high) & mask;
    }
    *out++ = value;
  }
}

#if defined(__SSE2__) || defined(_M_IX86) || defined(_M_X64)

#include <emmintrin.h>

#define PackedArray_uint32x4_t                      __m128i
#define PackedArray_uint32x4_zero                   _mm_setzero_si128()
#define PackedArray_set_uint32x4(i)                 _mm_set1_epi32(i)
#define PackedArray_sub_uint32x4(lhs, rhs)          _mm_sub_epi32(lhs, rhs)
#define PackedArray_loadu_uint32x4(ptr)             _mm_loadu_si128((const __m128i*)ptr)
#define PackedArray_storeu_uint32x4(ptr, v)         _mm_storeu_si128((__m128i*)ptr, v)
#define PackedArray_load_uint32x4(ptr)              _mm_load_si128((const __m128i*)ptr)
#define PackedArray_store_uint32x4(ptr, v)          _mm_store_si128((__m128i*)ptr, v)
#define PackedArray_shl_uint32x4(v, shift)          _mm_slli_epi32(v, shift)
#define PackedArray_shr_uint32x4(v, shift)          _mm_srli_epi32(v, shift)
#define PackedArray_or_uint32x4(lhs, rhs)           _mm_or_si128(lhs, rhs)
#define PackedArray_and_uint32x4(lhs, rhs)          _mm_and_si128(lhs, rhs)
#define PackedArray_andnot_uint32x4(lhs, rhs)       _mm_andnot_si128(rhs, lhs)
// assumes dst bits are cleared at inserted bit positions
#define PackedArray_vsli0_uint32x4(dst, src, shift) PackedArray_or_uint32x4(dst, PackedArray_shl_uint32x4(src, shift))

#elif defined(__ARM_NEON__) || defined(_M_ARM)

#include <arm_neon.h>

#define PackedArray_uint32x4_t                      uint32x4_t
#define PackedArray_uint32x4_zero                   vdupq_n_u32(0)
#define PackedArray_set_uint32x4(i)                 vdupq_n_u32(i)
#define PackedArray_sub_uint32x4(lhs, rhs)          vsubq_u32(lhs, rhs)
#define PackedArray_loadu_uint32x4(ptr)             vld1q_u32((const uint32_t*)ptr)
#define PackedArray_storeu_uint32x4(ptr, v)         vst1q_u32(ptr, v)
#if defined(__GNUC__)
// because  __builtin_assume_aligned isn't always available...

typedef uint32_t __attribute__((aligned(16)))       PackedArray_aligned_uint32_t;
#define PackedArray_load_uint32x4(ptr)              vld1q_u32((const PackedArray_aligned_uint32_t*)ptr)
#define PackedArray_store_uint32x4(ptr, v)          vst1q_u32((PackedArray_aligned_uint32_t*)ptr, v)
#else
#define PackedArray_load_uint32x4(ptr)              vld1q_u32((const uint32_t*)ptr)
#define PackedArray_store_uint32x4(ptr, v)          vst1q_u32(ptr, v)
#endif
#define PackedArray_shl_uint32x4(v, shift)          vshlq_u32(v, vdupq_n_s32(shift))
#define PackedArray_shr_uint32x4(v, shift)          vshlq_u32(v, vdupq_n_s32(-shift))
#define PackedArray_or_uint32x4(lhs, rhs)           vorrq_u32(lhs, rhs)
#define PackedArray_and_uint32x4(lhs, rhs)          vandq_u32(lhs, rhs)
#define PackedArray_andnot_uint32x4(lhs, rhs)       vbicq_u32(lhs, rhs)
// assumes dst bits are cleared at inserted bit positions
#define PackedArray_vsli0_uint32x4(dst, src, shift) vsliq_n_u32(dst, src, shift)

#else

#error unsupported SIMD platform

#endif

#define PACKEDARRAY_IMPL
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 1
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 2
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 3
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 4
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 5
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 6
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 7
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 8
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 9
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 10
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 11
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 12
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 13
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 14
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 15
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 16
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 17
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 18
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 19
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 20
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 21
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 22
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 23
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 24
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 25
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 26
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 27
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 28
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 29
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 30
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 31
#include PACKEDARRAY_SELF
#define PACKEDARRAY_IMPL_BITS_PER_ITEM 32
#include PACKEDARRAY_SELF
#undef PACKEDARRAY_IMPL


#if !defined(PACKEDARRAY_ALIGNED_MALLOC) || !defined(PACKEDARRAY_FREE)
#include <stdlib.h>
#endif

#if !defined(PACKEDARRAY_ALIGNED_MALLOC)
#if defined (_MSC_VER)
#define PACKEDARRAY_ALIGNED_MALLOC(alignment, size) _aligned_malloc(size, alignment)
#elif defined (ANDROID) || defined (__ANDROID__)
#define PACKEDARRAY_ALIGNED_MALLOC(alignment, size) memalign(alignment, size)
#else
static void* __PackedArray_aligned_malloc(size_t alignment, size_t size)
{
  void* p = NULL;
  posix_memalign(&p, alignment, size);
  return p;
}
#define PACKEDARRAY_ALIGNED_MALLOC(alignment, size) __PackedArray_aligned_malloc(alignment, size)
#endif
#endif

#if !defined(PACKEDARRAY_FREE)
#if defined (_MSC_VER)
#define PACKEDARRAY_FREE(p) _aligned_free(p)
#else
#define PACKEDARRAY_FREE(p) free(p)
#endif
#endif

PackedArray* PackedArray_create(uint32_t bitsPerItem, uint32_t count)
{
  PackedArray* a;
  size_t bufferSize;

  PACKEDARRAY_ASSERT(bitsPerItem > 0);
  PACKEDARRAY_ASSERT(bitsPerItem <= 32);

  bufferSize = sizeof(uint32_t) * (((uint64_t)count / 4 * (uint64_t)bitsPerItem + 31) / 32 * 4);
  bufferSize += count < 4 ? sizeof(uint32_t) * count : sizeof(uint32_t) * 4;
  a = (PackedArray*)PACKEDARRAY_ALIGNED_MALLOC(16, sizeof(PackedArray) + bufferSize);
  PACKEDARRAY_ASSERT((uint64_t)a->buffer % 16 == 0);

  if (a != NULL)
  {
    a->buffer[bufferSize / sizeof(uint32_t) - 1] = 0;
    a->bitsPerItem = bitsPerItem;
    a->count = count;
  }

  return a;
}

void PackedArray_destroy(PackedArray* a)
{
  PACKEDARRAY_ASSERT(a);
  PACKEDARRAY_FREE(a);
}

void PackedArray_pack(PackedArray* a, const uint32_t offset, const uint32_t* in, uint32_t count)
{
  PACKEDARRAY_ASSERT(a != NULL);
  PACKEDARRAY_ASSERT(in != NULL);

  switch (a->bitsPerItem)
  {
    case 1:   __PackedArray_pack_1(a->buffer, offset, in, count); break;
    case 2:   __PackedArray_pack_2(a->buffer, offset, in, count); break;
    case 3:   __PackedArray_pack_3(a->buffer, offset, in, count); break;
    case 4:   __PackedArray_pack_4(a->buffer, offset, in, count); break;
    case 5:   __PackedArray_pack_5(a->buffer, offset, in, count); break;
    case 6:   __PackedArray_pack_6(a->buffer, offset, in, count); break;
    case 7:   __PackedArray_pack_7(a->buffer, offset, in, count); break;
    case 8:   __PackedArray_pack_8(a->buffer, offset, in, count); break;
    case 9:   __PackedArray_pack_9(a->buffer, offset, in, count); break;
    case 10:  __PackedArray_pack_10(a->buffer, offset, in, count); break;
    case 11:  __PackedArray_pack_11(a->buffer, offset, in, count); break;
    case 12:  __PackedArray_pack_12(a->buffer, offset, in, count); break;
    case 13:  __PackedArray_pack_13(a->buffer, offset, in, count); break;
    case 14:  __PackedArray_pack_14(a->buffer, offset, in, count); break;
    case 15:  __PackedArray_pack_15(a->buffer, offset, in, count); break;
    case 16:  __PackedArray_pack_16(a->buffer, offset, in, count); break;
    case 17:  __PackedArray_pack_17(a->buffer, offset, in, count); break;
    case 18:  __PackedArray_pack_18(a->buffer, offset, in, count); break;
    case 19:  __PackedArray_pack_19(a->buffer, offset, in, count); break;
    case 20:  __PackedArray_pack_20(a->buffer, offset, in, count); break;
    case 21:  __PackedArray_pack_21(a->buffer, offset, in, count); break;
    case 22:  __PackedArray_pack_22(a->buffer, offset, in, count); break;
    case 23:  __PackedArray_pack_23(a->buffer, offset, in, count); break;
    case 24:  __PackedArray_pack_24(a->buffer, offset, in, count); break;
    case 25:  __PackedArray_pack_25(a->buffer, offset, in, count); break;
    case 26:  __PackedArray_pack_26(a->buffer, offset, in, count); break;
    case 27:  __PackedArray_pack_27(a->buffer, offset, in, count); break;
    case 28:  __PackedArray_pack_28(a->buffer, offset, in, count); break;
    case 29:  __PackedArray_pack_29(a->buffer, offset, in, count); break;
    case 30:  __PackedArray_pack_30(a->buffer, offset, in, count); break;
    case 31:  __PackedArray_pack_31(a->buffer, offset, in, count); break;
    case 32:  __PackedArray_pack_32(a->buffer, offset, in, count); break;
  }
}

void PackedArray_unpack(const PackedArray* a, const uint32_t offset, uint32_t* out, uint32_t count)
{
  PACKEDARRAY_ASSERT(a != NULL);
  PACKEDARRAY_ASSERT(out != NULL);

  switch (a->bitsPerItem)
  {
    case 1:   __PackedArray_unpack_1(a->buffer, offset, out, count); break;
    case 2:   __PackedArray_unpack_2(a->buffer, offset, out, count); break;
    case 3:   __PackedArray_unpack_3(a->buffer, offset, out, count); break;
    case 4:   __PackedArray_unpack_4(a->buffer, offset, out, count); break;
    case 5:   __PackedArray_unpack_5(a->buffer, offset, out, count); break;
    case 6:   __PackedArray_unpack_6(a->buffer, offset, out, count); break;
    case 7:   __PackedArray_unpack_7(a->buffer, offset, out, count); break;
    case 8:   __PackedArray_unpack_8(a->buffer, offset, out, count); break;
    case 9:   __PackedArray_unpack_9(a->buffer, offset, out, count); break;
    case 10:  __PackedArray_unpack_10(a->buffer, offset, out, count); break;
    case 11:  __PackedArray_unpack_11(a->buffer, offset, out, count); break;
    case 12:  __PackedArray_unpack_12(a->buffer, offset, out, count); break;
    case 13:  __PackedArray_unpack_13(a->buffer, offset, out, count); break;
    case 14:  __PackedArray_unpack_14(a->buffer, offset, out, count); break;
    case 15:  __PackedArray_unpack_15(a->buffer, offset, out, count); break;
    case 16:  __PackedArray_unpack_16(a->buffer, offset, out, count); break;
    case 17:  __PackedArray_unpack_17(a->buffer, offset, out, count); break;
    case 18:  __PackedArray_unpack_18(a->buffer, offset, out, count); break;
    case 19:  __PackedArray_unpack_19(a->buffer, offset, out, count); break;
    case 20:  __PackedArray_unpack_20(a->buffer, offset, out, count); break;
    case 21:  __PackedArray_unpack_21(a->buffer, offset, out, count); break;
    case 22:  __PackedArray_unpack_22(a->buffer, offset, out, count); break;
    case 23:  __PackedArray_unpack_23(a->buffer, offset, out, count); break;
    case 24:  __PackedArray_unpack_24(a->buffer, offset, out, count); break;
    case 25:  __PackedArray_unpack_25(a->buffer, offset, out, count); break;
    case 26:  __PackedArray_unpack_26(a->buffer, offset, out, count); break;
    case 27:  __PackedArray_unpack_27(a->buffer, offset, out, count); break;
    case 28:  __PackedArray_unpack_28(a->buffer, offset, out, count); break;
    case 29:  __PackedArray_unpack_29(a->buffer, offset, out, count); break;
    case 30:  __PackedArray_unpack_30(a->buffer, offset, out, count); break;
    case 31:  __PackedArray_unpack_31(a->buffer, offset, out, count); break;
    case 32:  __PackedArray_unpack_32(a->buffer, offset, out, count); break;
  }
}

void PackedArray_set(PackedArray* a, const uint32_t offset, const uint32_t in)
{
  uint32_t* __restrict out;
  uint32_t bitsPerItem;
  uint32_t startBit;
  uint32_t bitsAvailable;
  uint32_t mask;

  PACKEDARRAY_ASSERT(a != NULL);

  bitsPerItem = a->bitsPerItem;

  out = &a->buffer[4 * (((uint64_t)offset / 4 * (uint64_t)bitsPerItem) / 32) + (offset % 4)];
  startBit = ((uint64_t)offset / 4 * (uint64_t)bitsPerItem) % 32;

  bitsAvailable = 32 - startBit;

  mask = (uint32_t)(1ULL << bitsPerItem) - 1;
  PACKEDARRAY_ASSERT(0 == (~mask & in));

  if (bitsPerItem <= bitsAvailable)
  {
    out[0] = (out[0] & ~(mask << startBit)) | (in << startBit);
  }
  else
  {
    // value spans 2 buffer cells
    uint32_t low, high;

    low = in << startBit;
    high = in >> bitsAvailable;

    out[0] = (out[0] & ~(mask << startBit)) | low;

    out[4] = (out[4] & ~(mask >> (32 - startBit))) | high;
  }
}

uint32_t PackedArray_get(const PackedArray* a, const uint32_t offset)
{
  const uint32_t* __restrict in;
  uint32_t bitsPerItem;
  uint32_t startBit;
  uint32_t bitsAvailable;
  uint32_t mask;
  uint32_t out;

  PACKEDARRAY_ASSERT(a != NULL);

  bitsPerItem = a->bitsPerItem;

  in = &a->buffer[4 * (((uint64_t)offset / 4 * (uint64_t)bitsPerItem) / 32) + (offset % 4)];
  startBit = ((uint64_t)offset / 4 * (uint64_t)bitsPerItem) % 32;

  bitsAvailable = 32 - startBit;

  mask = (uint32_t)(1ULL << bitsPerItem) - 1;

  if (bitsPerItem <= bitsAvailable)
  {
    out = (in[0] >> startBit) & mask;
  }
  else
  {
    // value spans 2 buffer cells
    uint32_t low, high;

    low = in[0] >> startBit;
    high = in[4] << (32 - startBit);

    out = low ^ ((low ^ high) & (mask >> bitsAvailable << bitsAvailable));
  }

  return out;
}

uint32_t PackedArray_bufferSize(const PackedArray* a)
{
  size_t bufferSize;
  uint32_t bitsPerItem, count;

  PACKEDARRAY_ASSERT(a != NULL);

  bitsPerItem = a->bitsPerItem;
  count = a->count;

  bufferSize = ((uint64_t)count / 4 * (uint64_t)bitsPerItem + 31) / 32 * 4;
  bufferSize += count < 4 ? count : 4;

  return (uint32_t)bufferSize;
}

#if !(defined(_MSC_VER) && _MSC_VER >= 1400) && !defined(__GNUC__)
// log base 2 of an integer, aka the position of the highest bit set
static uint32_t __PackedArray_log2(uint32_t v)
{
  // references
  // http://aggregate.org/MAGIC
  // http://graphics.stanford.edu/~seander/bithacks.html

  static const uint32_t multiplyDeBruijnBitPosition[32] =
  {
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
  };

  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;

  return multiplyDeBruijnBitPosition[(uint32_t)(v * 0x7C4ACDDU) >> 27];
}
#endif

// position of the highest bit set
static int __PackedArray_highestBitSet(uint32_t v)
{
#if defined(_MSC_VER) && _MSC_VER >= 1400
  unsigned long index;
  return _BitScanReverse(&index, v) ? index : -1;
#elif defined(__GNUC__)
  return v == 0 ? -1 : 31 - __builtin_clz(v);
#else
  return v != 0 ? __PackedArray_log2(v) : -1;
#endif
}

uint32_t PackedArray_computeBitsPerItem(const uint32_t* in, uint32_t count)
{
  uint32_t i, in_max, bitsPerItem;

  in_max = 0;
  for (i = 0; i < count; ++i)
    in_max = in[i] > in_max ? in[i] : in_max;

  bitsPerItem = __PackedArray_highestBitSet(in_max) + 1;
  return bitsPerItem == 0 ? 1 : bitsPerItem;
}


// - 8< ------------------------------------------------------------------------

#if defined(PACKEDARRAY_SELF_TEST) && defined(PACKEDARRAY_SELF_BENCH)
  #error choose either PACKEDARRAY_SELF_TEST or PACKEDARRAY_SELF_BENCH
#endif

#if defined(PACKEDARRAY_SELF_TEST)

#undef NDEBUG // we want asserts
#include <assert.h>

#include <stdio.h>
#include <string.h> // memcmp

static void PackedArray_pack_reference(PackedArray* a, uint32_t offset, const uint32_t* in, uint32_t count)
{
  uint32_t* __restrict out;
  uint32_t bitsPerItem;
  uint32_t startBit;
  uint32_t bitsAvailable;
  uint32_t mask;

  assert(a != NULL);
  assert(in != NULL);

  bitsPerItem = a->bitsPerItem;

  mask = (uint32_t)((1ULL << bitsPerItem) - 1);

  while (count--)
  {
    uint32_t value = *in++;

    assert(0 == (~mask & value));

    out = &a->buffer[((uint64_t)offset / 4 * (uint64_t)bitsPerItem) / 32 * 4 + offset % 4];
    startBit = ((uint64_t)offset / 4 * (uint64_t)bitsPerItem) % 32;
    ++offset;

    bitsAvailable = 32 - startBit;

    if (bitsPerItem <= bitsAvailable)
    {
      out[0] = (out[0] & ~(mask << startBit)) | (value << startBit);
    }
    else
    {
      // value spans 2 buffer cells
      uint32_t low, high;

      low = value << startBit;
      high = value >> bitsAvailable;

      out[0] = (out[0] & ~(mask << startBit)) | low;

      out[4] = (out[4] & ~(mask >> bitsAvailable)) | high;
    }
  }
}

static void PackedArray_unpack_reference(const PackedArray* a, uint32_t offset, uint32_t* out, uint32_t count)
{
  const uint32_t* __restrict in;
  uint32_t bitsPerItem;
  uint32_t startBit;
  uint32_t bitsAvailable;
  uint32_t mask;

  assert(a != NULL);
  assert(out != NULL);

  bitsPerItem = a->bitsPerItem;

  mask = (uint32_t)((1ULL << bitsPerItem) - 1);

  while (count--)
  {
    uint32_t value;

    in = &a->buffer[((uint64_t)offset / 4 * (uint64_t)bitsPerItem) / 32 * 4 + offset % 4];
    startBit = ((uint64_t)offset / 4 * (uint64_t)bitsPerItem) % 32;
    ++offset;

    bitsAvailable = 32 - startBit;

    if (bitsPerItem <= bitsAvailable)
    {
      value = (in[0] >> startBit) & mask;
    }
    else
    {
      // value spans 2 buffer cells
      uint32_t low, high;

      low = in[0] >> startBit;
      high = in[4] << bitsAvailable;

      value = low ^ ((low ^ high) & (mask >> bitsAvailable << bitsAvailable));
    }
    *out++ = value;
  }
}

int main(void)
{
  uint32_t bitsPerItem;

  printf("-- PackedArray self test -------------------------------------------------------\n");
  printf("\n");

  printf("sizeof(PackedArray) = %d\n", (int)sizeof(PackedArray));
  printf("\n");

  printf("1 by 1 packing / unpacking:\n");
  for (bitsPerItem = 1; bitsPerItem <= 32; ++bitsPerItem)
  {
    uint32_t mask = (uint32_t)(1ULL << bitsPerItem) - 1;
    int count;

    for (count = 1; count <= 1024; ++count)
    {
      PackedArray* a1 = PackedArray_create(bitsPerItem, count);
      PackedArray* a2 = PackedArray_create(bitsPerItem, count);
      int i;

      assert(a1->count == (uint32_t)count);
      assert(a2->count == (uint32_t)count);
      assert(a1->bitsPerItem == bitsPerItem);
      assert(a2->bitsPerItem == bitsPerItem);
      assert(PackedArray_bufferSize(a1) == PackedArray_bufferSize(a2));

      for (i = 0; i < (int)PackedArray_bufferSize(a1); ++i)
        a1->buffer[i] = a2->buffer[i] = rand();

      for (i = 0; i < count; ++i)
      {
        uint32_t v1, v2;

        v1 = rand() & mask;
        v2 = v1 + 1;
        PackedArray_pack(a1, i, &v1, 1);
        PackedArray_pack_reference(a2, i, &v1, 1);
        assert(memcmp(a1->buffer, a2->buffer, sizeof(a1->buffer[0]) * PackedArray_bufferSize(a1)) == 0);

        PackedArray_set(a1, i, v1);
        assert(memcmp(a1->buffer, a2->buffer, sizeof(a1->buffer[0]) * PackedArray_bufferSize(a1)) == 0);

        PackedArray_unpack(a1, i, &v2, 1);
        assert(v1 == v2);
        PackedArray_unpack_reference(a2, i, &v2, 1);
        assert(v1 == v2);
        v2 = PackedArray_get(a2, i);
        assert(v1 == v2);
      }

      for (i = count - 1; i >= 0; --i)
      {
        uint32_t v1, v2;

        v1 = rand() & mask;
        v2 = v1 + 1;
        PackedArray_pack(a1, i, &v1, 1);
        PackedArray_pack_reference(a2, i, &v1, 1);
        assert(memcmp(a1->buffer, a2->buffer, sizeof(a1->buffer[0]) * PackedArray_bufferSize(a1)) == 0);

        PackedArray_set(a1, i, v1);
        assert(memcmp(a1->buffer, a2->buffer, sizeof(a1->buffer[0]) * PackedArray_bufferSize(a1)) == 0);

        PackedArray_unpack(a1, i, &v2, 1);
        assert(v1 == v2);
        PackedArray_unpack_reference(a2, i, &v2, 1);
        assert(v1 == v2);
        v2 = PackedArray_get(a2, i);
        assert(v1 == v2);
      }

      PackedArray_destroy(a1);
      PackedArray_destroy(a2);
    }
    printf("  %2d bits per item -- success.\n", bitsPerItem);
  }

  printf("\n");
  printf("bulk packing / unpacking:\n");
  for (bitsPerItem = 1; bitsPerItem <= 32; ++bitsPerItem)
  {
    uint32_t mask = (uint32_t)(1ULL << bitsPerItem) - 1;
    int count;

    for (count = 1; count <= 128; ++count)
    {
      uint32_t* v1;
      uint32_t* v2;
      PackedArray* a1;
      PackedArray* a2;
      int i, j;

      v1 = (uint32_t*)malloc(sizeof(uint32_t) * count);
      assert(v1 != NULL);
      v2 = (uint32_t*)malloc(sizeof(uint32_t) * count);
      assert(v2 != NULL);

      a1 = PackedArray_create(bitsPerItem, count);
      assert(a1 != NULL);
      a2 = PackedArray_create(bitsPerItem, count);
      assert(a2 != NULL);

      for (i = 0; i < (int)PackedArray_bufferSize(a1); ++i)
        a1->buffer[i] = a2->buffer[i] = rand();

      for (i = 0; i < count; ++i)
        v1[i] = rand() & mask;

      assert(bitsPerItem >= PackedArray_computeBitsPerItem(v1, count));

      for (i = 0; i < count; ++i)
      {
        for (j = 1; j <= count - i; ++j)
        {
          PackedArray_pack(a1, i, v1, j);
          PackedArray_pack_reference(a2, i, v1, j);
          assert(memcmp(a1->buffer, a2->buffer, sizeof(a1->buffer[0]) * PackedArray_bufferSize(a1)) == 0);

          PackedArray_unpack(a1, i, v2, j);
          assert(bitsPerItem >= PackedArray_computeBitsPerItem(v2, j));
          assert(memcmp(v1, v2, j * sizeof(uint32_t)) == 0);
          PackedArray_unpack_reference(a2, i, v2, j);
          assert(memcmp(v1, v2, j * sizeof(uint32_t)) == 0);
        }
      }

      PackedArray_destroy(a1);
      PackedArray_destroy(a2);
      free(v1);
      free(v2);
    }
    printf("  %2d bits per item -- success.\n", bitsPerItem);
  }

  return 0;
}

#elif defined(PACKEDARRAY_SELF_BENCH)  // #if defined(PACKEDARRAY_SELF_TEST)

#ifndef NDEBUG
#error please define NDEBUG to inhibit asserts when compiling the benchmark
#endif

#include <stdio.h>
#include <string.h>
#include <float.h>

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) (((x) < (y)) ? (y) : (x))
#endif

#ifdef _MSC_VER
#pragma warning(push, 3)
#include <windows.h>
#pragma warning(pop)
static double getChronometerTime(void)
{
  LARGE_INTEGER frequency;
  LARGE_INTEGER t;

  QueryPerformanceFrequency(&frequency);
  QueryPerformanceCounter(&t);

  return (double)t.QuadPart / (double)frequency.QuadPart * 1000;
}

#else
#include <sys/time.h>

static double getChronometerTime()
{
  struct timeval now = { 0 };
  gettimeofday(&now, NULL);

  return (double)now.tv_sec + (double)now.tv_usec * 1e-6;
}
#endif

#define LOOP_COUNT 1000
static double bench_memcpy(uint32_t* in, uint32_t* out, uint32_t count)
{
  double start, end;
  uint32_t i;

  start = getChronometerTime();

  for (i = 0; i < LOOP_COUNT; ++i)
    memcpy(out, in, count * sizeof(uint32_t));

  end = getChronometerTime();

  return 1e6 * (end - start) / LOOP_COUNT;
}

static double bench_loopcpy(uint32_t* in, uint32_t* out, uint32_t count)
{
  double start, end;
  uint32_t i;

  start = getChronometerTime();

  for (i = 0; i < LOOP_COUNT; ++i)
  {
    uint32_t j;

    for (j = 0; j < count; ++j)
      out[j] = in[j];
  }

  end = getChronometerTime();

  return 1e6 * (end - start) / LOOP_COUNT;
}

static double bench_pack(uint32_t* in, PackedArray* out, uint32_t count)
{
  double start, end;
  int i;

  start = getChronometerTime();

  for (i = 0; i < LOOP_COUNT; ++i)
    PackedArray_pack(out, 0, in, count);

  end = getChronometerTime();

  return 1e6 * (end - start) / LOOP_COUNT;
}

static double bench_unpack(PackedArray* in, uint32_t* out, uint32_t count)
{
  double start, end;
  int i;

  start = getChronometerTime();

  for (i = 0; i < LOOP_COUNT; ++i)
    PackedArray_unpack(in, 0, out, count);

  end = getChronometerTime();

  return 1e6 * (end - start) / LOOP_COUNT;
}

#define MAX_ELEMENT_COUNT (1 << 18)
#define LOG2_MAX_ELEMENT_COUNT 18
int main(void)
{
  double start, end;
  uint32_t* b1;
  uint32_t* b2;
  uint32_t count, bitsPerItem;
  PackedArray** packed;
  uint32_t i;
  double* speed_memcpy;
  double avg_memcpy, min_memcpy, max_memcpy;
  double* speed_loopcpy;
  double avg_loopcpy, min_loopcpy, max_loopcpy;
  double* speed_pack[32];
  double avg_pack, min_pack, max_pack;
  double* speed_unpack[32];
  double avg_unpack, min_unpack, max_unpack;

  printf("-- PackedArray self bench ------------------------------------------------------\n");

  start = getChronometerTime();

  b1 = (uint32_t*)malloc(sizeof(uint32_t) * MAX_ELEMENT_COUNT);
  assert(b1 != NULL);
  b2 = (uint32_t*)malloc(sizeof(uint32_t) * MAX_ELEMENT_COUNT);
  assert(b2 != NULL);

  packed = (PackedArray**)malloc(sizeof(PackedArray*) * 32);
  assert(packed != NULL);
  for (i = 0; i < 32; ++i)
    packed[i] = PackedArray_create(i + 1, MAX_ELEMENT_COUNT);

  for (i = 0; i < MAX_ELEMENT_COUNT; ++i)
    b1[i] = rand();

  speed_memcpy = (double*)malloc(sizeof(double) * (LOG2_MAX_ELEMENT_COUNT + 1));
  assert(speed_memcpy != NULL);
  avg_memcpy = 0;
  min_memcpy = DBL_MAX;
  max_memcpy = 0;

  printf("memcpy:\n");
  printf("bits\tsize (B)\ttime (µs)\tspeed (B/µs)\n");

  for (count = 1, i = 0; count <= MAX_ELEMENT_COUNT; count *= 2, ++i)
  {
    double elapsed = bench_memcpy(b1, b2, count);
    double speed = count * sizeof(uint32_t) / elapsed;
    printf("%4d\t%8d\t%9.3f\t%12.3f\n", 32, (uint32_t)(count * sizeof(uint32_t)), elapsed, speed);

    avg_memcpy += speed;
    min_memcpy = MIN(min_memcpy, speed);
    max_memcpy = MAX(max_memcpy, speed);

    speed_memcpy[i] = speed;
  }

  avg_memcpy /= i;

  printf("avg (B/µs)\tmin (B/µs)\tmax (B/µs)\n");
  printf("%10.3f\t%10.3f\t%10.3f\n", avg_memcpy, min_memcpy, max_memcpy);
  printf("\n");

  speed_loopcpy = (double*)malloc(sizeof(double) * (LOG2_MAX_ELEMENT_COUNT + 1));
  assert(speed_loopcpy != NULL);
  avg_loopcpy = 0;
  min_loopcpy = DBL_MAX;
  max_loopcpy = 0;

  printf("loopcpy:\n");
  printf("bits\tsize (B)\ttime (µs)\tspeed (B/µs)\n");

  for (count = 1, i = 0; count <= MAX_ELEMENT_COUNT; count *= 2, ++i)
  {
    double elapsed = bench_loopcpy(b1, b2, count);
    double speed = count * sizeof(uint32_t) / elapsed;
    printf("%4d\t%8d\t%9.3f\t%12.3f\n", 32, (uint32_t)(count * sizeof(uint32_t)), elapsed, speed);

    avg_loopcpy += speed;
    min_loopcpy = MIN(min_loopcpy, speed);
    max_loopcpy = MAX(max_loopcpy, speed);

    speed_loopcpy[i] = speed;
  }

  avg_loopcpy /= i;

  printf("avg (B/µs)\tmin (B/µs)\tmax (B/µs)\n");
  printf("%10.3f\t%10.3f\t%10.3f\n", avg_loopcpy, min_loopcpy, max_loopcpy);
  printf("\n");

  for (bitsPerItem = 1; bitsPerItem <= 32; ++bitsPerItem)
  {
    avg_pack = 0;
    min_pack = DBL_MAX;
    max_pack = 0;
    avg_unpack = 0;
    min_unpack = DBL_MAX;
    max_unpack = 0;

    printf("pack:\t        \t         \t            \t");
    printf("unpack:\t        \t         \t            \t");
    printf("\n");
    printf("bits\tsize (B)\ttime (µs)\tspeed (B/µs)");
    printf("\t");
    printf("bits\tsize (B)\ttime (µs)\tspeed (B/µs)");
    printf("\n");

    speed_pack[bitsPerItem - 1] = (double*)malloc(sizeof(double) * (LOG2_MAX_ELEMENT_COUNT + 1));
    assert(speed_pack[bitsPerItem - 1] != NULL);
    speed_unpack[bitsPerItem - 1] = (double*)malloc(sizeof(double) * (LOG2_MAX_ELEMENT_COUNT + 1));
    assert(speed_unpack[bitsPerItem - 1] != NULL);
    for (count = 1, i = 0; count <= MAX_ELEMENT_COUNT; count *= 2, ++i)
    {
      uint32_t mask = (uint32_t)(1ULL << bitsPerItem) - 1;
      uint32_t j;
      double elapsed, speed;

      for (j = 0; j < count; ++j)
        b2[j] = b1[j] & mask;

      elapsed = bench_pack(b2, packed[bitsPerItem - 1], count);
      speed = count * sizeof(uint32_t) / elapsed;
      printf("%4d\t%8d\t%9.3f\t%12.3f", bitsPerItem, (uint32_t)(count * sizeof(uint32_t)), elapsed, speed);

      avg_pack += speed;
      min_pack = MIN(min_pack, speed);
      max_pack = MAX(max_pack, speed);

      speed_pack[bitsPerItem - 1][i] = speed;

      printf("\t");

      elapsed = bench_unpack(packed[bitsPerItem - 1], b2, count);
      speed = count * sizeof(uint32_t) / elapsed;
      printf("%4d\t%8d\t%9.3f\t%12.3f", bitsPerItem, (uint32_t)(count * sizeof(uint32_t)), elapsed, speed);

      avg_unpack += speed;
      min_unpack = MIN(min_unpack, speed);
      max_unpack = MAX(max_unpack, speed);

      speed_unpack[bitsPerItem - 1][i] = speed;

      printf("\n");
    }
    printf("%d\n",i);
    assert(i == LOG2_MAX_ELEMENT_COUNT + 1);

    printf("avg (B/µs)\tmin (B/µs)\tmax (B/µs)");
    printf("\t\t");
    printf("avg (B/µs)\tmin (B/µs)\tmax (B/µs)");
    printf("\n");

    avg_pack /= i;
    printf("%10.3f\t%10.3f\t%10.3f", avg_pack, min_pack, max_pack);

    printf("\t\t");

    avg_unpack /= i;
    printf("%10.3f\t%10.3f\t%10.3f", avg_unpack, min_unpack, max_unpack);
    printf("\n");
    printf("\n");
  }

  printf("\n");

  printf("stats by bits per item\n");
  printf("pack:\t          \t          \t          \t");
  printf("unpack:\t          \t          \t          \t");
  printf("\n");
  printf("bits\tavg (B/µs)\tmin (B/µs)\tmax (B/µs)");
  printf("\t");
  printf("avg (B/µs)\tmin (B/µs)\tmax (B/µs)");
  printf("\n");
  for (bitsPerItem = 1; bitsPerItem <= 32; ++bitsPerItem)
  {
    avg_pack = 0;
    min_pack = DBL_MAX;
    max_pack = 0;
    avg_unpack = 0;
    min_unpack = DBL_MAX;
    max_unpack = 0;

    for (count = 1, i = 0; count <= MAX_ELEMENT_COUNT; count *= 2, ++i)
    {
      double speed;

      speed = speed_pack[bitsPerItem - 1][i];

      avg_pack += speed;
      min_pack = MIN(min_pack, speed);
      max_pack = MAX(max_pack, speed);

      speed = speed_unpack[bitsPerItem - 1][i];

      avg_unpack += speed;
      min_unpack = MIN(min_unpack, speed);
      max_unpack = MAX(max_unpack, speed);
    }
    assert(i == LOG2_MAX_ELEMENT_COUNT + 1);

    printf("%4d\t", bitsPerItem);

    avg_pack /= i;
    printf("%10.3f\t%10.3f\t%10.3f", avg_pack, min_pack, max_pack);
    printf("\t");

    avg_unpack /= i;
    printf("%10.3f\t%10.3f\t%10.3f", avg_unpack, min_unpack, max_unpack);
    printf("\n");
  }

  printf("\n");

  printf("stats by size\n");
  printf("pack:\t          \t          \t          \t");
  printf("unpack:\t          \t          \t          \t");
  printf("\n");
  printf("size(B)\tavg (B/µs)\tmin (B/µs)\tmax (B/µs)");
  printf("\t");
  printf("avg (B/µs)\tmin (B/µs)\tmax (B/µs)");
  printf("\n");
  for (count = 1, i = 0; count <= MAX_ELEMENT_COUNT; count *= 2, ++i)
  {
    avg_pack = 0;
    min_pack = DBL_MAX;
    max_pack = 0;
    avg_unpack = 0;
    min_unpack = DBL_MAX;
    max_unpack = 0;

    for (bitsPerItem = 1; bitsPerItem <= 32; ++bitsPerItem)
    {
      double speed;

      speed = speed_pack[bitsPerItem - 1][i];

      avg_pack += speed;
      min_pack = MIN(min_pack, speed);
      max_pack = MAX(max_pack, speed);

      speed = speed_unpack[bitsPerItem - 1][i];

      avg_unpack += speed;
      min_unpack = MIN(min_unpack, speed);
      max_unpack = MAX(max_unpack, speed);
    }

    printf("%7d\t", (uint32_t)sizeof(uint32_t) * count);

    avg_pack /= 32;
    printf("%10.3f\t%10.3f\t%10.3f", avg_pack, min_pack, max_pack);
    printf("\t");

    avg_unpack /= 32;
    printf("%10.3f\t%10.3f\t%10.3f", avg_unpack, min_unpack, max_unpack);
    printf("\n");
  }

  printf("\n");

  free(b1);
  free(b2);
  free(speed_memcpy);
  free(speed_loopcpy);

  for (i = 0; i < 32; ++i)
  {
    PackedArray_destroy(packed[i]);
    free(speed_pack[i]);
    free(speed_unpack[i]);
  }

  free(packed);

  end = getChronometerTime();
  printf("total time (s): %f\n", (end - start));
  printf("\n");

  return 0;
}

#endif // #elif defined(PACKEDARRAY_SELF_BENCH)

#endif // #ifdef PACKEDARRAY_IMPL
