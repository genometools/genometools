/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
   copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef BITPACKSTRINGSIMPLEOP_H
#define BITPACKSTRINGSIMPLEOP_H

#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "core/assert_api.h"
#include "core/bitpackstring.h"
#include "core/dynalloc.h"
#include "core/safearith.h"

/**
 * \if INTERNAL \file bitpackstringsimpleop.h \endif
 * \brief Trivial (i.e. inlined) operations on bitstrings.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>

 * The contents of this file is to be considered private
 * implementation detail but is exposed to the compiler solely for
 * performance reasons. Just pretend you didn't know how these
 * functions work.
 */

static inline size_t
bitElemsAllocSize(BitOffset numBits)
{
  BitOffset size = numBits/bitElemBits + ((numBits%bitElemBits)?1:0);
  return ((size * sizeof (BitElem) > SIZE_MAX)?
          SIZE_MAX / sizeof (BitElem) : size);
}

static inline int8_t
gt_bsGetInt8(constBitString str, BitOffset offset, unsigned numBits)
{
  /* requires sign extension */
  int8_t m = 1 << (numBits - 1);
  return (gt_bsGetUInt8(str, offset, numBits) ^ m) - m;
}

static inline int16_t
gt_bsGetInt16(constBitString str, BitOffset offset, unsigned numBits)
{
  /* requires sign extension */
  int16_t m = 1 << (numBits - 1);
  return (gt_bsGetUInt16(str, offset, numBits) ^ m) - m;
}

static inline int32_t
gt_bsGetInt32(constBitString str, BitOffset offset, unsigned numBits)
{
  /* requires sign extension */
  int32_t m = 1 << (numBits - 1);
  return (gt_bsGetUInt32(str, offset, numBits) ^ m) - m;
}

static inline int64_t
gt_bsGetInt64(constBitString str, BitOffset offset, unsigned numBits)
{
  /* requires sign extension */
  int64_t m = (int64_t)1 << (numBits - 1);
  return (gt_bsGetUInt64(str, offset, numBits) ^ m) - m;
}

static inline void
gt_bsStoreInt8(BitString str, BitOffset offset, unsigned numBits,
            uint8_t val)
{
  gt_bsStoreUInt8(str, offset, numBits,  val);
}

static inline void
gt_bsStoreInt16(BitString str, BitOffset offset, unsigned numBits,
             uint16_t val)
{
  gt_bsStoreUInt16(str, offset, numBits,  val);
}

static inline void
gt_bsStoreInt32(BitString str, BitOffset offset, unsigned numBits, uint32_t val)
{
  gt_bsStoreUInt32(str, offset, numBits,  val);
}

static inline void
gt_bsStoreInt64(BitString str, BitOffset offset, unsigned numBits, uint64_t val)
{
  gt_bsStoreUInt64(str, offset, numBits,  val);
}

static inline int
gt_requiredUInt8Bits(uint8_t v)
{
  return gt_requiredUInt32Bits((uint32_t)v);
}

static inline int
gt_requiredUInt16Bits(uint16_t v)
{
  return gt_requiredUInt32Bits((uint32_t)v);
}

static inline int
gt_requiredInt8Bits(int8_t v)
{
  if (v == INT8_MIN)
    return sizeof (v)*CHAR_BIT;
  else
    /* one extra for sign bit */
    return gt_requiredUInt32Bits(gt_safe_abs(v)) + 1;
}

static inline int
gt_requiredInt16Bits(int16_t v)
{
  if (v == INT16_MIN)
    return sizeof (v)*CHAR_BIT;
  else
    /* one extra for sign bit */
    return gt_requiredUInt32Bits(gt_safe_abs(v)) + 1;
}

static inline int
gt_requiredInt32Bits(int32_t v)
{
  if (v == INT32_MIN)
    return sizeof (v)*CHAR_BIT;
  else
    /* one extra for sign bit */
    return gt_requiredUInt32Bits(gt_safe_labs(v)) + 1;
}

static inline int
gt_requiredInt64Bits(int64_t v)
{
  if (v == INT64_MIN)
    return sizeof (v)*CHAR_BIT;
  else
    /* one extra for sign bit */
    return gt_requiredUInt64Bits(gt_safe_llabs(v)) + 1;
}

static inline void
gt_bsStoreUniformInt8Array(BitString str, BitOffset offset, unsigned numBits,
                           size_t numValues, const int8_t val[])
{
  gt_bsStoreUniformUInt8Array(str, offset, numBits, numValues,
                              (const uint8_t *)val);
}

static inline void
gt_bsStoreUniformInt16Array(BitString str, BitOffset offset, unsigned numBits,
                            size_t numValues, const int16_t val[])
{
  gt_bsStoreUniformUInt16Array(str, offset, numBits, numValues,
                               (const uint16_t *)val);
}

static inline void
gt_bsStoreUniformInt32Array(BitString str, BitOffset offset, unsigned numBits,
                            size_t numValues, const int32_t val[])
{
  gt_bsStoreUniformUInt32Array(str, offset, numBits, numValues,
                               (const uint32_t *)val);
}

static inline void
gt_bsStoreUniformInt64Array(BitString str, BitOffset offset, unsigned numBits,
                            size_t numValues, const int64_t val[])
{
  gt_bsStoreUniformUInt64Array(str, offset, numBits, numValues,
                               (const uint64_t *)val);
}

static inline void
gt_bsStoreNonUniformInt8Array(BitString str, BitOffset offset, size_t numValues,
                              BitOffset totalBits, unsigned *numBitsList,
                              const int8_t val[])
{
  gt_bsStoreNonUniformUInt8Array(str, offset, numValues, totalBits, numBitsList,
                                 (const uint8_t *)val);
}

static inline void
gt_bsStoreNonUniformInt16Array(BitString str, BitOffset offset,
                               size_t numValues, BitOffset totalBits,
                               unsigned *numBitsList, const int16_t val[])
{
  gt_bsStoreNonUniformUInt16Array(str, offset, numValues, totalBits,
                                  numBitsList, (const uint16_t *)val);
}

static inline void
gt_bsStoreNonUniformInt32Array(BitString str, BitOffset offset,
                               size_t numValues, BitOffset totalBits,
                               unsigned *numBitsList, const int32_t val[])
{
  gt_bsStoreNonUniformUInt32Array(str, offset, numValues, totalBits,
                                  numBitsList, (const uint32_t *)val);
}

static inline void
gt_bsStoreNonUniformInt64Array(BitString str, BitOffset offset,
                               size_t numValues, BitOffset totalBits,
                               unsigned *numBitsList, const int64_t val[])
{
  gt_bsStoreNonUniformUInt64Array(str, offset, numValues, totalBits,
                                  numBitsList, (const uint64_t *)val);
}

static inline void
bsSignExpandArray8(int8_t val[], size_t numValues, unsigned numBits)
{
  int8_t m = 1 << (numBits - 1);
  size_t i;
  for (i = 0; i < numValues; ++i)
    val[i] = (val[i] ^ m) - m;
}

static inline void
bsSignExpandArray16(int16_t val[], size_t numValues, unsigned numBits)
{
  int16_t m = 1 << (numBits - 1);
  size_t i;
  for (i = 0; i < numValues; ++i)
    val[i] = (val[i] ^ m) - m;
}

static inline void
bsSignExpandArray32(int32_t val[], size_t numValues, unsigned numBits)
{
  int32_t m = ((int32_t)1) << (numBits - 1);
  size_t i;
  for (i = 0; i < numValues; ++i)
    val[i] = (val[i] ^ m) - m;
}

static inline void
bsSignExpandArray64(int64_t val[], size_t numValues, unsigned numBits)
{
  int64_t m = ((int64_t)1) << (numBits - 1);
  size_t i;
  for (i = 0; i < numValues; ++i)
    val[i] = (val[i] ^ m) - m;
}

static inline void
gt_bsGetUniformInt8Array(constBitString str, BitOffset offset, unsigned numBits,
                         size_t numValues, int8_t val[])
{
  /* read blocksize many ints at once, for later sign expansion */
  int blockSize = 16 * getpagesize() / sizeof (val[0]);
  /* the factor (16) is completely arbitrary and needs some profiling */
  int8_t *blockPtr = val;
  size_t revIndex = numValues;
  BitOffset offsetTemp = offset;
  while (revIndex > blockSize)
  {
    gt_bsGetUniformUInt8Array(str, offsetTemp, numBits,
                              blockSize, (uint8_t *)blockPtr);
    bsSignExpandArray8(blockPtr, blockSize, numBits);
    blockPtr += blockSize;
    revIndex -= blockSize;
    offsetTemp += blockSize * numBits;
  }
  gt_bsGetUniformUInt8Array(str, offsetTemp, numBits,
                            revIndex, (uint8_t *)blockPtr);
  bsSignExpandArray8(blockPtr, revIndex, numBits);
}

static inline void
gt_bsGetUniformInt16Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, int16_t val[])
{
  /* read blocksize many ints at once, for later sign expansion */
  int blockSize = 16 * getpagesize() / sizeof (val[0]);
  /* the factor (16) is completely arbitrary and needs some profiling */
  int16_t *blockPtr = val;
  size_t revIndex = numValues;
  BitOffset offsetTemp = offset;
  while (revIndex > blockSize)
  {
    gt_bsGetUniformUInt16Array(str, offsetTemp, numBits,
                               blockSize, (uint16_t *)blockPtr);
    bsSignExpandArray16(blockPtr, blockSize, numBits);
    blockPtr += blockSize;
    revIndex -= blockSize;
    offsetTemp += blockSize * numBits;
  }
  gt_bsGetUniformUInt16Array(str, offsetTemp, numBits,
                             revIndex, (uint16_t *)blockPtr);
  bsSignExpandArray16(blockPtr, revIndex, numBits);
}

static inline void
gt_bsGetUniformInt32Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, int32_t val[])
{
  /* read blocksize many ints at once, for later sign expansion */
  int blockSize = 16 * getpagesize() / sizeof (val[0]);
  /* the factor (16) is completely arbitrary and needs some profiling */
  int32_t *blockPtr = val;
  size_t revIndex = numValues;
  BitOffset offsetTemp = offset;
  while (revIndex > blockSize)
  {
    gt_bsGetUniformUInt32Array(str, offsetTemp, numBits,
                               blockSize, (uint32_t *)blockPtr);
    bsSignExpandArray32(blockPtr, blockSize, numBits);
    blockPtr += blockSize;
    revIndex -= blockSize;
    offsetTemp += blockSize * numBits;
  }
  gt_bsGetUniformUInt32Array(str, offsetTemp, numBits,
                             revIndex, (uint32_t *)blockPtr);
  bsSignExpandArray32(blockPtr, revIndex, numBits);
}

static inline void
gt_bsGetUniformInt64Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, int64_t val[])
{
  /* read blocksize many ints at once, for later sign expansion */
  int blockSize = 16 * getpagesize() / sizeof (val[0]);
  /* the factor (16) is completely arbitrary and needs some profiling */
  int64_t *blockPtr = val;
  size_t revIndex = numValues;
  BitOffset offsetTemp = offset;
  while (revIndex > blockSize)
  {
    gt_bsGetUniformUInt64Array(str, offsetTemp, numBits,
                               blockSize, (uint64_t *)blockPtr);
    bsSignExpandArray64(blockPtr, blockSize, numBits);
    blockPtr += blockSize;
    revIndex -= blockSize;
    offsetTemp += blockSize * numBits;
  }
  gt_bsGetUniformUInt64Array(str, offsetTemp, numBits,
                             revIndex, (uint64_t *)blockPtr);
  bsSignExpandArray64(blockPtr, revIndex, numBits);
}

static inline void
bsSetBit(BitString str, BitOffset pos)
{
  unsigned bitTop = pos%bitElemBits;
  BitElem *p = str + pos/bitElemBits;
  gt_assert(str);
  *p |= 1 << (bitElemBits - bitTop - 1);
}

static inline void
gt_bsClearBit(BitString str, BitOffset pos)
{
  unsigned bitTop = pos%bitElemBits;
  BitElem *p = str + pos/bitElemBits;
  gt_assert(str);
  *p &= ~(1 << (bitElemBits - bitTop - 1));
}

static inline void
bsToggleBit(BitString str, BitOffset pos)
{
  unsigned bitTop = pos%bitElemBits;
  BitElem *p = str + pos/bitElemBits;
  gt_assert(str);
  *p ^= 1 << (bitElemBits - bitTop - 1);
}

static inline int
gt_bsGetBit(constBitString str, BitOffset pos)
{
  unsigned bitTop = pos%bitElemBits;
  const BitElem *p = str + pos/bitElemBits;
  gt_assert(str);
  return (*p & 1 << (bitElemBits - bitTop - 1))?1:0;
}

/*
 * the following helper function counts the bits set in a 32 bit int
 * see http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
 * for an exhaustive discussion
 */
static inline unsigned
bitCountUInt32(uint32_t v)
{
  v = v - ((v >> 1) & (uint32_t)0x55555555UL);   /* reuse input as temporary */
  v = (v & (uint32_t)0x33333333UL) +
      ((v >> 2) & (uint32_t)0x33333333UL);       /* temp  */
  return (((v + (v >> 4)) & (uint32_t)0xF0F0F0FUL) *
         (uint32_t)0x1010101UL) >> 24;           /* count */
}

#endif
