/*
** autogenerated content - DO NOT EDIT
*/
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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include <time.h>
#include <sys/time.h>

#include "core/assert_api.h"
#include "core/bitpackstring.h"
#include "core/error.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/yarandom.h"

enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS_uint8_t = 100000,
};

static inline int
icmp_uint8_t(uint8_t a, uint8_t b)
{
  if (a > b)
    return 1;
  else if (a < b)
    return -1;
  else /* if (a == b) */
    return 0;
}

/**
 * \brief bit count reference
 * @param v count the number of bits set in v
 */
static inline int
genBitCount_uint8_t(uint8_t v)
{
  unsigned c; /* c accumulates the total bits set in v */
  for (c = 0; v; c++)
    v &= v - 1; /* clear the least significant bit set */
  return c;
}

#define freeResourcesAndReturn(retval) \
  do {                                 \
    gt_free(numBitsList);              \
    gt_free(randSrc);                  \
    gt_free(randCmp);                  \
    gt_free(bitStore);                 \
    gt_free(bitStoreCopy);             \
    return retval;                     \
  } while (0)

int
gt_bitPackStringInt8_unit_test(GtError *err)
{
  BitString bitStore = NULL;
  BitString bitStoreCopy = NULL;
  uint8_t *randSrc = NULL; /*< create random ints here for input as bit
                                *  store */
  uint8_t *randCmp = NULL; /*< used for random ints read back */
  unsigned *numBitsList = NULL;
  size_t i, numRnd;
  BitOffset offsetStart, offset;
  int had_err = 0;
  offset = offsetStart = random()%(sizeof (uint8_t) * CHAR_BIT);
  numRnd = random() % (MAX_RND_NUMS_uint8_t + 1);
  gt_log_log("offset=%lu, numRnd=%lu\n",
          (long unsigned)offsetStart, (long unsigned)numRnd);
  {
    BitOffset numBits = sizeof (uint8_t) * CHAR_BIT * numRnd + offsetStart;
    randSrc = gt_malloc(sizeof (uint8_t)*numRnd);
    bitStore = gt_malloc(bitElemsAllocSize(numBits) * sizeof (BitElem));
    bitStoreCopy = gt_calloc(bitElemsAllocSize(numBits), sizeof (BitElem));
    randCmp = gt_malloc(sizeof (uint8_t)*numRnd);
  }
  /* first test unsigned types */
  gt_log_log("gt_bsStoreUInt8/gt_bsGetUInt8: ");
  for (i = 0; i < numRnd; ++i)
  {
#if 8 > 32 && LONG_BIT < 8
    uint8_t v = randSrc[i] = (uint8_t)random() << 32 | random();
#else /* 8 > 32 && LONG_BIT < 8 */
    uint8_t v = randSrc[i] = random();
#endif /* 8 > 32 && LONG_BIT < 8 */
    int bits = gt_requiredUInt8Bits(v);
    gt_bsStoreUInt8(bitStore, offset, bits, v);
    offset += bits;
  }
  offset = offsetStart;
  for (i = 0; i < numRnd; ++i)
  {
    uint8_t v = randSrc[i];
    int bits = gt_requiredUInt8Bits(v);
    uint8_t r = gt_bsGetUInt8(bitStore, offset, bits);
    gt_ensure(r == v);
    if (had_err)
    {
      gt_log_log("Expected %"PRIu8", got %"PRIu8", i = %lu\n",
              v, r, (unsigned long)i);
      freeResourcesAndReturn(had_err);
    }
    offset += bits;
  }
  gt_log_log("passed\n");
  if (numRnd > 0)
  {
    uint8_t v = randSrc[0], r = 0;
    unsigned numBits = gt_requiredUInt8Bits(v);
    BitOffset i = offsetStart + numBits;
    uint8_t mask = ~(uint8_t)0;
    if (numBits < 8)
      mask = ~(mask << numBits);
    gt_log_log("bsSetBit, gt_bsClearBit, bsToggleBit, gt_bsGetBit: ");
    while (v)
    {
      int lowBit = v & 1;
      v >>= 1;
      gt_ensure(lowBit == (r = gt_bsGetBit(bitStore, --i)));
      if (had_err)
      {
        gt_log_log("Expected %d, got %d, i = %llu\n",
                lowBit, (int)r, (unsigned long long)i);
        freeResourcesAndReturn(had_err);
      }
    }
    i = offsetStart + numBits;
    gt_bsClear(bitStoreCopy, offsetStart, numBits, random()&1);
    v = randSrc[0];
    while (i)
    {
      int lowBit = v & 1;
      v >>= 1;
      if (lowBit)
        bsSetBit(bitStoreCopy, --i);
      else
        gt_bsClearBit(bitStoreCopy, --i);
    }
    v = randSrc[0];
    r = gt_bsGetUInt8(bitStoreCopy, offsetStart, numBits);
    gt_ensure(r == v);
    if (had_err)
    {
      gt_log_log("Expected %"PRIu8", got %"PRIu8"\n", v, r);
      freeResourcesAndReturn(had_err);
    }
    for (i = 0; i < numBits; ++i)
      bsToggleBit(bitStoreCopy, offsetStart + i);
    r = gt_bsGetUInt8(bitStoreCopy, offsetStart, numBits);
    gt_ensure(r == (v = (~v & mask)));
    if (had_err)
    {
      gt_log_log("Expected %"PRIu8", got %"PRIu8"\n", v, r);
      freeResourcesAndReturn(had_err);
    }
    gt_log_log("passed\n");
  }
  if (numRnd > 1)
  {
    gt_log_log("gt_bsCompare: ");
    {
      uint8_t v0 = randSrc[0];
      int bits0 = gt_requiredUInt8Bits(v0);
      uint8_t r0;
      offset = offsetStart;
      r0 = gt_bsGetUInt8(bitStore, offset, bits0);
      for (i = 1; i < numRnd; ++i)
      {
        uint8_t v1 = randSrc[i];
        int bits1 = gt_requiredUInt8Bits(v1);
        uint8_t r1 = gt_bsGetUInt8(bitStore, offset + bits0, bits1);
        int result = -2;   /*< -2 is not a return value of gt_bsCompare, thus
                            *   if it is displayed, there was an earlier
                            *   error. */
        gt_ensure(r0 == v0 && r1 == v1);
        gt_ensure(icmp_uint8_t(v0, v1) ==
               (result = gt_bsCompare(bitStore, offset, bits0,
                                   bitStore, offset + bits0, bits1)));
        if (had_err)
        {
          gt_log_log("Expected v0 %s v1, got v0 %s v1,\n for v0=%"
                  PRIu8" and v1=%"PRIu8",\n"
                  "i = %lu, bits0=%u, bits1=%u\n",
                  (v0 > v1?">":(v0 < v1?"<":"==")),
                  (result > 0?">":(result < 0?"<":"==")), v0, v1,
                  (unsigned long)i, bits0, bits1);
          freeResourcesAndReturn(had_err);
        }
        offset += bits0;
        bits0 = bits1;
        v0 = v1;
        r0 = r1;
      }
    }
    gt_log_log("passed\n");
  }
  gt_log_log("gt_bsStoreUniformUInt8Array/gt_bsGetUInt8: ");
  if (numRnd > 0) {
    unsigned numBits = random()%8 + 1;
    uint8_t mask = ~(uint8_t)0;
    if (numBits < 8)
      mask = ~(mask << numBits);
    offset = offsetStart;
    gt_bsStoreUniformUInt8Array(bitStore, offset, numBits, numRnd,
                                     randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      uint8_t v = randSrc[i] & mask;
      uint8_t r = gt_bsGetUInt8(bitStore, offset, numBits);
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRIu8", got %"PRIu8",\n"
                "i = %lu, bits=%u\n", v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    gt_log_log("passed\n");
    gt_log_log("gt_bsStoreUniformUInt8Array/"
               "gt_bsGetUniformUInt8Array: ");
    gt_bsGetUniformUInt8Array(bitStore, offset = offsetStart,
                               numBits, numRnd, randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      uint8_t v = randSrc[i] & mask;
      uint8_t r = randCmp[i];
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log( "Expected %"PRIu8", got %"PRIu8",\n"
                " i = %lu, bits=%u\n",
                v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 1)
    {
      uint8_t v = randSrc[0] & mask;
      uint8_t r = 0;
      gt_bsGetUniformUInt8Array(bitStore, offsetStart,
                                 numBits, 1, &r);
      if (r != v)
      {
        gt_log_log("Expected %"PRIu8", got %"PRIu8","
                " one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    gt_log_log(" passed\n");
  }
  /* int types */
  gt_log_log("gt_bsStoreInt8/gt_bsGetInt8: ");
  for (i = 0; i < numRnd; ++i)
  {
    int8_t v = (int8_t)randSrc[i];
    unsigned bits = gt_requiredInt8Bits(v);
    gt_bsStoreInt8(bitStore, offset, bits, v);
    offset += bits;
  }
  offset = offsetStart;
  for (i = 0; i < numRnd; ++i)
  {
    int8_t v = randSrc[i];
    unsigned bits = gt_requiredInt8Bits(v);
    int8_t r = gt_bsGetInt8(bitStore, offset, bits);
    gt_ensure(r == v);
    if (had_err)
    {
      gt_log_log("Expected %"PRId8", got %"PRId8",\n"
                  "i = %lu, bits=%u\n",
                  v, r, (unsigned long)i, bits);
      freeResourcesAndReturn(had_err);
    }
    offset += bits;
  }
  gt_log_log("passed\n");
  gt_log_log("gt_bsStoreUniformInt8Array/gt_bsGetInt8: ");
  if (numRnd > 0) {
    unsigned numBits = random()%8 + 1;
    int8_t mask = ~(int8_t)0;
    if (numBits < 8)
      mask = ~(mask << numBits);
    offset = offsetStart;
    gt_bsStoreUniformInt8Array(bitStore, offset, numBits, numRnd,
                                (int8_t *)randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      int8_t m = (int8_t)1 << (numBits - 1);
      int8_t v = (int8_t)((randSrc[i] & mask) ^ m) - m;
      int8_t r = gt_bsGetInt8(bitStore, offset, numBits);
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRId8", got %"PRId8",\n"
                    "i = %lu, numBits=%u\n",
                    v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    gt_log_log("passed\n");
    gt_log_log("gt_bsStoreUniformInt8Array/"
               "gt_bsGetUniformInt8Array: ");
    gt_bsGetUniformInt8Array(bitStore, offset = offsetStart,
                              numBits, numRnd, (int8_t *)randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      int8_t m = (int8_t)1 << (numBits - 1);
      int8_t v = (int8_t)((randSrc[i] & mask) ^ m) - m;
      int8_t r = randCmp[i];
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRId8", got %"PRId8", i = %lu\n",
                v, r, (unsigned long)i);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 0)
    {
      int8_t m = (int8_t)1 << (numBits - 1);
      int8_t v = (int8_t)((randSrc[0] & mask) ^ m) - m;
      int8_t r = 0;
      gt_bsGetUniformInt8Array(bitStore, offsetStart,
                                numBits, 1, &r);
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRId8", got %"PRId8
                ", one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    gt_log_log("passed\n");
  }

  gt_log_log("gt_bsStoreNonUniformUInt8Array/gt_bsGetUInt8: ");
  if (numRnd != 0) {
    BitOffset bitsTotal = 0;
    numBitsList = gt_malloc(sizeof (unsigned) * numRnd);
    for (i = 0; i < numRnd; ++i)
      bitsTotal += (numBitsList[i] = random()%8 + 1);
    offset = offsetStart;
    gt_bsStoreNonUniformUInt8Array(bitStore, offset, numRnd, bitsTotal,
                                     numBitsList, randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      uint8_t mask = (numBits < 8)?
        ~((~(uint8_t)0) << numBits):~(uint8_t)0;
      uint8_t v = randSrc[i] & mask;
      uint8_t r = gt_bsGetUInt8(bitStore, offset, numBits);
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRIu8", got %"PRIu8",\n"
                "i = %lu, bits=%u\n",
                v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    gt_log_log("passed\n");
    gt_log_log("gt_bsStoreNonUniformUInt8Array/"
            "gt_bsGetNonUniformUInt8Array: ");
    gt_bsGetNonUniformUInt8Array(bitStore, offset = offsetStart,
                                   numRnd, bitsTotal, numBitsList, randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      uint8_t mask = (numBits < 8)?
        ~((~(uint8_t)0) << numBits):~(uint8_t)0;
      uint8_t v = randSrc[i] & mask,
        r = randCmp[i];
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log( "Expected %"PRIu8", got %"PRIu8",\n"
                " i = %lu, bits=%u\n",
                v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 1)
    {
      unsigned numBits = numBitsList[0];
      uint8_t mask = (numBits < 8)?
        ~((~(uint8_t)0) << numBits):~(uint8_t)0;
      uint8_t v = randSrc[0] & mask;
      uint8_t r = 0;
      gt_bsGetNonUniformUInt8Array(bitStore, offsetStart, 1, numBits,
                                     numBitsList, &r);
      if (r != v)
      {
        gt_log_log("Expected %"PRIu8", got %"PRIu8", "
                " one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    gt_log_log(" passed\n");
    gt_free(numBitsList);
    numBitsList = NULL;
  }
  gt_log_log("bsNonStoreUniformInt8Array/gt_bsGetInt8: ");
  if (numRnd != 0) {
    BitOffset bitsTotal = 0;
    numBitsList = gt_malloc(sizeof (unsigned) * numRnd);
    for (i = 0; i < numRnd; ++i)
      bitsTotal += (numBitsList[i] = random()%8 + 1);
    offset = offsetStart;
    gt_bsStoreNonUniformInt8Array(bitStore, offset, numRnd, bitsTotal,
                                     numBitsList, (int8_t *)randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      int8_t mask = (numBits < 8)
        ? ~((~(int8_t)0) << numBits) : ~(int8_t)0;
      int8_t m = (int8_t)1 << (numBits - 1);
      int8_t v = (int8_t)((randSrc[i] & mask) ^ m) - m;
      int8_t r = gt_bsGetInt8(bitStore, offset, numBits);
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRId8", got %"PRId8",\n"
                    "i = %lu, numBits=%u\n",
                    v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    gt_log_log("passed\n");
    gt_log_log("gt_bsStoreNonUniformInt8Array/"
            "gt_bsGetNonUniformInt8Array: ");
    gt_bsGetNonUniformInt8Array(bitStore, offset = offsetStart, numRnd,
                                   bitsTotal, numBitsList,
                                   (int8_t *)randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      int8_t mask = (numBits < 8)
        ? ~((~(int8_t)0) << numBits) : ~(int8_t)0;
      int8_t m = (int8_t)1 << (numBits - 1);
      int8_t v = (int8_t)((randSrc[i] & mask) ^ m) - m;
      int8_t r = randCmp[i];
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRId8", got %"PRId8", i = %lu\n",
                v, r, (unsigned long)i);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 0)
    {
      unsigned numBits = numBitsList[0];
      int8_t mask = (numBits < 8)
        ? ~((~(int8_t)0) << numBits) : ~(int8_t)0;
      int8_t m = (int8_t)1 << (numBits - 1);
      int8_t v = (int8_t)((randSrc[0] & mask) ^ m) - m;
      int8_t r = 0;
      gt_bsGetNonUniformInt8Array(bitStore, offsetStart,
                                     1, numBits, numBitsList, &r);
      gt_ensure(r == v);
      if (had_err)
      {
        gt_log_log("Expected %"PRId8", got %"PRId8
                ", one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    gt_log_log("passed\n");
    gt_free(numBitsList);
    numBitsList = NULL;
  }

  if (numRnd > 0)
  {
    gt_log_log("gt_bsCopy: ");
    {
      /* first decide how many of the values to use and at which to start */
      size_t numValueCopies, copyStart;
      BitOffset numCopyBits = 0, destOffset;
      unsigned numBits = random()%8 + 1;
      uint8_t mask = ~(uint8_t)0;
      if (numBits < 8)
        mask = ~(mask << numBits);
      if (random()&1)
      {
        numValueCopies = random()%(numRnd + 1);
        copyStart = random()%(numRnd - numValueCopies + 1);
      }
      else
      {
        copyStart = random() % numRnd;
        numValueCopies = random()%(numRnd - copyStart) + 1;
      }
      gt_assert(copyStart + numValueCopies <= numRnd);
      offset = offsetStart + (BitOffset)copyStart * numBits;
      gt_bsStoreUniformUInt8Array(bitStore, offset, numBits,
                                       numValueCopies, randSrc);
      destOffset = random()%(offsetStart + 8
                             * (BitOffset)(numRnd - numValueCopies) + 1);
      numCopyBits = (BitOffset)numBits * numValueCopies;
      /* the following gt_bsCopy should be equivalent to:
       * gt_bsStoreUniformUInt8Array(bitStoreCopy, destOffset,
       *                              numBits, numValueCopies, randSrc); */
      gt_bsCopy(bitStore, offset, bitStoreCopy, destOffset, numCopyBits);
      gt_ensure(
                gt_bsCompare(bitStore, offset, numCopyBits,
                             bitStoreCopy, destOffset, numCopyBits) == 0);
      if (had_err)
      {
        gt_log_log("Expected equality on bitstrings\n"
                    "offset = %llu, destOffset = %llu,"
                    " numCopyBits=%llu\n",
                    (unsigned long long)offset,
                    (unsigned long long)destOffset,
                    (unsigned long long)numCopyBits);
        /* FIXME: implement bitstring output function */
        freeResourcesAndReturn(had_err);
      }
      gt_log_log("passed\n");
    }
  }
  if (numRnd > 0)
  {
    gt_log_log("gt_bsClear: ");
    {
      /* first decide how many of the values to use and at which to start */
      size_t numResetValues, resetStart;
      BitOffset numResetBits = 0;
      unsigned numBits = random()%8 + 1;
      int bitVal = random()&1;
      int8_t cmpVal = bitVal?-1:0;
      uint8_t mask = ~(uint8_t)0;
      if (numBits < 8)
        mask = ~(mask << numBits);
      if (random()&1)
      {
        numResetValues = random()%(numRnd + 1);
        resetStart = random()%(numRnd - numResetValues + 1);
      }
      else
      {
        resetStart = random() % numRnd;
        numResetValues = random()%(numRnd - resetStart) + 1;
      }
      gt_assert(resetStart + numResetValues <= numRnd);
      offset = offsetStart;
      gt_bsStoreUniformInt8Array(bitStore, offset, numBits, numRnd,
                                    (int8_t *)randSrc);
      numResetBits = (BitOffset)numBits * numResetValues;
      gt_bsClear(bitStore, offset + (BitOffset)resetStart * numBits,
              numResetBits, bitVal);
      {
        int8_t m = (int8_t)1 << (numBits - 1);
        for (i = 0; i < resetStart; ++i)
        {
          int8_t v = (int8_t)((randSrc[i] & mask) ^ m) - m;
          int8_t r = gt_bsGetInt8(bitStore, offset, numBits);
          gt_ensure(r == v);
          if (had_err)
          {
            gt_log_log( "Expected %"PRId8", got %"PRId8",\n"
                     "i = %lu, numBits=%u\n",
                     v, r, (unsigned long)i, numBits);
            freeResourcesAndReturn(had_err);
          }
          offset += numBits;
        }
        for (; i < resetStart + numResetValues; ++i)
        {
          int8_t r = gt_bsGetInt8(bitStore, offset, numBits);
          gt_ensure(r == cmpVal);
          if (had_err)
          {
            gt_log_log("Expected %"PRId8", got %"PRId8",\n"
                    "i = %lu, numBits=%u\n",
                    cmpVal, r, (unsigned long)i, numBits);
            freeResourcesAndReturn(had_err);
          }
          offset += numBits;
        }
        for (; i < numRnd; ++i)
        {
          int8_t v = (int8_t)((randSrc[i] & mask) ^ m) - m;
          int8_t r = gt_bsGetInt8(bitStore, offset, numBits);
          gt_ensure(r == v);
          if (had_err)
          {
            gt_log_log("Expected %"PRId8", got %"PRId8",\n"
                    "i = %lu, numBits=%u\n",
                    v, r, (unsigned long)i, numBits);
            freeResourcesAndReturn(had_err);
          }
          offset += numBits;
        }
      }
    }
    gt_log_log("passed\n");
  }
  if (numRnd > 0)
  {
    gt_log_log("gt_bs1BitsCount: ");
    {
      /* first decide how many of the values to use and at which to start */
      size_t numCountValues, countStart;
      BitOffset numCountBits = 0, bitCountRef = 0, bitCountCmp;
      unsigned numBits = random()%8 + 1;
      uint8_t mask = ~(uint8_t)0;
      if (numBits < 8)
        mask = ~(mask << numBits);
      if (random()&1)
      {
        numCountValues = random()%(numRnd + 1);
        countStart = random()%(numRnd - numCountValues + 1);
      }
      else
      {
        countStart = random() % numRnd;
        numCountValues = random()%(numRnd - countStart) + 1;
      }
      gt_assert(countStart + numCountValues <= numRnd);
      offset = offsetStart;
      gt_bsStoreUniformUInt8Array(bitStore, offset, numBits, numRnd,
                                       randSrc);
      numCountBits = (BitOffset)numBits * numCountValues;
      bitCountCmp = gt_bs1BitsCount(bitStore,
                                 offset + (BitOffset)countStart * numBits,
                                 numCountBits);
      for (i = countStart; i < countStart + numCountValues; ++i)
      {
        uint8_t v = (uint8_t)randSrc[i] & mask;
        bitCountRef += genBitCount_uint8_t(v);
      }
      gt_ensure(bitCountRef == bitCountCmp);
      if (had_err)
      {
        gt_log_log("Expected %llu, got %llu,\n"
                "numBits=%u\n", (unsigned long long)bitCountRef,
                (unsigned long long)bitCountCmp, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    gt_log_log("passed\n");
  }
  freeResourcesAndReturn(had_err);
}
/*
vim: ft=c
 */
