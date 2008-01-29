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

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "libgtcore/bitpackstring.h"
#include "libgtcore/minmax.h"
/**
 * \if INTERNAL \file bitpackstringop.c \endif
 * Involved (i.e. not inlined) operations on bitstrings.
 */

/*
 * Both requiredUInt{32|64}Bits functions are based on the concepts
 * presented at
 * http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogDeBruijn
 * the method has two steps:
 * 1. isolate the highest bit set by first copying the highest bit set
 * via the shift and or instructions, then the final divide and
 * increment result in the highest bit being set.
 * 2. lookup the 5/6 top bits resulting from multiplication with a
 * DeBruijn bit sequence (the long unsigned constant), since a
 * DeBruijn sequence has all q-words differ by at least one bit, any
 * bit set in v results in a corresponding table lookup.
 */
int
requiredUInt32Bits(uint32_t v)
{
  int r;
  static const int MultiplyDeBruijnBitPosition[32] = {
    1, 2, 29, 3, 30, 15, 25, 4, 31, 23, 21, 16, 26, 18, 5, 9,
    32, 28, 14, 24, 22, 20, 17, 8, 27, 13, 19, 7, 12, 6, 11, 10
  };
  v |= v >> 1; /* first round down to power of 2 */
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v = (v >> 1) + 1;
  r = MultiplyDeBruijnBitPosition[(v * (uint32_t)0x077CB531UL) >> 27];
  return r;
}

int
requiredUInt64Bits(uint64_t v)
{
  int r;
  static const int MultiplyDeBruijnBitPosition[64] = {
    1, 2, 3, 57, 4, 33, 58, 47, 30, 5, 21, 34, 8, 59, 12, 48,
    63, 31, 19, 6, 17, 22, 35, 24, 54, 9, 60, 37, 26, 13, 49, 40,
    64, 56, 32, 46, 29, 20, 7, 11, 62, 18, 16, 23, 53, 36, 25, 39,
    55, 45, 28, 10, 61, 15, 52, 38, 44, 27, 14, 51, 43, 50, 42, 41
  };
  v |= v >> 1; /* first round down to power of 2 */
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v = (v >> 1) + 1;
  r = MultiplyDeBruijnBitPosition[(v * (uint64_t)0x26752B916FC7B0DULL) >> 58];
  return r;
}

#if 0
/**
 * This function was used to compute the 64-bit DeBruijn table shown
 * above. The multiplicator can be gathered from diverse web sources.
 */
static void
computeDeBruijn()
{
  unsigned prod, i;
  unsigned long long v;
  int MultiplyDeBruijnBitPosition[64];
  for (v = 1, i = 1; i <= 64; ++i, v <<= 1)
  {
    prod = ((unsigned long long)v * 0x26752B916FC7B0DULL) >> 58;
    printf("v = %llu, i = %u, prod = %u\n", v, i, prod);
    assert(prod < 64);
    MultiplyDeBruijnBitPosition[prod] = i;
  }
  fputs("int MultiplyDeBruijnBitPosition[64] = { ", stdout);
  for (i = 0; i < 63; ++i)
    printf("%d, ", MultiplyDeBruijnBitPosition[i]);
  printf("%d };\n", MultiplyDeBruijnBitPosition[i]);
}
#endif

extern int
bsCompare(constBitString a, BitOffset offsetA, BitOffset numBitsA,
          constBitString b, BitOffset offsetB, BitOffset numBitsB)
{
  BitOffset totalBitsLeftA = numBitsA, totalBitsLeftB = numBitsB;
  size_t elemStartA = offsetA/bitElemBits, elemStartB = offsetB/bitElemBits;
  unsigned bitTopA = offsetA%bitElemBits, bitTopB = offsetB%bitElemBits;
  const BitElem *pA = a + elemStartA, *pB = b + elemStartB;
  unsigned long accumA = 0, accumB = 0;
  unsigned bitsInAccumA, bitsInAccumB;
  assert(a && b);
  /* user requested zero length comparison, treat as equality */
  if (!numBitsA && !numBitsB)
  {
    return 0;
  }
  if (numBitsA > numBitsB)
    return -1 * bsCompare(b, offsetB, numBitsB, a, offsetA, numBitsA);
  if (numBitsB > numBitsA)
  {
    /* B is longer and thus compared with virtual zeros prepended to A */
    unsigned comparePreBits = numBitsB - numBitsA;
    do {
      bitsInAccumB = 0;
      if (bitTopB)
      {
        unsigned long mask; /*< all of the bits we want to get from *pB */
        unsigned bits2Read = MIN(bitElemBits - bitTopB, comparePreBits);
        unsigned unreadRightBits = (bitElemBits - bitTopB - bits2Read);
        mask = (~((~(unsigned long)0) << bits2Read)) << unreadRightBits;
        accumB = ((*pB++) & mask) >> unreadRightBits;
        bitsInAccumB += bits2Read;
        totalBitsLeftB -= bits2Read;
        comparePreBits -= bits2Read;
      }
      while (bitsInAccumB < CHAR_BIT * sizeof (accumB) && comparePreBits)
      {
        unsigned bits2Read,
          bitsFree = (CHAR_BIT * sizeof (accumA)) - bitsInAccumB;
        unsigned long mask;
        bits2Read = MIN3(bitsFree, bitElemBits, comparePreBits);
        mask = ~((~(unsigned long)0) << bits2Read);
        accumB = accumB << bits2Read
          | (((*pB) >> (bitElemBits - bits2Read)) & mask);
        bitsInAccumB += bits2Read;
        totalBitsLeftB -= bits2Read;
        comparePreBits -= bits2Read;
        /* all of *pB consumed? */
        if (bits2Read == bitElemBits)
          ++pB, bitTopB = 0;
        else
          bitTopB = bits2Read;
      }
    } while (accumB == 0 && comparePreBits);
    if (accumB > 0)
      return -1;
  }
  do {
    bitsInAccumB = bitsInAccumA = 0;
    /* get bits of first element if not aligned */
    if (bitTopA)
    {
      unsigned long mask; /*< all of the bits we want to get from *pA */
      unsigned bits2Read = MIN(bitElemBits - bitTopA, totalBitsLeftA);
      unsigned unreadRightBits = (bitElemBits - bitTopA - bits2Read);
      mask = (~((~(unsigned long)0) << bits2Read)) << unreadRightBits;
      accumA = ((*pA++) & mask) >> unreadRightBits;
      bitsInAccumA += bits2Read;
      totalBitsLeftA -= bits2Read;
    }
    else
      accumA = 0;
    /* get bits of first element if not aligned */
    if (bitTopB)
    {
      unsigned long mask; /*< all of the bits we want to get from *pB */
      unsigned bits2Read = MIN(bitElemBits - bitTopB, totalBitsLeftB);
      unsigned unreadRightBits = (bitElemBits - bitTopB - bits2Read);
      mask = (~((~(unsigned long)0) << bits2Read)) << unreadRightBits;
      accumB = ((*pB++) & mask) >> unreadRightBits;
      bitsInAccumB += bits2Read;
      totalBitsLeftB -= bits2Read;
    }
    else
      accumB = 0;
    while (bitsInAccumA < (CHAR_BIT * sizeof (accumA)) && totalBitsLeftA)
    {
      unsigned bits2Read,
               bitsFree = (CHAR_BIT * sizeof (accumA)) - bitsInAccumA;
      unsigned long mask;
      bits2Read = MIN3(bitsFree, bitElemBits, totalBitsLeftA);
      mask = (~((~(unsigned long)0) << bits2Read));
      accumA = accumA << bits2Read
        | (((*pA) >> (bitElemBits - bits2Read)) & mask);
      bitsInAccumA += bits2Read;
      totalBitsLeftA -= bits2Read;
      /* all of *pA consumed? */
      if (bits2Read == bitElemBits)
        ++pA, bitTopA = 0;
      else
        bitTopA = bits2Read;
    }
    while (bitsInAccumB < (CHAR_BIT * sizeof (accumA)) && totalBitsLeftB)
    {
      unsigned bits2Read,
        bitsFree = (CHAR_BIT * sizeof (accumA)) - bitsInAccumB;
      unsigned long mask;
      bits2Read = MIN3(bitsFree, bitElemBits, totalBitsLeftB);
      mask = (~((~(unsigned long)0) << bits2Read));
      accumB = accumB << bits2Read
        | (((*pB) >> (bitElemBits - bits2Read)) & mask);
      bitsInAccumB += bits2Read;
      totalBitsLeftB -= bits2Read;
      /* all of *pB consumed? */
      if (bits2Read == bitElemBits)
        ++pB, bitTopB = 0;
      else
        bitTopB = bits2Read;
    }
  } while (accumA == accumB && totalBitsLeftA);
  return accumA > accumB?1:(accumA < accumB?-1:0);
}

void
bsCopy(constBitString src, BitOffset offsetSrc,
       BitString dest, BitOffset offsetDest, BitOffset numBits)
{
  size_t elemStartSrc = offsetSrc/bitElemBits,
    elemStartDest = offsetDest/bitElemBits;
  unsigned bitTopSrc = offsetSrc%bitElemBits,
    bitTopDest = offsetDest%bitElemBits;
  BitOffset bitsLeft = numBits;
  const BitElem *p = src + elemStartSrc;
  BitElem *q = dest + elemStartDest;
  assert(src && dest);
  /* special optimization if equally aligned data will be copied */
  if (bitTopSrc == bitTopDest)
  {
    if (bitTopSrc)
    {
      BitElem mask = (~(BitElem)0) >> bitTopSrc;
      if (numBits < bitElemBits - bitTopSrc)
      {
        unsigned backShift = bitElemBits - numBits - bitTopSrc;
        mask &= ~(BitElem)0 << backShift;
        *q = (*q & ~mask) | (*p & mask);
        /* TODO: try wether  r = a ^ ((a ^ b) & mask) is faster, see above */
        return;
      }
      else
      {
        *q = (*q & ~mask) | (*p & mask);
        ++p, ++q;
        bitsLeft -= bitElemBits - bitTopSrc;
      }
    }
    if (bitsLeft)
    {
      size_t completeElems = bitsLeft/bitElemBits;
      memcpy(q, p, completeElems);
      p += completeElems, q += completeElems;
      bitsLeft %= bitElemBits;
    }
    if (bitsLeft)
    {
      BitElem mask = (~(BitElem)0) << (bitElemBits - bitsLeft);
      *q = (*q & ~mask) | (*p & mask);
    }
  }
  else
  {
    unsigned long accum = 0;
    unsigned bitsInAccum = 0;
    while (bitsLeft && (bitTopSrc || bitTopDest))
    {
      if (bitTopSrc)
      {
        unsigned long mask;
        unsigned bits2Read = MIN3(bitElemBits - bitTopSrc, bitsLeft,
                                  sizeof (accum) * CHAR_BIT - bitsInAccum);
        unsigned unreadRightBits = (bitElemBits - bitTopSrc - bits2Read);
        mask = (~((~(unsigned long)0) << bits2Read));
        accum = (accum << bits2Read) | (((*p) >> unreadRightBits) & mask);
        bitsLeft -= bits2Read;
        bitsInAccum += bits2Read;
        if ((bitTopSrc += bits2Read) == bitElemBits)
          bitTopSrc = 0, ++p;
      }
      if (bitTopDest)
      {
        unsigned bits2Write = MIN(bitsLeft + bitsInAccum,
                                  bitElemBits - bitTopDest);
        while (bitsLeft >= bitElemBits
               && sizeof (accum) * CHAR_BIT - bitsInAccum > bitElemBits)
        {
          accum = accum << bitElemBits | (*p++);
          bitsLeft -= bitElemBits;
          bitsInAccum += bitElemBits;
        }
        if (bits2Write > bitsInAccum)
        {
          /* very inconvinient: we have to take all the bits we can get
           * just to fill the first incomplete element */
          while (bitsInAccum < bits2Write)
          {
            unsigned long mask;
            unsigned bits2Read = MIN3(sizeof (accum) * CHAR_BIT - bitsInAccum,
                                      bitsLeft, bitElemBits);
            unsigned unreadRightBits = (bitElemBits - bits2Read);
            mask = (~((~(unsigned long)0) << bits2Read)) << unreadRightBits;
            accum = (accum << bits2Read) | ((*p) & mask) >> unreadRightBits;
            bitsLeft -= bits2Read;
            bitsInAccum += bits2Read;
            bitTopSrc = bits2Read;
          }
        }
        /* accum holds enough bits to fill incomplete region at dest start */
        {
          unsigned unwrittenRightBits = bitElemBits
            - bitTopDest - bits2Write;
          unsigned long mask =
            (~((~(unsigned long)0) << bits2Write)) << unwrittenRightBits;
          *q = (*q & ~mask)
            | (((accum >> (bitsInAccum -= bits2Write))
                << unwrittenRightBits) & mask);
          if ((bitTopDest += bits2Write) == bitElemBits)
            ++q, bitTopDest = 0;
        }
        while (bitsInAccum >= bitElemBits)
        {
          *q++ = accum >> (bitsInAccum -= bitElemBits);
        }
      }
    }
    /* all reads and writes are aligned on BitElems in this loop */
    do
    {
      /* fill accumulator */
      while (bitsLeft >= bitElemBits
             && sizeof (accum) * CHAR_BIT - bitsInAccum > bitElemBits)
      {
        accum = accum << bitElemBits | (*p++);
        bitsInAccum += bitElemBits;
        bitsLeft -= bitElemBits;
      }
      /* write out accum */
      while (bitsInAccum >= bitElemBits)
      {
        *q++ = accum >> (bitsInAccum -= bitElemBits);
      }
    }
    while (bitsLeft >= bitElemBits);
    /* write remaining (trailing) bits */
    while (bitsLeft || bitsInAccum)
    {
      while (bitsLeft && bitsInAccum < sizeof (accum) * CHAR_BIT)
      {
        unsigned bits2Read = MIN3(bitsLeft, bitElemBits - bitTopSrc,
                                  sizeof (accum) * CHAR_BIT - bitsInAccum);
        unsigned unreadRightBits = (bitElemBits - bitTopSrc - bits2Read);
        unsigned long mask =
          ~((~(unsigned long)0) << bits2Read);
        accum = (accum << bits2Read) | (((*p) >> unreadRightBits) & mask);
        bitsLeft -= bits2Read;
        bitsInAccum += bits2Read;
        if ((bitTopSrc += bits2Read) == bitElemBits)
          ++p, bitTopSrc = 0;
      }
      while (bitsInAccum)
      {
        unsigned bits2Write = MIN(bitsInAccum, bitElemBits - bitTopDest),
          unwrittenRightBits = bitElemBits - bits2Write - bitTopDest;
        unsigned long mask = (~(unsigned long)0);
        if (bits2Write != bitElemBits)
        {
          mask = (~(mask << bits2Write)) << unwrittenRightBits;
        }
        *q = (*q & ~mask) |
            (((accum >> (bitsInAccum -= bits2Write))
              << unwrittenRightBits) & mask);
        if ((bitTopDest += bits2Write) == bitElemBits)
          bitTopDest = 0, ++q;
      }
    }
  }
}

void
bsClear(BitString str, BitOffset offset, BitOffset numBits, int bitVal)
{
  unsigned bitsLeft = numBits,
    bitTop = offset%bitElemBits;
  size_t elemStart = offset/bitElemBits;
  BitElem *p = str + elemStart;
  unsigned long bitPatSource = 0UL;
  assert(str);
  if (bitVal)
    bitPatSource = ~0UL;
  if (bitTop)
  {
    unsigned long mask = ~0UL;
    if (bitElemBits < (sizeof (unsigned long)*CHAR_BIT))
    {
      mask <<= bitElemBits;
    }
    else
    {
      mask = 0;
    }
    mask = (~mask) >> bitTop;
    if (numBits < bitElemBits - bitTop)
    {
      unsigned backShift = bitElemBits - numBits - bitTop;
      mask &= ~0UL << backShift;
      *p = (*p & ~mask) | (bitPatSource & mask);
      return;
    }
    else
    {
      bitsLeft -= bitElemBits - bitTop;
      *p = (*p & ~mask) | (bitPatSource & mask);
      ++p;
    }
  }
  {
    size_t fullElems = bitsLeft / bitElemBits;
    memset(p, bitPatSource, sizeof (BitElem) * fullElems);
    p += fullElems;
    bitsLeft -= fullElems * bitElemBits;
  }
  if (bitsLeft)
  {
    unsigned long mask = ((~0UL)<<(bitElemBits - bitsLeft));
    if (bitElemBits < (sizeof (unsigned long)*CHAR_BIT))
      mask &= (~(~0UL<<bitElemBits));
    *p = (*p & ~mask) | (bitPatSource & mask);
  }
}

extern BitOffset
bs1BitsCount(constBitString str, BitOffset offset, BitOffset numBits)
{
  uint32_t accum = 0;
  BitOffset weight = 0, bitsLeft = numBits;
  unsigned bitTop = offset%bitElemBits, bitsInAccum = 0;
  size_t elemStart = offset/bitElemBits;
  const BitElem *p = str + elemStart;
  assert(str);
  if (bitTop)
  {
    uint32_t mask;
    unsigned bits2Read = MIN(bitElemBits - bitTop, bitsLeft);
    unsigned unreadRightBits = (bitElemBits - bitTop - bits2Read);
    mask = (~((~(uint32_t)0) << bits2Read)) << unreadRightBits;
    weight += bitCountUInt32(((*p++) & mask) >> unreadRightBits);
    bitsLeft -= bits2Read;
  }
  /* get bits from intervening elems */
  while (bitsLeft >= bitElemBits)
  {
    while (bitsLeft >= bitElemBits
           && sizeof (accum) * CHAR_BIT - bitElemBits >= bitsInAccum)
    {
      accum = accum << bitElemBits | (*p++);
      bitsLeft -= bitElemBits;
      bitsInAccum += bitElemBits;
    }
    weight += bitCountUInt32(accum);
    accum = 0; bitsInAccum = 0;
  }
  /* get bits from last elem */
  if (bitsLeft)
  {
    accum = ((*p) & ((~(uint32_t)0)<<(bitElemBits - bitsLeft)));
    weight += bitCountUInt32(accum);
  }
  return weight;
}

static inline void
bits2buf(char *buf, uint32_t v, unsigned numBits)
{
  unsigned i = numBits;
  uint32_t mask = 1;
  buf[i] = '\0';
  while (i)
  {
    --i;
    buf[i] = ((v & mask)?'1':'0');
    mask <<= 1;
  }
}

#define ACCUM2FP(accum, bitCount)               \
  bits2buf(buf, (accum), bitCount);             \
  if (fputs(buf, fp) == EOF)                    \
  {                                             \
    ioError = 1;                                \
    break;                                      \
  }

extern int
bsPrint(FILE *fp, constBitString str, BitOffset offset, BitOffset numBits)
{
  uint32_t accum = 0;
  unsigned bitsLeft = numBits, bitTop = offset%bitElemBits, bitsInAccum = 0;
  size_t elemStart = offset/bitElemBits;
  const BitElem *p = str + elemStart;
  char buf[sizeof(accum) * CHAR_BIT];
  int ioError = 0;
  assert(str);
  do {
    if (bitTop)
    {
      uint32_t mask;
      unsigned bits2Read = MIN(bitElemBits - bitTop, bitsLeft);
      unsigned unreadRightBits = (bitElemBits - bitTop - bits2Read);
      mask = (~((~(uint32_t)0) << bits2Read)) << unreadRightBits;
      ACCUM2FP(((*p++) & mask) >> unreadRightBits, bits2Read);
      bitsLeft -= bits2Read;
    }
    /* get bits from intervening elems */
    while (bitsLeft >= bitElemBits && !ioError)
    {
      while (bitsLeft >= bitElemBits
             && sizeof (accum) * CHAR_BIT - bitElemBits >= bitsInAccum)
      {
        accum = accum << bitElemBits | (*p++);
        bitsLeft -= bitElemBits;
        bitsInAccum += bitElemBits;
      }
      ACCUM2FP(accum, bitsInAccum);
      accum = 0; bitsInAccum = 0;
    }
    if (ioError)
      break;
    /* get bits from last elem */
    if (bitsLeft)
    {
      accum = ((*p) & ((~(uint32_t)0)<<(bitElemBits - bitsLeft)))
        >> (bitElemBits - bitsLeft);
      ACCUM2FP(accum, bitsLeft);
    }
  } while (0);
  return ioError?-1:0;
}
