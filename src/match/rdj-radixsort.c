/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include <string.h>
#include "core/intbits.h"
#include "core/ma.h"
#include "core/stack-inlined.h"
#include "match/rdj-radixsort.h"

typedef uint8_t gt_radixsort_kmercode_t;
/* note: if kmercode_t size is changed, then the RADIXSORT_REVCOMPL
 * macro must be changed to handle the new size too */

typedef uint16_t gt_radixsort_bucketnum_t;

typedef struct {
  unsigned long *suffixes;
  unsigned long width;
  unsigned long depth;
} GtRadixsortBucketInfo;

#define GT_RADIXSORT_INSERTION_SORT_MAX 31UL

#define GT_RADIXSORT_UINT8_REVCOMPL(CODE) \
  (GT_RADIXSORT_KMERCODE_MAX ^ \
   (((((CODE) & (3 << 6)) >> 6) | (((CODE) & 3) << 6)) | \
    ((((CODE) & (3 << 4)) >> 2) | (((CODE) & (3 << 2)) << 2))))

#define GT_RADIXSORT_REVCOMPL(CODE) GT_RADIXSORT_UINT8_REVCOMPL(CODE)

#define GT_RADIXSORT_KMERCODE_BITS    \
  (sizeof (gt_radixsort_kmercode_t) * CHAR_BIT)
#define GT_RADIXSORT_KMERSIZE         (GT_RADIXSORT_KMERCODE_BITS >> 1)
#define GT_RADIXSORT_KMERSIZELOG      2
#define GT_RADIXSORT_NOFKMERCODES     (1 << GT_RADIXSORT_KMERCODE_BITS)
#define GT_RADIXSORT_KMERCODE_MAX     \
  ((gt_radixsort_kmercode_t)(GT_RADIXSORT_NOFKMERCODES - 1))

#define GT_RADIXSORT_NOFBUCKETS \
  (size_t)(((1 << (GT_RADIXSORT_KMERSIZE << 1)) * GT_RADIXSORT_KMERSIZE) + 1)

#define GT_RADIXSORT_BUCKETNUM(CODE, OVERFLOW) \
  (((gt_radixsort_bucketnum_t)(CODE) << GT_RADIXSORT_KMERSIZELOG) + (OVERFLOW))

static inline gt_radixsort_kmercode_t gt_radixsort_code_at_position(
    const GtTwobitencoding *twobitencoding, unsigned long pos)
{
  unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  unsigned long unitindex = GT_DIVBYUNITSIN2BITENC(pos);

  if (unitoffset <=
      (unsigned int)(GT_UNITSIN2BITENC - GT_RADIXSORT_KMERSIZE))
  {
    return (gt_radixsort_kmercode_t)
      (twobitencoding[unitindex] >>
       GT_MULT2(GT_UNITSIN2BITENC - GT_RADIXSORT_KMERSIZE -
         unitoffset))
      & GT_RADIXSORT_KMERCODE_MAX;
  }
  else
  {
    unsigned int shiftleft =
      GT_MULT2(unitoffset + (unsigned int)GT_RADIXSORT_KMERSIZE -
          GT_UNITSIN2BITENC);
    return (gt_radixsort_kmercode_t)
      ((twobitencoding[unitindex] << shiftleft) |
       (twobitencoding[unitindex + 1] >>
           (GT_MULT2(GT_UNITSIN2BITENC) - shiftleft)))
      & GT_RADIXSORT_KMERCODE_MAX;
  }
}

static inline gt_radixsort_bucketnum_t gt_radixsort_get_code(
    GtTwobitencoding *twobitencoding, unsigned long suffixnum,
    unsigned long depth, unsigned long seqlen, unsigned long totallength)
{
  uint8_t overflow = 0;
  gt_radixsort_kmercode_t code;
  unsigned long remaining, pos = suffixnum + depth;
  if (suffixnum % seqlen + depth > seqlen - 2)
    return GT_RADIXSORT_NOFBUCKETS - 1;
  if (suffixnum <= totallength)
  {
    remaining = seqlen - 1 - pos % seqlen;
    code = gt_radixsort_code_at_position(twobitencoding, pos);
    if (remaining < (unsigned long)GT_RADIXSORT_KMERSIZE)
    {
      overflow = GT_RADIXSORT_KMERSIZE - remaining;
      code |= (1 << (overflow << 1)) - 1;
    }
  }
  else
  {
    pos = ((totallength + 1) << 1) - pos - 1;
    remaining = pos % seqlen;
    pos -= (remaining > (unsigned long)GT_RADIXSORT_KMERSIZE) ?
      (unsigned long)GT_RADIXSORT_KMERSIZE : remaining;
    code = GT_RADIXSORT_REVCOMPL(gt_radixsort_code_at_position(twobitencoding,
          pos));
    if (remaining < (unsigned long)GT_RADIXSORT_KMERSIZE)
    {
      overflow = GT_RADIXSORT_KMERSIZE - remaining;
      code = (code << (overflow << 1)) | ((1 << (overflow << 1)) - 1);
    }
  }
  return GT_RADIXSORT_BUCKETNUM(code, overflow);
}

GT_STACK_DECLARESTRUCT(GtRadixsortBucketInfo, 1024);

static inline void gt_radixsort_insertionsort(GtTwobitencoding *twobitencoding,
    unsigned long seqlen, unsigned long totallength,
    const GtRadixsortBucketInfo *bucket)
{
  unsigned long i, j;
  for (i = 1UL; i < bucket->width; i++)
  {
    const unsigned long u = bucket->suffixes[i];
    for (j = i; j > 0; j--)
    {
      const unsigned long v = bucket->suffixes[j - 1];
      unsigned long depth;
      gt_radixsort_bucketnum_t unk = 0, vnk = 0;
      int uvcmp = 0;
      for (depth = bucket->depth; uvcmp == 0 && (unk & 3) == 0
          && (vnk & 3) == 0 && unk != GT_RADIXSORT_NOFBUCKETS - 1 &&
          vnk != GT_RADIXSORT_NOFBUCKETS - 1; depth += GT_RADIXSORT_KMERSIZE)
      {
        unk = gt_radixsort_get_code(twobitencoding, u, depth,
          seqlen, totallength);
        vnk = gt_radixsort_get_code(twobitencoding, v, depth,
          seqlen, totallength);
        uvcmp = (int)vnk - (int)unk;
      }
      if (uvcmp <= 0)
        break;
      bucket->suffixes[j] = v;
    }
    bucket->suffixes[j] = u;
  }
}

void gt_radixsort_eqlen(GtTwobitencoding *twobitencoding,
    unsigned long *suffixes, unsigned long depth, unsigned long width,
    unsigned long seqlen, unsigned long totallength)
{
  unsigned long i;
  GtStackGtRadixsortBucketInfo stack;
  unsigned long bucketindex[GT_RADIXSORT_NOFBUCKETS];
  unsigned long bucketsize[GT_RADIXSORT_NOFBUCKETS];
  unsigned long *sorted;
  gt_radixsort_bucketnum_t *oracle;
  GtRadixsortBucketInfo bucket;
  GtRadixsortBucketInfo subbucket;

  oracle = gt_malloc(sizeof (*oracle) * width);
  sorted = gt_malloc(sizeof (*sorted) * width);

  GT_STACK_INIT(&stack, 1024UL);

  bucket.suffixes = suffixes;
  bucket.width = width;
  bucket.depth = depth;
  GT_STACK_PUSH(&stack, bucket);

  while (!GT_STACK_ISEMPTY(&stack))
  {
    bucket = GT_STACK_POP(&stack);
    memset(&bucketsize, 0, sizeof (bucketsize));

    /* Loop A */
    for (i = 0; i < bucket.width; ++i)
      oracle[i] = gt_radixsort_get_code(twobitencoding, bucket.suffixes[i],
          bucket.depth, seqlen, totallength);

    for (i = 0; i < bucket.width; ++i)
      ++bucketsize[oracle[i]];

    bucketindex[0] = 0;
    for (i = 1UL; i < (unsigned long)GT_RADIXSORT_NOFBUCKETS; ++i)
      bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

    /* Loop B */
    for (i = 0; i < bucket.width; ++i)
      sorted[bucketindex[oracle[i]]++] = bucket.suffixes[i];

    memcpy(bucket.suffixes, sorted, sizeof (unsigned long) * bucket.width);

    subbucket.suffixes = bucket.suffixes;
    subbucket.depth = bucket.depth + GT_RADIXSORT_KMERSIZE;
    if (subbucket.depth < seqlen)
    {
      for (i = 0; i < (unsigned long)GT_RADIXSORT_NOFBUCKETS; ++i)
      {
        subbucket.width = bucketsize[i];
        if (subbucket.width > 1UL)
        {
          if (subbucket.width <= GT_RADIXSORT_INSERTION_SORT_MAX)
            gt_radixsort_insertionsort(twobitencoding, seqlen, totallength,
                &subbucket);
          else
          {
            GT_STACK_PUSH(&stack, subbucket);
          }
        }
        subbucket.suffixes += subbucket.width;
        gt_assert(subbucket.suffixes <= bucket.suffixes + bucket.width);
      }
      gt_assert(bucket.suffixes + bucket.width == subbucket.suffixes);
    }
  }

  GT_STACK_DELETE(&stack);
  gt_free(sorted);
  gt_free(oracle);
}
