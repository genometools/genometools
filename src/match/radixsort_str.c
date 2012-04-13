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
#include "match/radixsort_str.h"

typedef uint8_t gt_radixsort_str_kmercode_t;
/* note: if kmercode_t size is changed, then the RADIXSORT_STR_REVCOMPL
 * macro must be changed to handle the new size too */

typedef uint16_t gt_radixsort_str_bucketnum_t;

typedef struct {
  unsigned long *suffixes;
  unsigned long width;
  unsigned long depth;
} GtRadixsortStrBucketInfo;

#define GT_RADIXSORT_STR_INSERTION_SORT_MAX 31UL

#define GT_RADIXSORT_STR_UINT8_REVCOMPL(CODE) \
  (GT_RADIXSORT_STR_KMERCODE_MAX ^ \
   (((((CODE) & (3 << 6)) >> 6) | (((CODE) & 3) << 6)) | \
    ((((CODE) & (3 << 4)) >> 2) | (((CODE) & (3 << 2)) << 2))))

#define GT_RADIXSORT_STR_REVCOMPL(CODE) GT_RADIXSORT_STR_UINT8_REVCOMPL(CODE)

#define GT_RADIXSORT_STR_KMERCODE_BITS    \
  (sizeof (gt_radixsort_str_kmercode_t) * CHAR_BIT)
#define GT_RADIXSORT_STR_KMERSIZE         (GT_RADIXSORT_STR_KMERCODE_BITS >> 1)
#define GT_RADIXSORT_STR_KMERSIZELOG      2
#define GT_RADIXSORT_STR_NOFKMERCODES     (1 << GT_RADIXSORT_STR_KMERCODE_BITS)
#define GT_RADIXSORT_STR_KMERCODE_MAX     \
  ((gt_radixsort_str_kmercode_t)(GT_RADIXSORT_STR_NOFKMERCODES - 1))

#define GT_RADIXSORT_STR_NOFBUCKETS \
  (size_t)(((1 << (GT_RADIXSORT_STR_KMERSIZE << 1)) * \
        GT_RADIXSORT_STR_KMERSIZE) + 1)

#define GT_RADIXSORT_STR_LAST_BUCKET (GT_RADIXSORT_STR_NOFBUCKETS - 1)

#define GT_RADIXSORT_STR_HAS_OVERFLOW(CODE) \
  (((CODE) & 3) > 0 || (CODE == GT_RADIXSORT_STR_LAST_BUCKET))

#define GT_RADIXSORT_STR_BUCKETNUM(CODE, OVERFLOW) \
  (((gt_radixsort_str_bucketnum_t)(CODE) << GT_RADIXSORT_STR_KMERSIZELOG) + \
   (OVERFLOW))

static inline gt_radixsort_str_kmercode_t gt_radixsort_str_code_at_position(
    const GtTwobitencoding *twobitencoding, unsigned long pos)
{
  unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  unsigned long unitindex = GT_DIVBYUNITSIN2BITENC(pos);

  if (unitoffset <=
      (unsigned int)(GT_UNITSIN2BITENC - GT_RADIXSORT_STR_KMERSIZE))
  {
    return (gt_radixsort_str_kmercode_t)
      (twobitencoding[unitindex] >>
       GT_MULT2(GT_UNITSIN2BITENC - GT_RADIXSORT_STR_KMERSIZE -
         unitoffset))
      & GT_RADIXSORT_STR_KMERCODE_MAX;
  }
  else
  {
    unsigned int shiftleft =
      GT_MULT2(unitoffset + (unsigned int)GT_RADIXSORT_STR_KMERSIZE -
          GT_UNITSIN2BITENC);
    return (gt_radixsort_str_kmercode_t)
      ((twobitencoding[unitindex] << shiftleft) |
       (twobitencoding[unitindex + 1] >>
           (GT_MULT2(GT_UNITSIN2BITENC) - shiftleft)))
      & GT_RADIXSORT_STR_KMERCODE_MAX;
  }
}

static inline gt_radixsort_str_bucketnum_t gt_radixsort_str_get_code(
    const GtTwobitencoding *twobitencoding, unsigned long suffixnum,
    unsigned long depth, unsigned long seqlen, unsigned long totallength)
{
  uint8_t overflow = 0;
  gt_radixsort_str_kmercode_t code;
  unsigned long remaining, pos = suffixnum + depth;
  if (suffixnum % seqlen + depth > seqlen - 2)
    return GT_RADIXSORT_STR_LAST_BUCKET;
  if (suffixnum <= totallength)
  {
    remaining = seqlen - 1 - pos % seqlen;
    code = gt_radixsort_str_code_at_position(twobitencoding, pos);
    if (remaining < (unsigned long)GT_RADIXSORT_STR_KMERSIZE)
    {
      overflow = GT_RADIXSORT_STR_KMERSIZE - remaining;
      code |= (1 << (overflow << 1)) - 1;
    }
  }
  else
  {
    pos = ((totallength + 1) << 1) - pos - 1;
    remaining = pos % seqlen;
    pos -= (remaining > (unsigned long)GT_RADIXSORT_STR_KMERSIZE) ?
      (unsigned long)GT_RADIXSORT_STR_KMERSIZE : remaining;
    code = GT_RADIXSORT_STR_REVCOMPL(gt_radixsort_str_code_at_position(
          twobitencoding, pos));
    if (remaining < (unsigned long)GT_RADIXSORT_STR_KMERSIZE)
    {
      overflow = GT_RADIXSORT_STR_KMERSIZE - remaining;
      code = (code << (overflow << 1)) | ((1 << (overflow << 1)) - 1);
    }
  }
  return GT_RADIXSORT_STR_BUCKETNUM(code, overflow);
}

GT_STACK_DECLARESTRUCT(GtRadixsortStrBucketInfo, 1024);

static inline void gt_radixsort_str_insertionsort(
    const GtTwobitencoding *twobitencoding, unsigned long seqlen,
    unsigned long maxdepth, unsigned long totallength,
    const GtRadixsortStrBucketInfo *bucket)
{
  unsigned long i, j;
  for (i = 1UL; i < bucket->width; i++)
  {
    const unsigned long u = bucket->suffixes[i];
    for (j = i; j > 0; j--)
    {
      const unsigned long v = bucket->suffixes[j - 1];
      unsigned long depth;
      gt_radixsort_str_bucketnum_t unk = 0, vnk = 0;
      gt_assert(maxdepth == 0 || bucket->depth <= maxdepth);
      int uvcmp = 0;
      for (depth = bucket->depth;
          (maxdepth == 0 || depth <= maxdepth)
          && uvcmp == 0 && !GT_RADIXSORT_STR_HAS_OVERFLOW(unk) &&
          !GT_RADIXSORT_STR_HAS_OVERFLOW(vnk);
          depth += GT_RADIXSORT_STR_KMERSIZE)
      {
        unk = gt_radixsort_str_get_code(twobitencoding, u, depth,
          seqlen, totallength);
        vnk = gt_radixsort_str_get_code(twobitencoding, v, depth,
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

void gt_radixsort_str_eqlen(const GtTwobitencoding *twobitencoding,
    unsigned long *suffixes, unsigned long offset, unsigned long depth,
    unsigned long maxdepth, unsigned long width, unsigned long seqlen,
    unsigned long totallength)
{
  unsigned long i;
  GtStackGtRadixsortStrBucketInfo stack;
  unsigned long bucketindex[GT_RADIXSORT_STR_NOFBUCKETS];
  unsigned long bucketsize[GT_RADIXSORT_STR_NOFBUCKETS];
  unsigned long *sorted;
  gt_radixsort_str_bucketnum_t *oracle;
  GtRadixsortStrBucketInfo bucket;
  GtRadixsortStrBucketInfo subbucket;

  gt_assert(maxdepth == 0 || depth <= maxdepth);
  gt_assert(maxdepth == 0 || maxdepth <= seqlen);

  oracle = gt_malloc(sizeof (*oracle) * width);
  sorted = gt_malloc(sizeof (*sorted) * width);

  GT_STACK_INIT(&stack, 1024UL);

  bucket.suffixes = suffixes + offset;
  bucket.width = width;
  bucket.depth = depth;
  GT_STACK_PUSH(&stack, bucket);

  while (!GT_STACK_ISEMPTY(&stack))
  {
    bucket = GT_STACK_POP(&stack);
    memset(&bucketsize, 0, sizeof (bucketsize));

    /* Loop A */
    for (i = 0; i < bucket.width; ++i)
      oracle[i] = gt_radixsort_str_get_code(twobitencoding, bucket.suffixes[i],
          bucket.depth, seqlen, totallength);

    for (i = 0; i < bucket.width; ++i)
      ++bucketsize[oracle[i]];

    bucketindex[0] = 0;
    for (i = 1UL; i < (unsigned long)GT_RADIXSORT_STR_NOFBUCKETS; ++i)
      bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

    /* Loop B */
    for (i = 0; i < bucket.width; ++i)
      sorted[bucketindex[oracle[i]]++] = bucket.suffixes[i];

    memcpy(bucket.suffixes, sorted, sizeof (unsigned long) * bucket.width);

    subbucket.suffixes = bucket.suffixes;
    subbucket.depth = bucket.depth + GT_RADIXSORT_STR_KMERSIZE;
    if (subbucket.depth < seqlen &&
        (maxdepth == 0 || subbucket.depth <= maxdepth))
    {
      for (i = 0; i < (unsigned long)GT_RADIXSORT_STR_NOFBUCKETS; ++i)
      {
        subbucket.width = bucketsize[i];
        if (subbucket.width > 1UL)
        {
          if (subbucket.width <= GT_RADIXSORT_STR_INSERTION_SORT_MAX)
            gt_radixsort_str_insertionsort(twobitencoding, seqlen, maxdepth,
                totallength, &subbucket);
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
