/*
  Copyright (c) 2011-2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#include <errno.h>
#include <string.h>
#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include "core/fa.h"
#include "core/error_api.h"
#include "core/compact_ulong_store.h"
#include "core/encseq.h"
#include "core/filelengthvalues.h"
#include "core/fileutils.h"
#include "core/intbits.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/stack-inlined.h"
#include "core/str_array_api.h"
#include "core/xansi_api.h"
#include "match/rdj-cntlist.h"
#include "match/rdj-filesuf-def.h"
#include "match/reads2twobit.h"
#include "match/rdj-contfinder.h"

#ifdef GT_READJOINER_LARGE_READSET /* > 2^32 reads */
typedef GtCompactUlongStore GtReadjoinerContfinderSeqnumsType;
typedef uint64_t gt_contfinder_seqnum_t;
#define GT_READJOINER_CONTFINDER_ALLOC_SEQNUMS(ARR, NOFSEQS)\
  (ARR) = gt_compact_ulong_store_new((NOFSEQS),\
      gt_determinebitspervalue(NOFSEQS))
#define GT_READJOINER_CONTFINDER_GET_SEQNUM(ARR, POS)\
  gt_compact_ulong_store_get((ARR), (POS))
#define GT_READJOINER_CONTFINDER_SET_SEQNUM(ARR, POS, VALUE) \
  gt_compact_ulong_store_update((ARR), (POS), (VALUE))
#define GT_READJOINER_CONTFINDER_FREE_SEQNUMS(ARR)\
  gt_compact_ulong_store_delete(ARR)
#else
typedef uint32_t GtReadjoinerContfinderSeqnumsType;
typedef uint32_t gt_contfinder_seqnum_t;
#define GT_READJOINER_CONTFINDER_ALLOC_SEQNUMS(ARR, NOFSEQS)\
  (ARR) = gt_malloc(sizeof (GtReadjoinerContfinderSeqnumsType) * (NOFSEQS))
#define GT_READJOINER_CONTFINDER_GET_SEQNUM(ARR, POS)\
  ((ARR)[POS])
#define GT_READJOINER_CONTFINDER_SET_SEQNUM(ARR, POS, VALUE)\
  ((ARR)[POS]=(VALUE))
#define GT_READJOINER_CONTFINDER_FREE_SEQNUMS(ARR)\
  gt_free(ARR)
#endif

typedef uint8_t gt_contfinder_kmercode_t;
typedef uint8_t gt_contfinder_overflow_t;

typedef uint8_t gt_contfinder_copynum_t;
#define GT_CONTFINDER_COPYNUM_MAX UINT8_MAX
#define GT_CONTFINDER_COPYNUM_INC(VALUE, INC) \
  if ((VALUE) < GT_CONTFINDER_COPYNUM_MAX - INC + 1)\
  {\
    (VALUE) += INC;\
  }\

#define GT_CONTFINDER_CONTAINS(CONTFINDERPTR, CONTAINER, CONTAINED)\
  GT_SETIBIT((CONTFINDERPTR)->contained, (CONTAINED));\
  if ((CONTFINDERPTR)->copynum != NULL)\
  {\
    gt_assert(!GT_ISIBITSET((CONTFINDERPTR)->contained, CONTAINER));\
    GT_CONTFINDER_COPYNUM_INC((CONTFINDERPTR)->copynum[(CONTAINER)],\
        (CONTFINDERPTR)->copynum[(CONTAINED)]);\
    (CONTFINDERPTR)->copynum[(CONTAINED)] = 0;\
  }

typedef struct {
  gt_contfinder_kmercode_t code;
  gt_contfinder_overflow_t overflow;
} GtContfinderKmercodeWithOverflow;

#define GT_CONTFINDER_INSERTION_SORT_MAX (size_t)31

#define GT_CONTFINDER_KMERCODE_BITS    \
  (sizeof (gt_contfinder_kmercode_t) * CHAR_BIT)
#define GT_CONTFINDER_KMERSIZE         (GT_CONTFINDER_KMERCODE_BITS >> 1)
#define GT_CONTFINDER_NOFKMERCODES     (1 << GT_CONTFINDER_KMERCODE_BITS)
#define GT_CONTFINDER_KMERCODE_MAX     \
  ((gt_contfinder_kmercode_t)(GT_CONTFINDER_NOFKMERCODES - 1))

#define GT_CONTFINDER_NOFBUCKETS(LEVEL) \
  (1 << ((GT_CONTFINDER_KMERSIZE - (LEVEL)) << 1))

#define GT_CONTFINDER_CODE_UNDEF     (GtTwobitencoding)ULONG_MAX

#define GT_CONTFINDER_UINT8_REVCOMPL(CODE) \
  (GT_CONTFINDER_KMERCODE_MAX ^ \
   (((((CODE) & (3 << 6)) >> 6) | (((CODE) & 3) << 6)) | \
    ((((CODE) & (3 << 4)) >> 2) | (((CODE) & (3 << 2)) << 2))))

#define GT_CONTFINDER_UINT16_REVCOMPL(CODE) \
  ((gt_contfinder_kmercode_t)\
  (GT_CONTFINDER_KMERCODE_MAX ^ \
    (\
      ((((CODE) & (3UL << 14)) >> 14) | (((CODE) & (3UL      )) << 14)) | \
      ((((CODE) & (3UL << 12)) >> 10) | (((CODE) & (3UL <<  2)) << 10)) | \
      ((((CODE) & (3UL << 10)) >> 6)  | (((CODE) & (3UL <<  4)) << 6 )) | \
      ((((CODE) & (3UL << 8 )) >> 2 ) | (((CODE) & (3UL <<  6)) << 2))   \
    )\
  ))

#define GT_CONTFINDER_REVCOMPL(CODE) GT_CONTFINDER_UINT8_REVCOMPL(CODE)

#define GT_CONTFINDER_SEPARATOR (gt_contfinder_kmercode_t)3

#define GT_CONTFINDER_READBUFFER_SIZE ((size_t)256)

#define GT_CONTFINDER_SEPPOS_INC ((size_t)(1 << 14))

#define GT_CONTFINDER_ALPHASIZE 4

struct GtContfinder {
  GtTwobitencoding *twobitencoding;
  gt_contfinder_seqnum_t logicalnofseqs;
  gt_contfinder_seqnum_t nofseqs;
  GtBitsequence *contained;
  /* eqlen: */
  size_t len;
  /* varlen: */
  size_t *seppos;
  /* temp buffers: */
  gt_contfinder_seqnum_t *sorted;
  gt_contfinder_kmercode_t *oracle;
  gt_contfinder_overflow_t *overflows;
  /* copynum buffer */
  gt_contfinder_copynum_t *copynum;
  /* seqnums to sort */
  GtReadjoinerContfinderSeqnumsType *seqnums;
};

static inline gt_contfinder_kmercode_t gt_contfinder_code_at_position(
    const GtTwobitencoding *twobitencoding, size_t pos)
{
  unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  size_t unitindex = GT_DIVBYUNITSIN2BITENC(pos);

  if (unitoffset <=
      (unsigned int)(GT_UNITSIN2BITENC - GT_CONTFINDER_KMERSIZE))
  {
    return (gt_contfinder_kmercode_t)
      (twobitencoding[unitindex] >>
       GT_MULT2(GT_UNITSIN2BITENC - GT_CONTFINDER_KMERSIZE -
         unitoffset))
      & GT_CONTFINDER_KMERCODE_MAX;
  }
  else
  {
    unsigned int shiftleft = GT_MULT2(unitoffset +
        (unsigned int)GT_CONTFINDER_KMERSIZE -
        (unsigned int)GT_UNITSIN2BITENC);
    return (gt_contfinder_kmercode_t)
      ((twobitencoding[unitindex] << shiftleft) |
       (twobitencoding[unitindex + 1] >>
           (GT_MULT2(GT_UNITSIN2BITENC) - shiftleft)))
      & GT_CONTFINDER_KMERCODE_MAX;
  }
}

static inline GtContfinderKmercodeWithOverflow gt_contfinder_get_code(
    const GtTwobitencoding *twobitencoding, gt_contfinder_seqnum_t seqnum,
    const size_t depth, const gt_contfinder_seqnum_t firstrevcompl, size_t len,
    const size_t *seppos)
{
  GtContfinderKmercodeWithOverflow k;
  size_t remaining, seqstart;

  if (seqnum < firstrevcompl)
  {
    if (len != 0)
    {
      seqstart = seqnum * len;
    }
    else
    {
      seqstart = (seqnum > 0) ? seppos[seqnum - 1] + 1 : 0;
      len = seppos[seqnum] - seqstart + 1;
    }
    if (depth < len - 1)
    {
      k.code = gt_contfinder_code_at_position(twobitencoding, seqstart + depth);
      k.overflow = 0;
      remaining = len - 1 - depth;
      if (remaining < GT_CONTFINDER_KMERSIZE)
      {
        k.overflow = GT_CONTFINDER_KMERSIZE -
          (gt_contfinder_overflow_t)remaining;
        k.code &= (GT_CONTFINDER_KMERCODE_MAX - ((1 << (k.overflow << 1)) - 1));
      }
    }
    else
    {
      k.code = 0;
      k.overflow = GT_CONTFINDER_KMERSIZE;
    }
    return k;
  }
  else
  {
    size_t pos;
    if (len != 0)
    {
      pos = ((firstrevcompl << 1) - seqnum) * len - 1 - depth;
    }
    else
    {
      seqnum = (firstrevcompl << 1) - seqnum - 1;
      pos = seppos[seqnum] - depth;
      seqstart = (seqnum > 0) ? seppos[seqnum - 1] + 1 : 0;
      len = seppos[seqnum] - seqstart + 1;
    }
    if (depth < len - 1)
    {
      remaining = len - 1 - depth;
      pos -= (remaining > GT_CONTFINDER_KMERSIZE) ? GT_CONTFINDER_KMERSIZE :
        remaining;
      k.overflow = 0;
      k.code = GT_CONTFINDER_REVCOMPL(
          gt_contfinder_code_at_position(twobitencoding, pos));
      if (remaining < GT_CONTFINDER_KMERSIZE)
      {
        k.overflow = GT_CONTFINDER_KMERSIZE - remaining;
        k.code <<= (k.overflow << 1);
      }
    }
    else
    {
      k.code = 0;
      k.overflow = GT_CONTFINDER_KMERSIZE;
    }
    return k;
  }
}

typedef struct {
  GtReadjoinerContfinderSeqnumsType *seqnums;
  unsigned long seqnums_offset;
  gt_contfinder_seqnum_t nofseqs;
  size_t depth;
} GtContfinderBucketInfo;

GT_STACK_DECLARESTRUCT(GtContfinderBucketInfo, 1024);

static inline void gt_contfinder_insertion_sort(
    const GtContfinder contfinder, GtContfinderBucketInfo bucket)
{
  size_t i, j;
  gt_contfinder_seqnum_t u, v, container;
  int uvcmp;

  gt_assert(bucket.nofseqs > (gt_contfinder_seqnum_t)1);

  for (i = (size_t)1; i < bucket.nofseqs; i++)
  {
    gt_contfinder_seqnum_t ucorrected;
    size_t ulen = 0;
    u = GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums,
        i + bucket.seqnums_offset);
    ucorrected = u;
    if (ucorrected >= contfinder.nofseqs)
      ucorrected = (contfinder.nofseqs << 1) - 1 - ucorrected;
    if (contfinder.len == 0)
    {
      ulen = (ucorrected > 0) ? contfinder.seppos[ucorrected] -
        contfinder.seppos[ucorrected - 1] : contfinder.seppos[ucorrected] + 1;
    }
    for (j = i; j > 0; j--)
    {
      size_t len = contfinder.len;
      size_t pos;
      size_t vlen = 0;
      gt_contfinder_seqnum_t vcorrected;
      GtContfinderKmercodeWithOverflow unk = {0, 0}, vnk = {0, 0};
      v = GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums, j - 1 +
          bucket.seqnums_offset);
      vcorrected = v;
      if (vcorrected >= contfinder.nofseqs)
        vcorrected = (contfinder.nofseqs << 1) - 1 - vcorrected;
      if (len == 0)
      {
        vlen = (vcorrected > 0) ? contfinder.seppos[vcorrected] -
          contfinder.seppos[vcorrected - 1] : contfinder.seppos[vcorrected] + 1;
        len = MIN(ulen, vlen);
      }
      for (pos = bucket.depth, uvcmp = 0; uvcmp == 0 && pos < len;
          pos += GT_CONTFINDER_KMERSIZE)
      {
        unk = gt_contfinder_get_code(contfinder.twobitencoding, u,
            pos, contfinder.nofseqs, contfinder.len,
            contfinder.seppos);
        vnk = gt_contfinder_get_code(contfinder.twobitencoding, v,
            pos, contfinder.nofseqs, contfinder.len,
            contfinder.seppos);
        uvcmp = (int)vnk.code - (int)unk.code;
      }
      if (ulen > vlen)
      {
        unsigned short shift = (unsigned short)
          ((vnk.overflow - unk.overflow) << 1);
        uvcmp = (int)(vnk.code >> shift) - (int)(unk.code >> shift);
        if (uvcmp == 0)
        {
          /* v is contained in u */
          GT_SETIBIT(contfinder.contained, vcorrected);
          break;
        }
      }
      else if (ulen < vlen)
      {
        unsigned short shift = (unsigned short)
          ((unk.overflow - vnk.overflow) << 1);
        uvcmp = (int)(vnk.code >> shift) - (int)(unk.code >> shift);
        if (uvcmp == 0)
        {
          /* u is contained in v */
          GT_SETIBIT(contfinder.contained, ucorrected);
          break;
        }
      }
      if (uvcmp < 0)
      {
        break;
      }
      if (uvcmp == 0)
      {
        if (ucorrected > vcorrected)
        {
          GT_SETIBIT(contfinder.contained, ucorrected);
          break;
        }
      }
      GT_READJOINER_CONTFINDER_SET_SEQNUM(bucket.seqnums,
          j + bucket.seqnums_offset, v);
    }
    GT_READJOINER_CONTFINDER_SET_SEQNUM(bucket.seqnums,
        j + bucket.seqnums_offset, u);
  }

  container = GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums,
        bucket.seqnums_offset);
  if (container >= contfinder.nofseqs)
    container = (contfinder.nofseqs << 1) - 1 - container;
  for (i = (size_t)1; i < bucket.nofseqs; i++)
  {
    u = GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums,
        i + bucket.seqnums_offset);
    if (u >= contfinder.nofseqs)
      u = (contfinder.nofseqs << 1) - 1 - u;
    if (GT_ISIBITSET(contfinder.contained, u))
    {
      GT_CONTFINDER_CONTAINS(&contfinder, container, u);
    }
    else
    {
      container = u;
    }
  }
}

/* the following assumes the seqnums are sorted in the bucket */
static inline void gt_contfinder_mark_as_contained(
    GtContfinder contfinder, GtContfinderBucketInfo bucket,
    bool except_lowest_seqnum)
{
  size_t i, from = 0, to = bucket.nofseqs;
  gt_contfinder_seqnum_t container = 0;
  gt_assert(bucket.nofseqs > 0);
  if (except_lowest_seqnum)
  {
    gt_contfinder_seqnum_t first, last;
    if (bucket.nofseqs == (gt_contfinder_seqnum_t)1)
      return;
    first = GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums, 0 +
        bucket.seqnums_offset);
    last = GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums, bucket.nofseqs
        - 1 + bucket.seqnums_offset);
    from = (size_t)1;
    if (last >= contfinder.nofseqs)
      last = (contfinder.nofseqs << 1) - last - 1;
    if (first >= contfinder.nofseqs || last < first)
    {
      from--;
      to--;
      container = last;
    }
    else
    {
      container = first;
      if (first == last) /* palindromic */
        to--;
    }
  }
  for (i = from; i < to; i++)
  {
    gt_contfinder_seqnum_t corrected =
      GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums,
          i + bucket.seqnums_offset);
    if (corrected >= contfinder.nofseqs)
      corrected = (contfinder.nofseqs << 1) - corrected - 1;
    GT_CONTFINDER_CONTAINS(&contfinder, container, corrected);
  }
}

static void gt_contfinder_radixsort_eqlen(GtContfinder contfinder,
    GtContfinderBucketInfo all)
{
  size_t i;
  GtStackGtContfinderBucketInfo stack;
  gt_contfinder_seqnum_t bucketindex[GT_CONTFINDER_NOFBUCKETS(0)];
#ifndef S_SPLINT_S
  gt_contfinder_seqnum_t bucketsize[GT_CONTFINDER_NOFBUCKETS(0)] = {};
#endif
  GtContfinderBucketInfo subbucket;

  /* Loop A */
  for (i = 0; i < all.nofseqs; ++i)
    ++bucketsize[gt_contfinder_get_code(contfinder.twobitencoding, i, 0,
        contfinder.nofseqs, contfinder.len, NULL).code];

  bucketindex[0] = 0;
  for (i = (size_t)1; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
    bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

  /* Loop B */
  for (i = 0; i < all.nofseqs; ++i)
    GT_READJOINER_CONTFINDER_SET_SEQNUM(all.seqnums, (bucketindex[
        gt_contfinder_get_code(contfinder.twobitencoding, i, 0,
          contfinder.nofseqs, contfinder.len, NULL).code]++) +
        all.seqnums_offset, i);

  if (contfinder.len <= GT_CONTFINDER_KMERSIZE)
  {
    subbucket.depth = 0;
    subbucket.seqnums = all.seqnums;
    subbucket.seqnums_offset = all.seqnums_offset;
    for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
    {
      subbucket.nofseqs = bucketsize[i];
      if (subbucket.nofseqs > (size_t)1)
        gt_contfinder_mark_as_contained(contfinder, subbucket, true);
      subbucket.seqnums_offset += subbucket.nofseqs;
    }
    return;
  }

  /* alloc oracle and sorted */
  {
    gt_contfinder_seqnum_t maxbsize = bucketsize[0];
    for (i = (size_t)1; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
      if (maxbsize < bucketsize[i])
        maxbsize = bucketsize[i];
    contfinder.oracle = gt_malloc(sizeof (gt_contfinder_kmercode_t) * maxbsize);
    contfinder.sorted = gt_malloc(sizeof (gt_contfinder_seqnum_t) * maxbsize);
  }

  GT_STACK_INIT(&stack, 1024UL);

  subbucket.depth = GT_CONTFINDER_KMERSIZE;
  subbucket.seqnums = all.seqnums;
  subbucket.seqnums_offset = all.seqnums_offset;
  for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
  {
    subbucket.nofseqs = bucketsize[i];
    if (subbucket.nofseqs > (size_t)1)
    {
      if (subbucket.nofseqs <= GT_CONTFINDER_INSERTION_SORT_MAX)
        gt_contfinder_insertion_sort(contfinder, subbucket);
      else
      {
        GT_STACK_PUSH(&stack, subbucket);
      }
    }
    subbucket.seqnums_offset += subbucket.nofseqs;
  }

  while (!GT_STACK_ISEMPTY(&stack))
  {
    GtContfinderBucketInfo bucket = GT_STACK_POP(&stack);

    memset(bucketsize, 0, sizeof (bucketsize));

    /* Loop A */
    for (i = 0; i < bucket.nofseqs; ++i)
    {
      gt_contfinder_seqnum_t ith_seqnum =
        GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums, i +
            bucket.seqnums_offset);
      contfinder.oracle[i] = gt_contfinder_get_code(contfinder.twobitencoding,
          ith_seqnum, bucket.depth, contfinder.nofseqs, contfinder.len,
          NULL).code;
    }

    for (i = 0; i < bucket.nofseqs; ++i)
      ++bucketsize[contfinder.oracle[i]];

    bucketindex[0] = 0;
    for (i = (size_t)1; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
      bucketindex[i] = bucketindex[i - 1] + bucketsize[i - 1];

    /* Loop B */
    for (i = 0; i < bucket.nofseqs; ++i)
      contfinder.sorted[bucketindex[contfinder.oracle[i]]++] =
        GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums,
            i + bucket.seqnums_offset);
    for (i = 0; i < bucket.nofseqs; ++i)
      GT_READJOINER_CONTFINDER_SET_SEQNUM(bucket.seqnums,
          i + bucket.seqnums_offset,
          contfinder.sorted[i]);

    subbucket.depth = bucket.depth + GT_CONTFINDER_KMERSIZE;
    subbucket.seqnums = bucket.seqnums;
    subbucket.seqnums_offset = bucket.seqnums_offset;
    if (subbucket.depth < contfinder.len)
    {
      for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
      {
        subbucket.nofseqs = bucketsize[i];
        if (subbucket.nofseqs > (size_t)1)
        {
          if (subbucket.nofseqs <= GT_CONTFINDER_INSERTION_SORT_MAX)
            gt_contfinder_insertion_sort(contfinder, subbucket);
          else
          {
            GT_STACK_PUSH(&stack, subbucket);
          }
        }
        subbucket.seqnums_offset += subbucket.nofseqs;
      }
    }
    else
    {
      for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
      {
        subbucket.nofseqs = bucketsize[i];
        if (subbucket.nofseqs > (size_t)1)
          gt_contfinder_mark_as_contained(contfinder, subbucket, true);
        subbucket.seqnums_offset += subbucket.nofseqs;
      }
    }
  }

  GT_STACK_DELETE(&stack);
  gt_free(contfinder.sorted);
  gt_free(contfinder.oracle);
}

static inline void gt_contfinder_process_buckets(
    const GtContfinder contfinder, const GtContfinderBucketInfo parentbucket,
    gt_contfinder_seqnum_t **bucketsize,
    GtStackGtContfinderBucketInfo *stackptr, GtBitsequence **bucketcontained)
{
  GtContfinderBucketInfo subbucket;
  gt_contfinder_seqnum_t i;
  gt_contfinder_overflow_t overflow;

  /* process buckets with overflow = 0 */
  subbucket.seqnums = parentbucket.seqnums;
  subbucket.seqnums_offset = parentbucket.seqnums_offset;
  subbucket.depth = parentbucket.depth + GT_CONTFINDER_KMERSIZE;
  for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
  {
    subbucket.nofseqs = bucketsize[0][i];
    if (subbucket.nofseqs > 0)
    {
      if (subbucket.nofseqs > (size_t)1)
      {
        if (subbucket.nofseqs <= GT_CONTFINDER_INSERTION_SORT_MAX)
          gt_contfinder_insertion_sort(contfinder, subbucket);
        else
        {
          GT_STACK_PUSH(stackptr, subbucket);
        }
      }
      for (overflow = (gt_contfinder_overflow_t)1;
          overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
      {
        GT_SETIBIT(bucketcontained[overflow - 1], i >> (overflow << 1));
      }
    }
    subbucket.seqnums_offset += subbucket.nofseqs;
  }

  /* process buckets with overflow > 0 */
  for (overflow = (gt_contfinder_overflow_t)1;
      overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
  {
    for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(overflow); ++i)
    {
      subbucket.nofseqs = bucketsize[overflow][i];
      if (subbucket.nofseqs > 0)
        gt_contfinder_mark_as_contained(contfinder, subbucket,
            (bool)GT_ISIBITSET(bucketcontained[overflow - 1], i));
      subbucket.seqnums_offset += subbucket.nofseqs;
    }
  }
}

static void gt_contfinder_radixsort(GtContfinder contfinder,
    GtContfinderBucketInfo all)
{
  gt_contfinder_seqnum_t i;
  GtStackGtContfinderBucketInfo stack;
  GtContfinderKmercodeWithOverflow k;
  gt_contfinder_overflow_t overflow;
  GtBitsequence *bucketcontained[GT_CONTFINDER_KMERSIZE];
  gt_contfinder_seqnum_t *bucketindex[GT_CONTFINDER_KMERSIZE + 1];
  gt_contfinder_seqnum_t *bucketsize[GT_CONTFINDER_KMERSIZE + 1];
  size_t baseindex;

  for (overflow = 0; overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
  {
    bucketindex[overflow] =
      gt_malloc(sizeof (gt_contfinder_seqnum_t) *
          GT_CONTFINDER_NOFBUCKETS(overflow));
    bucketsize[overflow] =
      gt_calloc(GT_CONTFINDER_NOFBUCKETS(overflow),
          sizeof (gt_contfinder_seqnum_t));
  }
  for (overflow = (gt_contfinder_overflow_t)1;
      overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
  {
    GT_INITBITTAB(bucketcontained[overflow - 1],
        GT_CONTFINDER_NOFBUCKETS(overflow));
  }

  /* Loop A */
  for (i = 0; i < all.nofseqs; ++i)
  {
    k = gt_contfinder_get_code(contfinder.twobitencoding, i, 0,
        contfinder.nofseqs, contfinder.len, contfinder.seppos);
    ++bucketsize[k.overflow][k.code >> GT_MULT2(k.overflow)];
  }

  baseindex = 0;
  for (overflow = 0; overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
  {
    bucketindex[overflow][0] = baseindex;
    for (i = (size_t)1; i < GT_CONTFINDER_NOFBUCKETS(overflow); ++i)
      bucketindex[overflow][i] = bucketindex[overflow][i - 1] +
        bucketsize[overflow][i - 1];
    baseindex = bucketindex[overflow][GT_CONTFINDER_NOFBUCKETS(overflow) - 1]
      + bucketsize[overflow][GT_CONTFINDER_NOFBUCKETS(overflow) - 1];
  }

  /* Loop B */
  for (i = 0; i < all.nofseqs; ++i)
  {
    k = gt_contfinder_get_code(contfinder.twobitencoding, i, 0,
        contfinder.nofseqs, contfinder.len, contfinder.seppos);
    GT_READJOINER_CONTFINDER_SET_SEQNUM(all.seqnums,
        (bucketindex[k.overflow][k.code >> GT_MULT2(k.overflow)]++)
        + all.seqnums_offset, i);
  }

  /* alloc temp buffers */
  {
    gt_contfinder_seqnum_t maxbucketsize = bucketsize[0][0];
    for (i = (size_t)1; i < GT_CONTFINDER_NOFBUCKETS(0); ++i)
      if (maxbucketsize < bucketsize[0][i])
        maxbucketsize = bucketsize[0][i];
    contfinder.oracle = gt_malloc(sizeof (gt_contfinder_kmercode_t) *
        maxbucketsize);
    contfinder.sorted = gt_malloc(sizeof (gt_contfinder_seqnum_t) *
        maxbucketsize);
    contfinder.overflows = gt_malloc(sizeof (unsigned char) * maxbucketsize);
  }

  GT_STACK_INIT(&stack, 1024UL);

  gt_contfinder_process_buckets(contfinder, all, bucketsize, &stack,
      bucketcontained);

  while (!GT_STACK_ISEMPTY(&stack))
  {
    GtContfinderBucketInfo bucket = GT_STACK_POP(&stack);

    for (overflow = 0; overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
    {
      for (i = 0; i < GT_CONTFINDER_NOFBUCKETS(overflow); ++i)
      {
        bucketsize[overflow][i] = 0;
      }
    }

    /* Loop A */
    for (i = 0; i < bucket.nofseqs; ++i)
    {
      k = gt_contfinder_get_code(
          contfinder.twobitencoding, GT_READJOINER_CONTFINDER_GET_SEQNUM(
            bucket.seqnums, i + bucket.seqnums_offset), bucket.depth,
          contfinder.nofseqs, contfinder.len, contfinder.seppos);
      contfinder.oracle[i] = k.code >> GT_MULT2(k.overflow);
      contfinder.overflows[i] = k.overflow;
    }

    for (i = 0; i < bucket.nofseqs; ++i)
      ++(bucketsize[contfinder.overflows[i]][contfinder.oracle[i]]);

    baseindex = 0;
    for (overflow = 0; overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
    {
      bucketindex[overflow][0] = baseindex;
      for (i = (size_t)1; i < GT_CONTFINDER_NOFBUCKETS(overflow); ++i)
        bucketindex[overflow][i] = bucketindex[overflow][i - 1] +
          bucketsize[overflow][i - 1];
      baseindex = bucketindex[overflow][GT_CONTFINDER_NOFBUCKETS(overflow) - 1];
    }

    /* Loop B */
    for (i = 0; i < bucket.nofseqs; ++i)
      contfinder.sorted[(bucketindex[contfinder.overflows[i]]
        [contfinder.oracle[i]])++] = GT_READJOINER_CONTFINDER_GET_SEQNUM(
            bucket.seqnums, i + bucket.seqnums_offset);
    for (i = 0; i < bucket.nofseqs; ++i)
      GT_READJOINER_CONTFINDER_SET_SEQNUM(bucket.seqnums,
          i + bucket.seqnums_offset, contfinder.sorted[i]);

    gt_contfinder_process_buckets(contfinder, bucket, bucketsize, &stack,
        bucketcontained);
  }

  GT_STACK_DELETE(&stack);
  gt_free(contfinder.sorted);
  gt_free(contfinder.oracle);
  gt_free(contfinder.overflows);
  for (overflow = 0; overflow <= GT_CONTFINDER_KMERSIZE; overflow++)
  {
    gt_free(bucketindex[overflow]);
    gt_free(bucketsize[overflow]);
  }
  for (overflow = 0; overflow < GT_CONTFINDER_KMERSIZE; overflow++)
  {
    gt_free(bucketcontained[overflow]);
  }
}

int gt_contfinder_write_seqnums(GtContfinder *contfinder, bool sorted,
    GtFile *outfp, GT_UNUSED GtError *err)
{
  int had_err = 0;
  gt_contfinder_seqnum_t i;
  if (sorted)
  {
    for (i = 0; i < contfinder->logicalnofseqs; i++)
      if (!GT_ISIBITSET(contfinder->contained,
            GT_READJOINER_CONTFINDER_GET_SEQNUM(contfinder->seqnums, i)))
        gt_file_xprintf(outfp, "%lu\n", (unsigned long)
            GT_READJOINER_CONTFINDER_GET_SEQNUM(contfinder->seqnums, i));
  }
  else
  {
    for (i = 0; i < contfinder->nofseqs; i++)
      if (!GT_ISIBITSET(contfinder->contained, i))
        gt_file_xprintf(outfp, "%lu\n", (unsigned long)i);
  }
  return had_err;
}

int gt_contfinder_write_sorted_seqnums(GtContfinder *contfinder, char* path,
    GtError *err)
{
  int had_err = 0;
  GtFile *file;
  file = gt_file_new(path, "wb", err);
  if (file == NULL)
    had_err = -1;
  if (!had_err)
  {
    had_err = gt_contfinder_write_seqnums(contfinder, true, file, err);
  }
  gt_file_delete(file);
  return had_err;
}

int gt_contfinder_write_cntlist(GtContfinder *contfinder, char* path,
    GtError *err)
{
  int had_err;
  had_err = gt_cntlist_show(contfinder->contained,
      (unsigned long)contfinder->nofseqs, path, true, err);
  return had_err;
}

int gt_contfinder_write_copynum(GtContfinder *contfinder, char* path,
    GtError *err)
{
  int had_err = 0;
  gt_contfinder_seqnum_t i;
  FILE *file;
  file = gt_fa_fopen(path, "wb", err);
  if (file == NULL)
    had_err = -1;
  else
  {
#ifndef NDEBUG
    gt_contfinder_seqnum_t n_noncontained = 0;
    unsigned long cnsum = 0;
    bool had_overflow = false;
#endif
    gt_contfinder_copynum_t cn;
    for (i = 0; i < contfinder->nofseqs; i++)
    {
      cn = contfinder->copynum[i];
      if (cn > 0)
      {
        gt_xfwrite(contfinder->copynum + i, sizeof (*contfinder->copynum),
            (size_t)1, file);
#ifndef NDEBUG
        n_noncontained++;
        if (cn == GT_CONTFINDER_COPYNUM_MAX)
        {
          had_overflow = true;
        }
        cnsum += cn;
#endif
      }
    }
    gt_assert(n_noncontained == contfinder->nofseqs);
    gt_assert(had_overflow || (cnsum == (unsigned long)contfinder->nofseqs));
    gt_fa_fclose(file);
  }
  return had_err;
}

GtBitsequence *gt_contfinder_contained(GtContfinder *contfinder)
{
  gt_assert(contfinder);
  return contfinder->contained;
}

static void gt_contfinder_init_copynum(GtContfinder *contfinder)
{
  gt_contfinder_seqnum_t idx;
  /* the mark_as_contained function is currently not able to set
   * correctly the copynum if sequences have variable length,
   * thus: */
  gt_assert(contfinder->len > 0);
  contfinder->copynum = gt_malloc(sizeof (*contfinder->copynum) *
      contfinder->nofseqs);
  for (idx = 0; idx < contfinder->nofseqs; idx++)
    contfinder->copynum[idx] = (gt_contfinder_seqnum_t)1;
}

void gt_contfinder_run(GtContfinder *contfinder, bool mirrored,
    bool calculate_copynum)
{
  GtContfinderBucketInfo all;

  if (contfinder->nofseqs == 0)
    return;

  GT_INITBITTAB(contfinder->contained, contfinder->logicalnofseqs);

  contfinder->nofseqs = contfinder->logicalnofseqs;

  if (mirrored)
    contfinder->logicalnofseqs <<= 1;

  GT_READJOINER_CONTFINDER_ALLOC_SEQNUMS(contfinder->seqnums,
      contfinder->logicalnofseqs);

  if (calculate_copynum)
    gt_contfinder_init_copynum(contfinder);

  all.nofseqs = contfinder->logicalnofseqs;
  all.seqnums = contfinder->seqnums;
  all.seqnums_offset = 0;
  all.depth = 0;

  if (contfinder->len > 0)
    gt_contfinder_radixsort_eqlen(*contfinder, all);
  else
    gt_contfinder_radixsort(*contfinder, all);
}

void gt_contfinder_delete(GtContfinder *contfinder)
{
  if (contfinder != NULL)
  {
    GT_READJOINER_CONTFINDER_FREE_SEQNUMS(contfinder->seqnums);
    gt_free(contfinder->contained);
    gt_free(contfinder->copynum);
    gt_free(contfinder);
  }
}

unsigned long gt_contfinder_nofcontained(GtContfinder *contfinder)
{
  size_t i;
  GtBitsequence v;
  unsigned long nofcontained;
  const size_t tabsize = GT_NUMOFINTSFORBITS(contfinder->nofseqs);
  if (contfinder->contained == NULL)
    return 0;
  nofcontained = 0;
  for (i = 0; i < tabsize; i++)
  {
    v = contfinder->contained[i];
    for (/**/; v; nofcontained++)
      v &= v - 1;
  }
  return nofcontained;
}

GtContfinder* gt_contfinder_new(GtReads2Twobit *r2t)
{
  GtContfinder *contfinder;
  contfinder = gt_malloc(sizeof (*contfinder));
  contfinder->contained = NULL;
  contfinder->copynum = NULL;
  contfinder->seqnums = NULL;
  contfinder->twobitencoding = gt_reads2twobit_export_twobitencoding(r2t);
  contfinder->seppos = (size_t*)gt_reads2twobit_export_seppos(r2t);
  contfinder->nofseqs = (gt_contfinder_seqnum_t)gt_reads2twobit_nofseqs(r2t);
  contfinder->logicalnofseqs = contfinder->nofseqs;
  contfinder->len = (size_t)gt_reads2twobit_seqlen_eqlen(r2t);
  return contfinder;
}

/* radixsort_str tester */

#include "match/radixsort_str.h"

void gt_contfinder_radixsort_str_eqlen_tester(GtContfinder *contfinder,
    bool mirrored, unsigned long depth,
    unsigned long maxdepth, bool print)
{
  unsigned long *suffixes, totallength, width, i;
  GtRadixsortstringinfo *rsi;

  totallength = (unsigned long) contfinder->nofseqs * contfinder->len - 1;
  width = mirrored ? GT_MULT2(totallength + 1) : (totallength + 1);
  rsi = gt_radixsort_str_new(contfinder->twobitencoding,
                             totallength,
                             (unsigned long) contfinder->len,
                             GT_MULT2(width));
  suffixes = gt_malloc(sizeof (unsigned long) * width);
  for (i = 0; i < width; i++)
  {
    suffixes[i] = i;
  }
  gt_radixsort_str_eqlen(rsi, suffixes, NULL, 0, depth, maxdepth, width);
  if (print)
  {
    for (i = 0; i < width; i++)
    {
      printf("%lu\n", suffixes[i]);
    }
  }
  gt_free(suffixes);
  gt_radixsort_str_delete(rsi);
}
