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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include <errno.h>
#include <string.h>
#include "core/alphabet.h"
#include "core/compactulongstore.h"
#include "core/encseq.h"
#include "core/fa.h"
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
#include "match/rdj-contfinder.h"

#ifdef GT_READJOINER_LARGE_READSET /* > 2^32 reads */
typedef GtCompactUlongstore GtReadjoinerContfinderSeqnumsType;
typedef uint64_t gt_contfinder_seqnum_t;
#define GT_READJOINER_CONTFINDER_ALLOC_SEQNUMS(ARR, NOFSEQS)\
  (ARR) = gt_GtCompactulongstore_new((NOFSEQS),\
      gt_determinebitspervalue(NOFSEQS))
#define GT_READJOINER_CONTFINDER_GET_SEQNUM(ARR, POS)\
  gt_GtCompactulongstore_get((ARR), (POS))
#define GT_READJOINER_CONTFINDER_SET_SEQNUM(ARR, POS, VALUE) \
  gt_GtCompactulongstore_update((ARR), (POS), (VALUE))
#define GT_READJOINER_CONTFINDER_FREE_SEQNUMS(ARR)\
  gt_GtCompactulongstore_delete(ARR)
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
  unsigned long nofcontained;
  unsigned long nofdiscarded;
  unsigned long discardedlength;
  /* eqlen: */
  size_t len;
  /* varlen: */
  size_t *seppos;
  size_t seppos_alloc, seppos_nextfree;
  size_t maxlen;
  /* temp buffers: */
  gt_contfinder_seqnum_t *sorted;
  gt_contfinder_kmercode_t *oracle;
  gt_contfinder_overflow_t *overflows;
  /* information needed for encseq output */
  const char *indexname;
  GtFilelengthvalues *filelengthtab;
  GtStrArray *filenametab;
  unsigned long characterdistribution[GT_CONTFINDER_ALPHASIZE];
  size_t totallength;
  bool contained_deleted;
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
  gt_contfinder_seqnum_t u, v;
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
        break;
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
}

/* the following assumes the seqnums are sorted in the bucket */
static inline void gt_contfinder_mark_as_contained(
    GtContfinder contfinder, GtContfinderBucketInfo bucket,
    bool except_lowest_seqnum)
{
  size_t i, from = 0, to = bucket.nofseqs;
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
    }
    else if (first == last) /* palindromic */
      to--;
  }
  for (i = from; i < to; i++)
  {
    gt_contfinder_seqnum_t corrected =
      GT_READJOINER_CONTFINDER_GET_SEQNUM(bucket.seqnums,
          i + bucket.seqnums_offset);
    if (corrected >= contfinder.nofseqs)
      corrected = (contfinder.nofseqs << 1) - corrected - 1;
    GT_SETIBIT(contfinder.contained, corrected);
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

static inline void gt_contfinder_sequential_decode(
    const GtContfinder contfinder, GtFile *outfp)
{
  size_t pos, nextsep;
  const char code2char[] = "acgt";
  gt_contfinder_seqnum_t seqnum;
  GtTwobitencoding code;
  const GtTwobitencoding *nextencoded = contfinder.twobitencoding;
  unsigned short charsincode;
  char *decoded;
  size_t nextdecoded;

  decoded = gt_malloc(sizeof (char) * (contfinder.maxlen + 3));
  decoded[0] = '>';
  decoded[1] = '\n';
  nextdecoded = (size_t)2;

  seqnum = 0;
  while (GT_ISIBITSET(contfinder.contained, seqnum))
    seqnum++;
  nextsep = (contfinder.len > 0) ? contfinder.len * (seqnum + 1) - 1 :
    contfinder.seppos[seqnum];
  pos = seqnum == 0 ? 0 : (contfinder.len > 0) ? seqnum * contfinder.len :
    contfinder.seppos[seqnum - 1] + 1;
  nextencoded = contfinder.twobitencoding +
    GT_DIVBYUNITSIN2BITENC(pos);
  code = *(nextencoded++);
  charsincode = (unsigned short)GT_UNITSIN2BITENC -
    (unsigned short)GT_MODBYUNITSIN2BITENC(pos);

  while (true)
  {
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (unsigned short)GT_UNITSIN2BITENC;
    }
    if (pos++ < nextsep)
      decoded[nextdecoded++] = code2char[code >> ((--charsincode) << 1) & 3];
    else
    {
      /* output sequence */
      charsincode--;
      decoded[nextdecoded++] = '\n';
      decoded[nextdecoded] = '\0';
      gt_file_xfputs(decoded, outfp);
      nextdecoded = (size_t)2;

      seqnum++;
      if (seqnum == contfinder.nofseqs)
        break;
      if (GT_ISIBITSET(contfinder.contained, seqnum))
      {
        /* jump to next non-contained sequence */
        while (seqnum < contfinder.nofseqs &&
            GT_ISIBITSET(contfinder.contained, seqnum))
          seqnum++;
        if (seqnum == contfinder.nofseqs)
          break;
        pos = contfinder.len > 0 ? seqnum * contfinder.len :
          contfinder.seppos[seqnum - 1] + 1;
        nextencoded = contfinder.twobitencoding +
          GT_DIVBYUNITSIN2BITENC(pos);
        code = *(nextencoded++);
        charsincode = (unsigned short)GT_UNITSIN2BITENC -
          (unsigned short)GT_MODBYUNITSIN2BITENC(pos);
      }
      nextsep = (contfinder.len > 0) ? contfinder.len * (seqnum + 1) - 1 :
        contfinder.seppos[seqnum];
    }
  }
  gt_free(decoded);
}

static inline unsigned short gt_contfinder_random_access_output(
    const GtContfinder *const contfinder, const gt_contfinder_seqnum_t seqnum,
    GtTwobitencoding *const outputbuffer,
    GtTwobitencoding *const outputoffsetptr)
{
  size_t firstpos, firstcodeidx, lastpos, lastcodeidx, seqlen;
  unsigned short retval;
  GtTwobitencoding  inputoffset, outputoffset = *outputoffsetptr;
  firstpos = seqnum == 0 ? 0 : contfinder->len > 0 ? contfinder->len * seqnum :
    contfinder->seppos[seqnum - 1] + 1;
  firstcodeidx = GT_DIVBYUNITSIN2BITENC(firstpos);
  lastpos = contfinder->len > 0 ? contfinder->len * (seqnum + 1) - 1 :
    contfinder->seppos[seqnum];
  lastcodeidx = GT_DIVBYUNITSIN2BITENC(lastpos);
  seqlen = contfinder->len > 0 ? contfinder->len : lastpos - firstpos + 1;

  inputoffset = (GtTwobitencoding)
    GT_MODBYUNITSIN2BITENC(firstpos);
  retval = (unsigned short)GT_DIVBYUNITSIN2BITENC(seqlen +
      outputoffset + GT_UNITSIN2BITENC - 1);

  if (inputoffset == outputoffset)
  {
    if (outputoffset == 0)
      outputbuffer[0] = contfinder->twobitencoding[firstcodeidx];
    else
    {
      GtTwobitencoding mask;
      mask = ((GtTwobitencoding)1
          << GT_MULT2(GT_UNITSIN2BITENC - outputoffset)) - 1;
      outputbuffer[0] =
        (contfinder->twobitencoding[firstcodeidx] & mask) |
        (outputbuffer[0] & ~mask);
    }
    if (lastcodeidx > firstcodeidx)
      memcpy(outputbuffer + 1, contfinder->twobitencoding + firstcodeidx + 1,
          sizeof (GtTwobitencoding) * (lastcodeidx - firstcodeidx));
    *outputoffsetptr = (GtTwobitencoding)GT_MODBYUNITSIN2BITENC(lastpos + 1);
  }
  else if (inputoffset > outputoffset)
  {
    const GtTwobitencoding netoffset = inputoffset - outputoffset;
    const GtTwobitencoding shiftright =
      GT_MULT2(GT_UNITSIN2BITENC - netoffset);
    const GtTwobitencoding shiftleft = GT_MULT2(netoffset);
    size_t i;
    GtTwobitencoding *nextinoutputbuffer = outputbuffer;
    if (outputoffset == 0)
      outputbuffer[0] = contfinder->twobitencoding[firstcodeidx] << shiftleft;
    else
    {
      GtTwobitencoding mask;
      mask = ((GtTwobitencoding)1 << GT_MULT2(GT_UNITSIN2BITENC -
            outputoffset)) - 1;
      outputbuffer[0] =
        ((contfinder->twobitencoding[firstcodeidx] << shiftleft) & mask) |
        (outputbuffer[0] & ~mask);
    }
    for (i = firstcodeidx + 1; i <= lastcodeidx; i++)
    {
      *(nextinoutputbuffer) |= (contfinder->twobitencoding[i] >> shiftright);
      *(++nextinoutputbuffer) = (contfinder->twobitencoding[i] << shiftleft);
    }
    *outputoffsetptr =
      GT_MODBYUNITSIN2BITENC(outputoffset + seqlen);
  }
  else
  {
    const GtTwobitencoding netoffset = outputoffset - inputoffset;
    const GtTwobitencoding shiftright = GT_MULT2(netoffset);
    const GtTwobitencoding shiftleft =
      GT_MULT2(GT_UNITSIN2BITENC - netoffset);
    size_t i;
    GtTwobitencoding *nextinoutputbuffer = outputbuffer + 1;
    GtTwobitencoding mask;
    mask = ((GtTwobitencoding)1 << GT_MULT2(GT_UNITSIN2BITENC -
          outputoffset)) - 1;
    outputbuffer[0] =
      ((contfinder->twobitencoding[firstcodeidx] >> shiftright) & mask) |
      (outputbuffer[0] & ~mask);
    outputbuffer[1] =
      contfinder->twobitencoding[firstcodeidx] << shiftleft;
    for (i = firstcodeidx + 1; i <= lastcodeidx; i++)
    {
      *(nextinoutputbuffer) |= (contfinder->twobitencoding[i] >> shiftright);
      *(++nextinoutputbuffer) = (contfinder->twobitencoding[i] << shiftleft);
    }
    *outputoffsetptr =
      GT_MODBYUNITSIN2BITENC(outputoffset + seqlen);
  }
  return retval;
}

static inline void gt_contfinder_random_access_decode(GtContfinder contfinder,
    gt_contfinder_seqnum_t seqnum, char *decoded)
{
  GtTwobitencoding code;
  size_t pos, seqlen;
  unsigned short charsincode;
  const char code2char[] = "acgt";
  char *nextdecoded = decoded;
  const GtTwobitencoding *nextencoded = contfinder.twobitencoding;

  seqlen = (contfinder.len > 0) ? contfinder.len :
      (seqnum > 0) ? contfinder.seppos[seqnum] - contfinder.seppos[seqnum - 1] :
      contfinder.seppos[0] + 1;
  *(nextdecoded++) = '>';
  *(nextdecoded++) = '\n';
  pos = seqnum == 0 ? 0 : (contfinder.len > 0) ? seqnum * contfinder.len :
    contfinder.seppos[seqnum - 1] + 1;
  nextencoded = contfinder.twobitencoding +
    GT_DIVBYUNITSIN2BITENC(pos);
  charsincode = (unsigned short)GT_UNITSIN2BITENC -
    (unsigned short)GT_MODBYUNITSIN2BITENC(pos);
  code = *(nextencoded++);
  for (pos = 0; pos < seqlen - 1; pos++)
  {
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (unsigned short)GT_UNITSIN2BITENC;
    }
    *(nextdecoded++) = code2char[code >> ((--charsincode) << 1) & 3];
  }
  *(nextdecoded++) = '\n';
  *(nextdecoded) = '\0';
}

static void gt_contfinder_subtract_chardistri_of_contained(
    GtContfinder *contfinder, gt_contfinder_seqnum_t i)
{
  size_t first, last;
  GtTwobitencoding *firstcode, *lastcode, *code, nextcode;
  unsigned short firstoffset, lastoffset, u;
  first = contfinder->len * i;
  gt_assert(contfinder->len >= (size_t)2);
  last = first + contfinder->len - 2;
  firstcode = contfinder->twobitencoding + GT_DIVBYUNITSIN2BITENC(first);
  lastcode = contfinder->twobitencoding + GT_DIVBYUNITSIN2BITENC(last);
  firstoffset = (unsigned short)GT_MODBYUNITSIN2BITENC(first);
  lastoffset = (unsigned short)GT_MODBYUNITSIN2BITENC(last);
  if (firstcode != lastcode)
  {
    nextcode = *firstcode;
    for (u = 0; u < (unsigned short)(GT_UNITSIN2BITENC - firstoffset); u++)
    {
      contfinder->characterdistribution[nextcode & 3]--;
      nextcode >>= 2;
    }
    for (code = firstcode + 1; code < lastcode; code++)
    {
      nextcode = *code;
      for (u = 0; u < (unsigned short)GT_UNITSIN2BITENC; u++)
      {
        contfinder->characterdistribution[nextcode & 3]--;
        nextcode >>= 2;
      }
    }
    nextcode = *lastcode;
    if (lastoffset > 0)
      nextcode >>= GT_MULT2(GT_UNITSIN2BITENC - lastoffset);
    for (u = 0; u < lastoffset; u++)
    {
      contfinder->characterdistribution[nextcode & 3]--;
      nextcode >>= 2;
    }
  }
  else
  {
    nextcode = *firstcode;
    if (lastoffset > 0)
      nextcode >>= GT_MULT2(GT_UNITSIN2BITENC - lastoffset);
    gt_assert(lastoffset > firstoffset);
    for (u = 0; u < lastoffset - firstoffset + 1; u++)
    {
      contfinder->characterdistribution[nextcode & 3]--;
      nextcode >>= 2;
    }
  }
}

void gt_contfinder_delete_contained(GtContfinder *contfinder)
{
  GtTwobitencoding outputoffset;
  gt_contfinder_seqnum_t i;
  contfinder->nofcontained = 0;
  outputoffset = 0;
  for (i = 0; i < contfinder->nofseqs; i++)
  {
    if (!GT_ISIBITSET(contfinder->contained, i))
    {
      (void)gt_contfinder_random_access_output(contfinder, i,
          contfinder->twobitencoding + GT_DIVBYUNITSIN2BITENC(contfinder->len
            * (i - contfinder->nofcontained)), &outputoffset);
    }
    else
    {
      gt_contfinder_subtract_chardistri_of_contained(contfinder, i);
      contfinder->nofcontained++;
    }
  }
  contfinder->totallength -= (contfinder->len * contfinder->nofcontained);
  contfinder->nofseqs -= contfinder->nofcontained;
  contfinder->contained_deleted = true;
}

static inline gt_contfinder_kmercode_t gt_contfinder_find_less_frequent_char(
    GtContfinder *contfinder)
{
  gt_contfinder_kmercode_t i, sepcode;
  unsigned long lowest_value, value;
  lowest_value = contfinder->characterdistribution[0];
  sepcode = 0;
  for (i = (gt_contfinder_kmercode_t)1; i < (gt_contfinder_kmercode_t)4; i++)
  {
    value = contfinder->characterdistribution[i];
    if (value < lowest_value)
    {
      lowest_value = value;
      sepcode = i;
    }
  }
  return sepcode;
}

static void gt_contfinder_set_separators_to_less_frequent_char(GtContfinder
    *contfinder)
{
  gt_contfinder_kmercode_t sepcode;
  gt_contfinder_seqnum_t seqnum;
  size_t pos, codenum, posincode, shift;
  GtTwobitencoding code, mask;
  sepcode = gt_contfinder_find_less_frequent_char(contfinder);
  if (sepcode != GT_CONTFINDER_SEPARATOR)
  {
    for (seqnum = (gt_contfinder_seqnum_t)1;
        seqnum < contfinder->nofseqs;
        ++seqnum)
    {
      pos = seqnum * contfinder->len - 1;
      codenum = GT_DIVBYUNITSIN2BITENC(pos);
      posincode = GT_MODBYUNITSIN2BITENC(pos);
      code = contfinder->twobitencoding[codenum];
      shift = GT_MULT2(GT_UNITSIN2BITENC - 1 - posincode);
      mask = ~((GtTwobitencoding)(3) << shift);
      gt_assert((code & ~mask) >> shift ==
          (GtTwobitencoding)GT_CONTFINDER_SEPARATOR);
      code = (code & mask) | ((GtTwobitencoding)sepcode << shift);
      contfinder->twobitencoding[codenum] = code;
    }
  }
  /* delete anything after the last char */
  pos = contfinder->nofseqs * contfinder->len - 2;
  codenum = GT_DIVBYUNITSIN2BITENC(pos);
  posincode = GT_MODBYUNITSIN2BITENC(pos);
  if (posincode < (size_t)(GT_UNITSIN2BITENC - 1))
  {
    shift = GT_MULT2(GT_UNITSIN2BITENC - 1 - posincode);
    contfinder->twobitencoding[codenum] =
      (contfinder->twobitencoding[codenum] >> shift) << shift;
  }
}

static int gt_contfinder_output_fasta(GtContfinder *contfinder, GtError *err)
{
  GtStr *fas_path;
  GtFile *fas_outfp;
  int had_err = 0;

  fas_path = gt_str_new_cstr(contfinder->indexname);
  gt_str_append_cstr(fas_path, GT_READJOINER_SUFFIX_PREFILTERED_FAS);
  gt_log_log("output prefiltered reads to %s", gt_str_get(fas_path));
  fas_outfp = gt_file_new(gt_str_get(fas_path), "w", err);
  if (fas_outfp == NULL)
    had_err = -1;
  else
  {
    gt_contfinder_sequential_decode(*contfinder, fas_outfp);
  }
  gt_file_delete(fas_outfp);
  gt_str_delete(fas_path);
  return had_err;
}

static int gt_contfinder_output_encseq(GtContfinder *contfinder, GtError *err)
{
  GtAlphabet *dna = gt_alphabet_new_dna();
  FILE *a_file = gt_fa_fopen_with_suffix(contfinder->indexname,
      GT_ALPHABETFILESUFFIX, "wb", err);
  int had_err = 0;
  gt_assert(contfinder->len > 0);
  gt_log_log("output prefiltered reads to encseq %s", contfinder->indexname);
  if (a_file == NULL)
    had_err = -1;
  gt_alphabet_output(dna, a_file);
  gt_fa_fclose(a_file);
  gt_alphabet_delete(dna);
  if (!had_err)
  {
    had_err = gt_encseq_write_twobitencoding_to_file(contfinder->indexname,
        (unsigned long)contfinder->totallength, (unsigned long)(contfinder->len
          - 1), contfinder->twobitencoding, (unsigned long)contfinder->nofseqs,
        gt_str_array_size(contfinder->filenametab), contfinder->filelengthtab,
        contfinder->filenametab, contfinder->characterdistribution, err);
  }
  return had_err;
}

int gt_contfinder_run(GtContfinder *contfinder, bool rev, GtFile *outfp,
    GtContfinderOutputFormat format, bool sorted, const char *cntlistfilename,
    const char *sepposfilename, bool output_encseq, GtError *err)
{
  gt_contfinder_seqnum_t i;
  GtReadjoinerContfinderSeqnumsType *seqnums;
  GtBitsequence *contained;
  GtContfinderBucketInfo all;
  int had_err = 0;

  if (contfinder->nofseqs == 0)
    return 0;

  GT_INITBITTAB(contained, contfinder->logicalnofseqs);

  contfinder->nofseqs = contfinder->logicalnofseqs;

  if (rev)
    contfinder->logicalnofseqs <<= 1;

  GT_READJOINER_CONTFINDER_ALLOC_SEQNUMS(seqnums, contfinder->logicalnofseqs);

  contfinder->contained = contained;
  all.nofseqs = contfinder->logicalnofseqs;
  all.seqnums = seqnums;
  all.seqnums_offset = 0;
  all.depth = 0;

  if (contfinder->len > 0)
    gt_contfinder_radixsort_eqlen(*contfinder, all);
  else
    gt_contfinder_radixsort(*contfinder, all);

  if (sorted)
  {
    switch (format)
    {
      case GT_CONTFINDER_SEQNUMS:
        for (i = 0; i < contfinder->logicalnofseqs; i++)
          if (!GT_ISIBITSET(contained,
                GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i)))
            gt_file_xprintf(outfp, "%lu\n",
                (unsigned long)GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i));
        break;
      case GT_CONTFINDER_FASTA:
        {
          char *decoded;
          decoded = gt_malloc(sizeof (char) * (contfinder->maxlen + 3));
          for (i = 0; i < contfinder->logicalnofseqs; i++)
          {
            if (GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i) <
                contfinder->nofseqs && !GT_ISIBITSET(contained,
                  GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i)))
            {
              gt_contfinder_random_access_decode(*contfinder,
                  GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i), decoded);
              gt_file_xfputs(decoded, outfp);
            }
          }
          gt_free(decoded);
        }
        break;
      case GT_CONTFINDER_2BIT:
        {
          GtTwobitencoding *outputbuffer, lastcode = 0,
                                         outputoffset;
          outputbuffer = gt_malloc(sizeof (GtTwobitencoding) *
              (GT_DIVBYUNITSIN2BITENC(contfinder->maxlen) + 3));
          outputoffset = 0;
          for (i = 0; i < contfinder->logicalnofseqs; i++)
          {
            if (GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i) <
                contfinder->nofseqs && !GT_ISIBITSET(contained,
                  GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i)))
            {
              unsigned short ncodes = gt_contfinder_random_access_output(
                  contfinder, GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i),
                  outputbuffer, &outputoffset);
              if (outputoffset > 0)
                lastcode = outputbuffer[--ncodes];
              gt_file_xwrite(outfp, outputbuffer,
                  sizeof (GtTwobitencoding) * ncodes);
              outputbuffer[0] = lastcode;
            }
          }
          if (outputoffset > 0)
          {
            gt_file_xwrite(outfp, &lastcode,
                sizeof (GtTwobitencoding));
          }
          gt_free(outputbuffer);
        }
        break;
      case GT_CONTFINDER_QUIET:
        break;
    }
  }
  else
  {
    switch (format)
    {
      case GT_CONTFINDER_SEQNUMS:
        for (i = 0; i < contfinder->nofseqs; i++)
          if (!GT_ISIBITSET(contained, i))
            gt_file_xprintf(outfp, "%lu\n", (unsigned long)i);
        break;
      case GT_CONTFINDER_FASTA:
        gt_contfinder_sequential_decode(*contfinder, outfp);
        break;
      case GT_CONTFINDER_2BIT:
        {
          GtTwobitencoding *outputbuffer, lastcode = 0,
                                         outputoffset;
          outputbuffer = gt_malloc(sizeof (GtTwobitencoding) *
              (GT_DIVBYUNITSIN2BITENC(contfinder->maxlen) + 3));
          outputoffset = 0;
          for (i = 0; i < contfinder->nofseqs; i++)
          {
            if (!GT_ISIBITSET(contained, i))
            {
              unsigned short ncodes = gt_contfinder_random_access_output(
                  contfinder, i, outputbuffer, &outputoffset);
              if (outputoffset > 0)
                lastcode = outputbuffer[--ncodes];
              gt_file_xwrite(outfp, outputbuffer,
                  sizeof (GtTwobitencoding) * ncodes);
              outputbuffer[0] = lastcode;
            }
          }
          if (outputoffset > 0)
          {
            gt_file_xwrite(outfp, &lastcode,
                sizeof (GtTwobitencoding));
          }
          gt_free(outputbuffer);
        }
        break;
      case GT_CONTFINDER_QUIET:
        break;
    }
  }

  if (!had_err && output_encseq)
  {
    if (contfinder->len > 0)
    {
      gt_log_log("read length: %lu", (unsigned long) contfinder->len);
      gt_contfinder_delete_contained(contfinder);
      gt_contfinder_set_separators_to_less_frequent_char(contfinder);
      had_err = gt_contfinder_output_encseq(contfinder, err);
    }
    else
    {
      /* currently actually writes fasta */
      had_err = gt_contfinder_output_fasta(contfinder, err);
    }
  }

  if (!had_err && cntlistfilename != NULL)
  {
    had_err = gt_cntlist_show(contained, (unsigned long)(contfinder->nofseqs +
          contfinder->nofcontained), cntlistfilename, true, err);
  }
  if (!had_err && sepposfilename != NULL && contfinder->seppos != NULL)
  {
    FILE *file;
    file = gt_fa_fopen(sepposfilename, "wb", err);
    if (file == NULL)
      had_err = -1;
    else
    {
      size_t pos;
      if (sorted)
      {
        pos = 0;
        for (i = 0; i < contfinder->logicalnofseqs; i++)
        {
          GtReadjoinerContfinderSeqnumsType seqnums_i =
            GT_READJOINER_CONTFINDER_GET_SEQNUM(seqnums, i);
          if (seqnums_i < contfinder->nofseqs)
          {
            pos += seqnums_i > 0 ? contfinder->seppos[seqnums_i] -
              contfinder->seppos[seqnums_i - 1] - 1 : contfinder->seppos[0];
            gt_xfwrite(&pos, sizeof (size_t), (size_t)1, file);
            pos++;
          }
        }
      }
      else
      {
        if (!GT_ISIBITSET(contained, 0))
        {
          gt_xfwrite(contfinder->seppos, sizeof (size_t), (size_t)1, file);
          pos = contfinder->seppos[0] + 1;
        }
        else
          pos = 0;
        for (i = (size_t)1; i < contfinder->nofseqs; i++)
        {
          if (!GT_ISIBITSET(contained, i))
          {
            pos += contfinder->seppos[i] - contfinder->seppos[i - 1] - 1;
            gt_xfwrite(&pos, sizeof (size_t), (size_t)1, file);
            pos++;
          }
        }
      }
      gt_fa_fclose(file);
    }
  }

  GT_READJOINER_CONTFINDER_FREE_SEQNUMS(seqnums);
  return had_err;
}

void gt_contfinder_delete(GtContfinder *contfinder)
{
  if (contfinder != NULL)
  {
    gt_free(contfinder->twobitencoding);
    gt_free(contfinder->seppos);
    gt_free(contfinder->filelengthtab);
    gt_free(contfinder->contained);
    gt_free(contfinder);
  }
}

static void gt_contfinder_append_seppos(GtContfinder *contfinder,
    size_t pos, unsigned long totallength)
{
  if (contfinder->seppos == NULL)
  {
    /* allocate seppos */
    gt_assert(pos > 0);
    gt_assert(totallength > (unsigned long)pos);
    contfinder->seppos_alloc = (size_t)(totallength / pos);
    gt_log_log("rough estimate of nofseqs = %lu",
        (unsigned long)(contfinder->seppos_alloc));
    contfinder->seppos = gt_malloc(sizeof (contfinder->seppos) *
        (contfinder->seppos_alloc));
    contfinder->seppos_nextfree = 0;
  }
  if (contfinder->seppos_nextfree == contfinder->seppos_alloc)
  {
    /* enlarge seppos */
    (contfinder->seppos_alloc) += GT_CONTFINDER_SEPPOS_INC;
    contfinder->seppos = gt_realloc(contfinder->seppos,
        sizeof (contfinder->seppos) * contfinder->seppos_alloc);
  }
  (contfinder->seppos)[(contfinder->seppos_nextfree)++] = pos;
}

static int gt_contfinder_encode_files(GtContfinder *contfinder,
    GtError *err)
{
  int had_err = 0;
  unsigned short codepos, codepos_seqstart = 0;
  GtTwobitencoding *twobitencoding_nextfree = NULL,
                   *twobitencoding_nextfree_seqstart = NULL;
  size_t twobitencoding_alloc = 0;
  GtTwobitencoding kmercode = 0, kmercode_seqstart = 0, nextcode;
  FILE *file;
  char line[GT_CONTFINDER_READBUFFER_SIZE], c, *fgetsretval;
  GtTwobitencoding char2code[UCHAR_MAX + 1];
  size_t j;
  unsigned long totallength = 0;
  size_t len, pos, maxlen;
  unsigned long nofseqs;
  const char *filename;
  unsigned long i;
  int a;
  bool valid_sequence = true;
  unsigned long characterdistribution[GT_CONTFINDER_ALPHASIZE];
  bool varlen = false;

  gt_error_check(err);

  for (i = 0; i < gt_str_array_size(contfinder->filenametab); i++)
  {
    uint64_t l = (uint64_t)gt_file_size(gt_str_array_get(
          contfinder->filenametab, (unsigned long)i));
    contfinder->filelengthtab[i].length = l;
    contfinder->filelengthtab[i].effectivelength = l;
    totallength += (unsigned long)l;
  }

  twobitencoding_alloc = (size_t)GT_DIVBYUNITSIN2BITENC(totallength +
      GT_UNITSIN2BITENC - 1);
  contfinder->twobitencoding =
    gt_malloc(sizeof (GtTwobitencoding) * twobitencoding_alloc);
  twobitencoding_nextfree = contfinder->twobitencoding;

  (void)memset(char2code, (int)UCHAR_MAX,
      sizeof (GtTwobitencoding) * (UCHAR_MAX + 1));
  char2code[(unsigned char)'A'] = (GtTwobitencoding)0;
  char2code[(unsigned char)'a'] = (GtTwobitencoding)0;
  char2code[(unsigned char)'C'] = (GtTwobitencoding)1;
  char2code[(unsigned char)'c'] = (GtTwobitencoding)1;
  char2code[(unsigned char)'G'] = (GtTwobitencoding)2;
  char2code[(unsigned char)'g'] = (GtTwobitencoding)2;
  char2code[(unsigned char)'T'] = (GtTwobitencoding)3;
  char2code[(unsigned char)'t'] = (GtTwobitencoding)3;

  len = 0;
  maxlen = 0;
  pos = 0;
  nofseqs = 0;

  for (a = 0; a < GT_CONTFINDER_ALPHASIZE; a++)
  {
    contfinder->characterdistribution[a] = 0;
    characterdistribution[a] = 0;
  }

  contfinder->len = 0;
  codepos = 0;
  for (i = 0; i < gt_str_array_size(contfinder->filenametab); i++)
  {
    filename = gt_str_array_get(contfinder->filenametab, i);
    file = gt_fa_fopen(filename, "r", err);
    if (file == NULL)
      had_err = -1;
    if (!had_err)
    {
      while ((fgetsretval = fgets(line, (int)GT_CONTFINDER_READBUFFER_SIZE,
              file)) == line)
      {
        if (line[0] == '>')
        {
          /* save encoding position in case next sequence must be skipped */
          twobitencoding_nextfree_seqstart = twobitencoding_nextfree;
          codepos_seqstart = codepos;
          kmercode_seqstart = kmercode;
          if (++nofseqs > 1UL)
          {
            kmercode = (kmercode << 2) | GT_CONTFINDER_SEPARATOR;
            if (++codepos == (unsigned short)GT_UNITSIN2BITENC)
            {
              *(twobitencoding_nextfree++) = kmercode;
              codepos = 0;
              kmercode = 0;
            }
            if (valid_sequence)
            {
              len++;
              if (varlen)
              {
                {
                  if (len > maxlen)
                    maxlen = len;
                  pos += len;
                  gt_contfinder_append_seppos(contfinder, pos - 1, totallength);
                }
              }
              else
              {
                if (nofseqs > 2UL)
                {
                  if (len != contfinder->len)
                  {
                    unsigned long seqnum;
                    gt_log_log("readset is varlen: sequences 0..%lu are "
                        "%lu bp long, sequence %lu is %lu bp long",
                        nofseqs - 2UL, (unsigned long)contfinder->len - 1,
                        nofseqs - 1UL, (unsigned long)len - 1);
                    varlen = true;
                    maxlen = (len > contfinder->len) ? len : contfinder->len;
                    for (seqnum = 0; seqnum < nofseqs - 2UL; seqnum++)
                    {
                      pos += contfinder->len;
                      gt_contfinder_append_seppos(contfinder, pos - 1,
                          totallength);
                    }
                    contfinder->len = 0;
                    pos += len;
                    gt_contfinder_append_seppos(contfinder, pos - 1,
                        totallength);
                  }
                }
                else
                {
                  contfinder->len = len;
                }
              }
            }
            len = 0;
          }
          /* handle the case in which a description
             is longer than the line array: */
          while (strlen(line) == GT_CONTFINDER_READBUFFER_SIZE - (size_t)1
              && line[GT_CONTFINDER_READBUFFER_SIZE - 2] != '\n'
              && fgetsretval != NULL)
            fgetsretval = fgets(line, (int)GT_CONTFINDER_READBUFFER_SIZE, file);
          /* update character distribution */
          for (a = 0; a < GT_CONTFINDER_ALPHASIZE; a++)
          {
            contfinder->characterdistribution[a] += characterdistribution[a];
            characterdistribution[a] = 0;
          }
          valid_sequence = true;
        }
        else if (valid_sequence)
        {
          j = 0;
          while (true)
          {
            c = line[j++];
            if (valid_sequence && (nextcode = char2code[(unsigned char)c])
                != GT_CONTFINDER_CODE_UNDEF)
            {
              characterdistribution[nextcode]++;
              kmercode = (kmercode << 2) | nextcode;
              if (++codepos == (unsigned short)GT_UNITSIN2BITENC)
              {
                *(twobitencoding_nextfree++) = kmercode;
                codepos = 0;
                kmercode = 0;
              }
              len++;
            }
            else
            {
              if (c == '\0')
                break;
              if (!isspace(c))
              {
                if (valid_sequence)
                {
                  valid_sequence = false;
                  contfinder->nofdiscarded++;

                  /* undo everything since last sequence start */
                  twobitencoding_nextfree = twobitencoding_nextfree_seqstart;
                  codepos = codepos_seqstart;
                  kmercode = kmercode_seqstart;
                  nofseqs--;
                  for (a = 0; a < GT_CONTFINDER_ALPHASIZE; a++)
                    characterdistribution[a] = 0;
                  contfinder->discardedlength += len;
                  len = contfinder->len;
                }
                contfinder->discardedlength++;
              }
            }
          }
        }
      }
      if (ferror(file) != 0)
      {
        gt_error_set(err, "Error by reading file %s: %s", filename,
            strerror(errno));
        had_err = -1;
      }
    }
    gt_fa_fclose(file);
    contfinder->filelengthtab[i].effectivelength = varlen ?
      contfinder->seppos[nofseqs - 1] : (len + 1) * nofseqs - 1;
  }
  /* add last sequence to char distri */
  for (a = 0; a < GT_CONTFINDER_ALPHASIZE; a++)
  {
    contfinder->characterdistribution[a] += characterdistribution[a];
    characterdistribution[a] = 0;
  }
  len++;
  if (varlen)
  {
    if (len > maxlen)
      maxlen = len;
    pos += len;
    gt_contfinder_append_seppos(contfinder, pos - 1, totallength);
  }
  else if (nofseqs == 1UL)
  {
    contfinder->len = len;
  }
  else if (nofseqs == 2UL && contfinder->len != len)
  {
    varlen = true;
    maxlen = (len > contfinder->len) ? len : contfinder->len;
    pos += contfinder->len;
    gt_contfinder_append_seppos(contfinder, pos - 1, totallength);
    contfinder->len = 0;
    pos += len;
    gt_contfinder_append_seppos(contfinder, pos - 1, totallength);
  }
  if (!had_err)
  {
    if (codepos > 0)
    {
      *(twobitencoding_nextfree++) = (kmercode <<
          ((GT_UNITSIN2BITENC - codepos) << 1));
    }
  }

#ifdef GT_READJOINER_LARGE_READSET
  gt_log_log("type for seqnums: GtCompactUlongstore");
#else
  if (nofseqs > UINT32_MAX)
  {
    gt_error_set(err,
        "too many reads: recompile "
        "using make CFLAGS+=\"-DGT_READJOINER_LARGE_READSET\"");
    return -1;
  }
  gt_log_log("type for seqnums: uint32_t");
#endif
  contfinder->logicalnofseqs = (gt_contfinder_seqnum_t)nofseqs;
  contfinder->nofseqs = contfinder->logicalnofseqs;
  if (contfinder->nofseqs > 0)
  {
    contfinder->twobitencoding = gt_realloc(contfinder->twobitencoding,
        sizeof (GtTwobitencoding) * (twobitencoding_nextfree -
          contfinder->twobitencoding));
  }
  else
  {
    gt_free(contfinder->twobitencoding);
    contfinder->twobitencoding = NULL;
  }
  if (varlen)
  {
    contfinder->seppos = gt_realloc(contfinder->seppos,
        sizeof (contfinder->seppos) * (contfinder->nofseqs));
    contfinder->len = 0;
    contfinder->maxlen = maxlen;
    contfinder->totallength = contfinder->seppos[contfinder->nofseqs - 1];
  }
  else
  {
    contfinder->maxlen = contfinder->len;
    contfinder->totallength = contfinder->len * contfinder->nofseqs - 1;
  }
  return had_err;
}

unsigned long gt_contfinder_totallength_without_sep(GtContfinder *contfinder)
{
  return (unsigned long)(contfinder->totallength - contfinder->nofseqs + 1);
}

unsigned long gt_contfinder_discarded_length(GtContfinder *contfinder)
{
  return contfinder->discardedlength;
}

unsigned long gt_contfinder_nofseqs(GtContfinder *contfinder)
{
  return (unsigned long)contfinder->nofseqs;
}

unsigned long gt_contfinder_nofdiscarded(GtContfinder *contfinder)
{
  return contfinder->nofdiscarded;
}

unsigned long gt_contfinder_read_length(GtContfinder *contfinder)
{
  return contfinder->len > 0 ? (unsigned long)(contfinder->len - 1) : 0;
}

unsigned long gt_contfinder_nofcontained(GtContfinder *contfinder)
{
  if (!contfinder->contained_deleted)
  {
    size_t i;
    GtBitsequence v;
    const size_t tabsize = GT_NUMOFINTSFORBITS(contfinder->nofseqs);
    if (contfinder->contained == NULL)
      return 0;
    contfinder->nofcontained = 0;
    for (i = 0; i < tabsize; i++)
    {
      v = contfinder->contained[i];
      for (/**/; v; contfinder->nofcontained++)
        v &= v - 1;
    }
  }
  return contfinder->nofcontained;
}

GtContfinder* gt_contfinder_new(GtStrArray *filenames, GtStr *indexname,
    bool output_encseq, GtError *err)
{
  GtContfinder *contfinder;
  int had_err = 0;
  contfinder = gt_malloc(sizeof (*contfinder));
  contfinder->twobitencoding = NULL;
  contfinder->seppos = NULL;
  contfinder->filenametab = filenames;
  contfinder->filelengthtab = gt_malloc(sizeof (*contfinder->filelengthtab) *
      gt_str_array_size(contfinder->filenametab));
  contfinder->indexname = gt_str_get(indexname);
  contfinder->nofcontained = 0;
  contfinder->nofdiscarded = 0;
  contfinder->discardedlength = 0;
  contfinder->contained = NULL;
  contfinder->contained_deleted = false;
  had_err = gt_contfinder_encode_files(contfinder, err);
  if (!had_err && contfinder->len > 0 &&
      output_encseq && contfinder->nofseqs > 0)
  {
    gt_contfinder_set_separators_to_less_frequent_char(contfinder);
    had_err = gt_contfinder_output_encseq(contfinder, err);
  }
  if (had_err)
  {
    gt_contfinder_delete(contfinder);
    return NULL;
  }
  return contfinder;
}

/* radixsort_str tester */

#include "match/radixsort_str.h"

void gt_contfinder_radixsort_str_eqlen_tester(GtContfinder *contfinder,
    bool mirrored, unsigned long offset, unsigned long depth,
    unsigned long maxdepth, bool print)
{
  unsigned long *suffixes, totallength, width, i;
  totallength = (unsigned long)contfinder->nofseqs * contfinder->len - 1;
  width = (mirrored ? ((totallength + 1) << 1) : (totallength + 1));
  suffixes = gt_malloc(sizeof (unsigned long) * width);
  for (i = 0; i < width; i++)
    suffixes[i] = i;
  gt_radixsort_str_eqlen(contfinder->twobitencoding, suffixes, offset, depth,
      maxdepth, width - offset, (unsigned long)contfinder->len, totallength);
  if (print)
    for (i = 0; i < width; i++)
      printf("%lu\n", suffixes[i]);
  gt_free(suffixes);
}
