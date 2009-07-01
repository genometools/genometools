/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/qsort_r.h"
#include "core/array2dim_api.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "sfx-copysort.h"

typedef struct
{
  bool sorted;
  unsigned long bucketend;
} Bucketinfo;

struct Bucketspec2
{
  Seqpos partwidth;
  unsigned int numofchars, *order;
  Bucketinfo *superbuckettab, **subbuckettab;
};

static unsigned long superbucketsize(const Bucketspec2 *bucketspec2,
                                     unsigned int bucketnum)
{
  if (bucketnum == 0)
  {
    return bucketspec2->superbuckettab[0].bucketend;
  }
  return bucketspec2->superbuckettab[bucketnum].bucketend -
         bucketspec2->superbuckettab[bucketnum-1].bucketend;
}

static int comparesuperbucketsizes(const void *a,const void *b,void *data)
{
  const Bucketspec2 *bucketspec2 = (const Bucketspec2 *) data;
  unsigned long size1 = superbucketsize(bucketspec2,
                                        *(const unsigned int *) a);
  unsigned long size2 = superbucketsize(bucketspec2,
                                        *(const unsigned int *) b);
  if (size1 < size2)
  {
    return -1;
  }
  if (size1 > size2)
  {
    return 1;
  }
  return 0;
}

static unsigned long getstartidx(const Bucketspec2 *bucketspec2,
                                 unsigned int first,
                                 unsigned int second)
{
  gt_assert(first < bucketspec2->numofchars);
  gt_assert(second <= bucketspec2->numofchars);
  if (second > 0)
  {
    return bucketspec2->subbuckettab[first][second-1].bucketend;
  }
  if (first > 0)
  {
    return bucketspec2->superbuckettab[first-1].bucketend;
  }
  return 0;
}

static unsigned long getendidx(const Bucketspec2 *bucketspec2,
                               unsigned int first,
                               unsigned int second)
{
  gt_assert(first < bucketspec2->numofchars);
  gt_assert(second <= bucketspec2->numofchars);
  if (second < bucketspec2->numofchars)
  {
    return bucketspec2->subbuckettab[first][second].bucketend;
  }
  return bucketspec2->superbuckettab[first].bucketend;
}

static void resetsorted(Bucketspec2 *bucketspec2)
{
  unsigned int idx, idx2;

  for (idx = 0; idx<bucketspec2->numofchars; idx++)
  {
    bucketspec2->superbuckettab[idx].sorted = false;
    for (idx2 = 0; idx2<bucketspec2->numofchars; idx2++)
    {
      unsigned long startidx = getstartidx(bucketspec2,idx,idx2),
                    endidx = getendidx(bucketspec2,idx,idx2);
      bucketspec2->subbuckettab[idx][idx2].sorted
        = (startidx < endidx) ? false : true;
    }
  }
}

static Codetype expandtwocharcode(Codetype twocharcode,
                                  unsigned int prefixlength,
                                  const Codetype *basepower,
                                  const Bcktab *bcktab)
{
  gt_assert(twocharcode < (Codetype) 16);
  return twocharcode * basepower[prefixlength-1] 
                     + bcktab_filltable(bcktab,2);
}

Bucketspec2 *bucketspec2_new(const Bcktab *bcktab,
                             Seqpos partwidth,
                             unsigned int numofchars)
{
  Codetype code, maxcode;
  Bucketspecification bucketspec;
  Bucketspec2 *bucketspec2;
  unsigned int idx, rightchar = 0, currentchar = 0, prefixlength;
  unsigned long accubucketsize = 0;

  gt_assert(numofchars > 0);
  bucketspec2 = gt_malloc(sizeof(*bucketspec2));
  bucketspec2->partwidth = partwidth;
  bucketspec2->numofchars = numofchars;
  bucketspec2->order = gt_malloc(sizeof(*bucketspec2->order) * numofchars);
  bucketspec2->superbuckettab
    = gt_malloc(sizeof(*bucketspec2->superbuckettab) * numofchars);
  gt_array2dim_malloc(bucketspec2->subbuckettab,(unsigned long) numofchars,
                      (unsigned long) numofchars);
  maxcode = bcktab_numofallcodes(bcktab) - 1;
  prefixlength = bcktab_prefixlength(bcktab);
  if (prefixlength == 2U)
  {
    for (code = 0; code <= maxcode; code++)
    {
      rightchar = calcbucketboundsparts(&bucketspec,
                                        bcktab,
                                        code,
                                        maxcode,
                                        partwidth,
                                        rightchar,
                                        numofchars);
      accubucketsize += bucketspec.nonspecialsinbucket;
      if (rightchar == 0)
      {
        bucketspec2->subbuckettab[currentchar][numofchars-1].bucketend 
          = accubucketsize;
        accubucketsize += bucketspec.specialsinbucket;
        bucketspec2->superbuckettab[currentchar].bucketend = accubucketsize;
        currentchar++;
      } else
      {
        gt_assert(bucketspec.specialsinbucket == 0);
        bucketspec2->subbuckettab[currentchar][rightchar-1].bucketend 
          = accubucketsize;
      }
    }
  } else
  {
    Codetype *basepower;

    basepower = initbasepower(numofchars,prefixlength);
    for (code = 0; code < 16; code++)
    {
      printf("code=%u = %lu\n",(unsigned int) code,
             (unsigned long) expandtwocharcode(code,prefixlength,
                                               basepower,bcktab));
    }
    gt_free(basepower);
    exit(EXIT_SUCCESS);
  }
  resetsorted(bucketspec2);
  for (idx = 0; idx<numofchars; idx++)
  {
    bucketspec2->order[idx] = idx;
  }
  gt_qsort_r(bucketspec2->order,(size_t) numofchars,
             sizeof (*bucketspec2->order),bucketspec2,
             comparesuperbucketsizes);
  return bucketspec2;
}

static void forwardderive(const Bucketspec2 *bucketspec2,
                          const Seqpos *suftab,
                          unsigned long *targetptr,
                          const Encodedsequence *encseq,
                          Readmode readmode,
                          unsigned int source,
                          unsigned long idx)
{
  Seqpos startpos;
  GtUchar cc;

  gt_assert (idx < targetptr[source]);
  for (; idx < targetptr[source]; idx++)
  {
    startpos = suftab[idx];
    if (startpos > 0)
    {
      cc = getencodedchar(encseq,startpos-1,readmode);
      if (ISNOTSPECIAL(cc) && !bucketspec2->superbuckettab[cc].sorted)
      {
        gt_assert(suftab[targetptr[cc]] == startpos - 1);
        targetptr[cc]++;
      }
    }
  }
}

static void backwardderive(const Bucketspec2 *bucketspec2,
                           const Seqpos *suftab,
                           unsigned long *targetptr,
                           const Encodedsequence *encseq,
                           Readmode readmode,
                           unsigned int source,
                           unsigned long idx)
{
  Seqpos startpos;
  GtUchar cc;

  gt_assert (idx > targetptr[source]);
  for (idx--; idx >= targetptr[source]; idx--)
  {
    startpos = suftab[idx];
    if (startpos > 0)
    {
      cc = getencodedchar(encseq,startpos-1,readmode);
      if (ISNOTSPECIAL(cc) && !bucketspec2->superbuckettab[cc].sorted)
      {
        targetptr[cc]--;
        gt_assert(suftab[targetptr[cc]] == startpos - 1);
      }
    }
  }
}

void gt_copysortsuffixes(const Bucketspec2 *bucketspec2, const Seqpos *suftab,
                         const Encodedsequence *encseq, Readmode readmode)
{
  unsigned long hardwork = 0, *targetptr;
  unsigned int idx, idxsource, source, second;

#ifdef WITHSUFFIXES
  {
  const Seqpos *ptr;
  for (ptr = suftab; ptr < suftab + bucketspec2->partwidth; ptr++)
  {
    showsequenceatstartpos(stdout,
                           ISDIRREVERSE(readmode) ? false : true,
                           ISDIRCOMPLEMENT(readmode) ? true : false,
                           encseq,
                           *ptr);
  }
  }
#endif
  source = bucketspec2->order[0];
  targetptr = gt_malloc(sizeof(*targetptr) * bucketspec2->numofchars);
  for (idxsource = 0; idxsource<bucketspec2->numofchars; idxsource++)
  {
    source = bucketspec2->order[idxsource];
    for (second = 0; second < bucketspec2->numofchars; second++)
    {
      if (!bucketspec2->subbuckettab[source][second].sorted && source != second)
      {
        printf("hard work for %u %u\n",source,second);
        hardwork += getendidx(bucketspec2,source,second) -
                    getstartidx(bucketspec2,source,second);
        bucketspec2->subbuckettab[source][second].sorted = true;
      }
    }
    bucketspec2->superbuckettab[source].sorted = true;
    if (getstartidx(bucketspec2,source,0) <
        getstartidx(bucketspec2,source,source))
    {
      for (idx = 0; idx < bucketspec2->numofchars; idx++)
      {
        targetptr[idx] = getstartidx(bucketspec2,idx,source);
      }
      forwardderive(bucketspec2,
                    suftab,
                    targetptr,
                    encseq,
                    readmode,
                    source,
                    getstartidx(bucketspec2,source,0));
    }
    if (getendidx(bucketspec2,source,source) <
        getendidx(bucketspec2,source,bucketspec2->numofchars))
    {
      for (idx = 0; idx < bucketspec2->numofchars; idx++)
      {
        targetptr[idx] = getendidx(bucketspec2,idx,source);
      }
      backwardderive(bucketspec2,
                     suftab,
                     targetptr,
                     encseq,
                     readmode,
                     source,
                     getendidx(bucketspec2,source,bucketspec2->numofchars));
    }
    for (idx = 0; idx < bucketspec2->numofchars; idx++)
    {
      bucketspec2->subbuckettab[idx][source].sorted = true;
    }
  }
  printf("# hardwork = %lu (%.2f)\n",
            hardwork,(double) hardwork/getencseqtotallength(encseq));
  gt_free(targetptr);
}

void bucketspec2_delete(Bucketspec2 *bucketspec2)
{
  gt_assert(bucketspec2 != NULL);
  gt_array2dim_delete(bucketspec2->subbuckettab);
  gt_free(bucketspec2->superbuckettab);
  gt_free(bucketspec2->order);
  gt_free(bucketspec2);
}
