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
  const Seqpos *suftab;
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

Bucketspec2 *bucketspec2_new(const Bcktab *bcktab,
                             Seqpos partwidth,
                             unsigned int numofchars)
{
  Codetype code, maxcode;
  Bucketspecification bucketspec;
  Bucketspec2 *bucketspec2;
  unsigned int idx, rightchar = 0, currentchar = 0;
  unsigned long accubucketsize = 0;

  gt_assert(numofchars > 0);
  bucketspec2 = gt_malloc(sizeof(*bucketspec2));
  bucketspec2->partwidth = partwidth;
  bucketspec2->numofchars = numofchars;
  bucketspec2->order = gt_malloc(sizeof(*bucketspec2->order) *
                                   numofchars);
  bucketspec2->superbuckettab
    = gt_malloc(sizeof(*bucketspec2->superbuckettab) * numofchars);
  gt_array2dim_malloc(bucketspec2->subbuckettab,(unsigned long) numofchars,
                      (unsigned long) numofchars);
  maxcode = bcktab_numofallcodes(bcktab) - 1;
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
    idx = (rightchar == 0) ? (numofchars-1) : (rightchar-1);
    bucketspec2->subbuckettab[currentchar][idx].bucketend = accubucketsize;
    bucketspec2->subbuckettab[currentchar][idx].sorted = false;
    printf("subbucket[%u][%u]=%lu\n", currentchar, idx, accubucketsize);
    if (rightchar == 0)
    {
      accubucketsize += bucketspec.specialsinbucket;
      bucketspec2->superbuckettab[currentchar].bucketend = accubucketsize;
      bucketspec2->superbuckettab[currentchar].sorted = false;
      printf("superbucket[%u].end=%lu\n", currentchar, accubucketsize);
      currentchar++;
    }
  }
  for (idx = 0; idx<numofchars; idx++)
  {
    bucketspec2->order[idx] = idx;
    printf("superbucketsize[%u]=%lu\n",
            idx,superbucketsize(bucketspec2,idx));
  }
  gt_qsort_r(bucketspec2->order,(size_t) numofchars,
             sizeof (*bucketspec2->order),bucketspec2,
             comparesuperbucketsizes);
  for (idx = 0; idx<numofchars; idx++)
  {
    printf("bucket %u: size %lu\n",bucketspec2->order[idx],
            superbucketsize(bucketspec2,bucketspec2->order[idx]));
  }
  return bucketspec2;
}

void gt_copysortsuffixes(Bucketspec2 *bucketspec2, const Seqpos *suftab,
                         const Encodedsequence *encseq, Readmode readmode)
{
  const Seqpos *ptr, *endptr;
  Seqpos startpos;
  unsigned long insertpos;
  GtUchar cc;

  bucketspec2->suftab = suftab;
  endptr = suftab+bucketspec2->superbuckettab[bucketspec2->order[0]].bucketend;
  for (ptr = endptr - 1; ptr >= suftab; ptr--)
  {
    startpos = *ptr;
    if (startpos > 0)
    {
      cc = getencodedchar(encseq,startpos-1,readmode);
      if (ISNOTSPECIAL(cc) && cc > 0)
      {
        bucketspec2->subbuckettab[cc][0].bucketend--;
        insertpos = bucketspec2->subbuckettab[cc][0].bucketend;
        gt_assert(suftab[insertpos] == startpos - 1);
      }
    }
  }
}

void bucketspec2_delete(Bucketspec2 *bucketspec2)
{
  gt_assert(bucketspec2 != NULL);
  gt_array2dim_delete(bucketspec2->subbuckettab);
  gt_free(bucketspec2->superbuckettab);
  gt_free(bucketspec2->order);
  gt_free(bucketspec2);
}
