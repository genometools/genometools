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
#include "bcktab.h"
#include "sfx-copysort.h"

struct Subbucketspec
{
  unsigned int numofchars, *order;
  unsigned long *bucketends,
                **subbucket;
};

static unsigned long superbucketsize(const Subbucketspec *subbucketspec,
                                     unsigned int bucketnum)
{
  if (bucketnum == 0)
  {
    return subbucketspec->bucketends[0];
  }
  return subbucketspec->bucketends[bucketnum] -
         subbucketspec->bucketends[bucketnum-1];
}

static int comparesuperbucketsizes(const void *a,const void *b,void *data)
{
  const Subbucketspec *subbucketspec = (const Subbucketspec *) data;

  unsigned long size1 = superbucketsize(subbucketspec,
                                        *(const unsigned int *) a);
  unsigned long size2 = superbucketsize(subbucketspec,
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

Subbucketspec *subbuckets_new(const Bcktab *bcktab,
                              Seqpos partwidth,
                              unsigned int numofchars)
{
  Codetype code, maxcode;
  Bucketspecification bucketspec;
  Subbucketspec *subbucketspec;
  unsigned int idx, rightchar = 0, currentchar = 0;
  unsigned long accubucketsize = 0;

  gt_assert(numofchars > 0);
  subbucketspec = gt_malloc(sizeof(*subbucketspec));
  subbucketspec->numofchars = numofchars;
  subbucketspec->order = gt_malloc(sizeof(*subbucketspec->order) *
                                   numofchars);
  subbucketspec->bucketends = gt_malloc(sizeof(*subbucketspec->bucketends) *
                                        numofchars);
  subbucketspec->subbucket = gt_malloc(sizeof(*subbucketspec->subbucket) *
                                       numofchars);
  subbucketspec->subbucket[0] = gt_malloc(sizeof(**subbucketspec->subbucket) *
                                          (numofchars * numofchars));
  for (idx = 1U; idx<numofchars; idx++)
  {
    subbucketspec->subbucket[idx] = subbucketspec->subbucket[idx-1] +
                                    numofchars;
  }
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
    subbucketspec->subbucket[currentchar]
                            [(rightchar == 0) ? (numofchars-1) : (rightchar-1)]
                            = accubucketsize;
    printf("subbucket[%u][%u]=%lu\n",
             currentchar,
             (unsigned int) (rightchar == 0) ? (numofchars-1) : (rightchar-1),
             accubucketsize);
    if (rightchar == 0)
    {
      accubucketsize += bucketspec.specialsinbucket;
      subbucketspec->bucketends[currentchar] = accubucketsize;
      printf("bucketends[%u]=%lu\n", currentchar, accubucketsize);
      currentchar++;
    }
  }
  for (idx = 0; idx<numofchars; idx++)
  {
    subbucketspec->order[idx] = idx;
    printf("superbucketsize[%u]=%lu\n",
            idx,superbucketsize(subbucketspec,idx));
  }
  gt_qsort_r(subbucketspec->order,(size_t) numofchars,
             sizeof (*subbucketspec->order),subbucketspec,
             comparesuperbucketsizes);
  for (idx = 0; idx<numofchars; idx++)
  {
    printf("bucket %u: size %lu\n",subbucketspec->order[idx],
            superbucketsize(subbucketspec,subbucketspec->order[idx]));
  }
  return subbucketspec;
}

void subbuckets_delete(Subbucketspec *subbucketspec)
{
  gt_free(subbucketspec->subbucket[0]);
  gt_free(subbucketspec->subbucket);
  gt_free(subbucketspec->bucketends);
  gt_free(subbucketspec->order);
  gt_free(subbucketspec);
}
