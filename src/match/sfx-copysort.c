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
#include "bcktab.h"
#include "sfx-copysort.h"

struct Subbucketspec
{
  unsigned long *bucketends,
                **subbucket;
};

Subbucketspec *subbuckets_new(const Bcktab *bcktab,
                              Seqpos partwidth,
                              unsigned int numofchars)
{
  Codetype code, maxcode;
  Bucketspecification bucketspec;
  unsigned int idx, rightchar = 0;
  unsigned long countbucketsize = 0;
  unsigned long currentchar = 0;
  Subbucketspec *subbucketspec;

  gt_assert(numofchars > 0);
  subbucketspec = gt_malloc(sizeof(*subbucketspec));
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
    countbucketsize += bucketspec.nonspecialsinbucket;
    printf("subbucket[%u][%u]=%lu\n",
             (unsigned int) currentchar,
             (unsigned int) (rightchar == 0) ? (numofchars-1) : rightchar-1,
             countbucketsize);
    if (rightchar == 0)
    {
      countbucketsize += bucketspec.specialsinbucket;
      printf("bucketend[%u]=%lu\n",
             (unsigned int) currentchar,
             countbucketsize);
      currentchar++;
    }
  }
  return subbucketspec;
}

void subbuckets_delete(Subbucketspec *subbucketspec)
{
  gt_free(subbucketspec->subbucket[0]);
  gt_free(subbucketspec->subbucket);
  gt_free(subbucketspec->bucketends);
  gt_free(subbucketspec);
}
