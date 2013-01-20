/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdbool.h>
#include "core/chardef.h"
#include "core/codetype.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "bcktab.h"
#include "cutendpfx.h"

struct Bucketenumerator
{
  const GtBcktab *bcktab; /* only need multimappower and filltable */
  unsigned int prefixlength,
               demandprefixlength;
  GtBucketspecification bucketspec;
  GtCodetype currentcode, lastcode;
};

static inline unsigned int prefixqgram2code(GtCodetype *code,
                                            const GtCodetype **multimappower,
                                            unsigned int qvalue,
                                            unsigned int qvalueprefix,
                                            const GtUchar *qgram)
{
  int i;
  GtCodetype tmpcode = 0;
  GtUchar a;

  gt_assert(qvalueprefix > 0);
  gt_assert(qvalue > qvalueprefix);
  for (i = (int) (qvalueprefix-1); i>=0; i--)
  {
    a = qgram[i];
    if (ISSPECIAL(a))
    {
      return (unsigned int) i;
    }
    tmpcode += multimappower[i][a];
  }
  *code = tmpcode;
  return qvalue;
}

Bucketenumerator *gt_newbucketenumerator(const GtBcktab *bcktab,
                                      unsigned int prefixlength,
                                      const GtUchar *demandprefix,
                                      unsigned int demandprefixlength)
{
  Bucketenumerator *bucketenumerator;
  GT_UNUSED unsigned int firstspecial;

  bucketenumerator = gt_malloc(sizeof *bucketenumerator);
  bucketenumerator->bcktab = bcktab;
  bucketenumerator->prefixlength = prefixlength;
  bucketenumerator->demandprefixlength = demandprefixlength;
  firstspecial = prefixqgram2code(&bucketenumerator->currentcode,
                                  gt_bcktab_multimappower(bcktab),
                                  prefixlength,
                                  demandprefixlength,
                                  demandprefix);
  gt_assert(firstspecial == prefixlength);
  bucketenumerator->lastcode
    = bucketenumerator->currentcode
      + gt_bcktab_filltable(bcktab,demandprefixlength);
  return bucketenumerator;
}

bool gt_nextbucketenumerator(Lcpinterval *itv,
                             Bucketenumerator *bucketenumerator)
{
  while (true)
  {
    if (bucketenumerator->currentcode > bucketenumerator->lastcode)
    {
      break;
    }
    gt_bcktab_calcboundaries(&bucketenumerator->bucketspec,
                             bucketenumerator->bcktab,
                             bucketenumerator->currentcode);
    bucketenumerator->currentcode++;
    if (bucketenumerator->bucketspec.nonspecialsinbucket > 0)
    {
      itv->left = bucketenumerator->bucketspec.left;
      itv->right = itv->left +
                   bucketenumerator->bucketspec.nonspecialsinbucket - 1;
      itv->offset = (unsigned long) bucketenumerator->demandprefixlength;
      return true;
    }
  }
  return false;
}

void gt_freebucketenumerator(Bucketenumerator *bucketenumerator)
{
  gt_free(bucketenumerator);
}
