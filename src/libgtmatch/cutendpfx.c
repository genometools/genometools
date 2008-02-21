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
#include "libgtcore/symboldef.h"
#include "libgtcore/chardef.h"
#include "intcode-def.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "bcktab.h"
#include "cutendpfx.h"

struct Bucketenumerator
{
  const Bcktab *bcktab; /* only need multimappower and filltable */
  unsigned int prefixlength,
               demandprefixlength;
  Bucketspecification bucketspec;
  Codetype currentcode, lastcode;
};

static inline unsigned int prefixqgram2code(Codetype *code,
                                            const Codetype **multimappower,
                                            unsigned int qvalue,
                                            unsigned int qvalueprefix,
                                            const Uchar *qgram)
{
  int i;
  Codetype tmpcode = 0;
  Uchar a;

  assert(qvalueprefix > 0);
  assert(qvalue > qvalueprefix);
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

Bucketenumerator *newbucketenumerator(const Bcktab *bcktab,
                                      unsigned int prefixlength,
                                      const Uchar *demandprefix,
                                      unsigned int demandprefixlength)
{
  Bucketenumerator *bucketenumerator;
  unsigned int firstspecial;

  ALLOCASSIGNSPACE(bucketenumerator,NULL,Bucketenumerator,1);
  bucketenumerator->bcktab = bcktab;
  bucketenumerator->prefixlength = prefixlength;
  bucketenumerator->demandprefixlength = demandprefixlength;
  firstspecial = prefixqgram2code(&bucketenumerator->currentcode,
                                  bcktab_multimappower(bcktab),
                                  prefixlength,
                                  demandprefixlength,
                                  demandprefix);
  assert(firstspecial == prefixlength);
  bucketenumerator->lastcode
    = bucketenumerator->currentcode
      + bcktab_filltable(bcktab,demandprefixlength);
  return bucketenumerator;
}

bool nextbucketenumerator(Lcpinterval *itv,Bucketenumerator *bucketenumerator)
{
  while (true)
  {
    if (bucketenumerator->currentcode > bucketenumerator->lastcode)
    {
      break;
    }
    calcbucketboundaries(&bucketenumerator->bucketspec,
                         bucketenumerator->bcktab,
                         bucketenumerator->currentcode);
    bucketenumerator->currentcode++;
    if (bucketenumerator->bucketspec.nonspecialsinbucket > 0)
    {
      itv->left = bucketenumerator->bucketspec.left;
      itv->right = itv->left +
                   bucketenumerator->bucketspec.nonspecialsinbucket - 1;
      itv->offset = (Seqpos) bucketenumerator->demandprefixlength;
      return true;
    }
  }
  return false;
}

void freebucketenumerator(Bucketenumerator **bucketenumerator)
{
  FREESPACE(*bucketenumerator);
}
