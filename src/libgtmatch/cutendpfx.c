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
#include "bckbound.h"
#include "cutendpfx.h"

struct Bucketenumerator
{
  const Seqpos *bcktab,
               *countspecialcodes;
  Codetype numofallcodes;
  const Codetype *filltable,
                 **multimappower;
  unsigned int prefixlength,
               numofchars,
               demandprefixlength;
  Seqpos totallength;
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

Bucketenumerator *newbucketenumerator(Seqpos totallength,
                                      const Seqpos *bcktab,
                                      const Seqpos *countspecialcodes,
                                      Codetype numofallcodes,
                                      const Codetype **multimappower,
                                      const Codetype *filltable,
                                      unsigned int prefixlength,
                                      const Uchar *demandprefix,
                                      unsigned int demandprefixlength,
                                      unsigned int numofchars)
{
  Bucketenumerator *bucketenumerator;
  unsigned int firstspecial;

  ALLOCASSIGNSPACE(bucketenumerator,NULL,Bucketenumerator,1);
  bucketenumerator->totallength = totallength;
  bucketenumerator->bcktab = bcktab;
  bucketenumerator->countspecialcodes = countspecialcodes;
  bucketenumerator->numofallcodes = numofallcodes;
  bucketenumerator->prefixlength = prefixlength;
  bucketenumerator->demandprefixlength = demandprefixlength;
  bucketenumerator->multimappower = multimappower;
  bucketenumerator->filltable = filltable;
  bucketenumerator->numofchars = numofchars;
  firstspecial = prefixqgram2code(&bucketenumerator->currentcode,
                                  bucketenumerator->multimappower,
                                  prefixlength,
                                  demandprefixlength,
                                  demandprefix);
  assert(firstspecial == prefixlength);
  bucketenumerator->lastcode = bucketenumerator->currentcode +
                               bucketenumerator->filltable[demandprefixlength];
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
    (void) calcbucketboundaries(&bucketenumerator->bucketspec,
                                bucketenumerator->bcktab,
                                bucketenumerator->countspecialcodes,
                                bucketenumerator->currentcode,
                                bucketenumerator->numofallcodes,
                                bucketenumerator->totallength,
                                bucketenumerator->currentcode %
                                  bucketenumerator->numofchars,
                                bucketenumerator->numofchars);
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
