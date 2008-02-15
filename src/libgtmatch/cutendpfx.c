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
#include <errno.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>
#include "libgtcore/symboldef.h"
#include "libgtcore/chardef.h"
#include "intcode-def.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "bckbound.h"

typedef struct
{
  const Seqpos *bcktab,
               *countspecialcodes;
  Codetype numofallcodes;
  const Codetype **multimappower;
  unsigned int prefixlength,
               demandprefixlength;
  Bucketspecification bucketspec;
  Codetype currentcode;
} Bucketenumerator;

static inline unsigned int extendqgram2code(Uchar extendchar,
                                            Codetype *code,
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
  for (i=(int) (qvalue-1); i>=(int) qvalueprefix; i--)
  {
    tmpcode += multimappower[i][extendchar];
  }
  for (i=(int) (qvalueprefix-1); i>=0; i--)
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

Bucketenumerator *newbucketenumerator(const Seqpos *bcktab,
                                      const Seqpos *countspecialcodes,
                                      Codetype numofallcodes,
                                      const Codetype **multimappower,
                                      unsigned int prefixlength,
                                      const Uchar *demandprefix,
                                      unsigned int demandprefixlength)
{
  Bucketenumerator *bucketenumerator;
  unsigned int firstspecial;

  ALLOCASSIGNSPACE(bucketenumerator,NULL,Bucketenumerator,1);
  bucketenumerator->bcktab = bcktab;
  bucketenumerator->countspecialcodes = countspecialcodes;
  bucketenumerator->numofallcodes = numofallcodes;
  bucketenumerator->prefixlength = prefixlength;
  bucketenumerator->demandprefixlength = demandprefixlength;
  bucketenumerator->multimappower = multimappower;
  firstspecial = extendqgram2code(0,
                                  &bucketenumerator->currentcode,
                                  bucketenumerator->multimappower,
                                  prefixlength,
                                  demandprefixlength,
                                  demandprefix);
  assert(firstspecial == demandprefixlength);
  return bucketenumerator;
}
