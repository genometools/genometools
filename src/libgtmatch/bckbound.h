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

#ifndef BCKBOUND_H
#define BCKBOUND_H

#include <assert.h>
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "seqpos-def.h"
#include "sfx-codespec.h"
#include "intcode-def.h"

typedef struct
{
  Seqpos left;
  unsigned long nonspecialsinbucket,
                specialsinbucket;
} Bucketspecification;

typedef struct
{
  Seqpos offset,
         left,
         right;
} Lcpinterval;

typedef struct
{
  Seqpos *leftborder,
         *countspecialcodes;
  Codetype numofallcodes,
           numofspecialcodes,
           **multimappower,
           *basepower,
           *filltable;
  unsigned long **distpfxidx;
} Bcktab;

int mapbcktab(Bcktab *bcktab,
              const Str *indexname,
              unsigned int numofchars,
              unsigned int prefixlength,
              Error *err);

void freebcktab(Bcktab *bcktab);

void initbcktabwithNULL(Bcktab *bcktab);

int bcktab2file(FILE *fp,
                const Bcktab *bcktab,
                Error *err);

/*@unused@*/ static inline unsigned int
              calcbucketboundaries(Bucketspecification *bucketspec,
                                   const Bcktab *bcktab,
                                   Codetype code,
                                   Codetype maxcode,
                                   Seqpos totalwidth,
                                   unsigned int rightchar,
                                   unsigned int numofchars)
{
  bucketspec->left = bcktab->leftborder[code];
  if (code == maxcode)
  {
    assert(totalwidth >= bucketspec->left);
    bucketspec->nonspecialsinbucket
      = (unsigned long) (totalwidth - bucketspec->left);
  } else
  {
    if (bcktab->leftborder[code+1] > 0)
    {
      bucketspec->nonspecialsinbucket
        = (unsigned long) (bcktab->leftborder[code+1] - bucketspec->left);
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
  }
  assert(rightchar == code % numofchars);
  if (rightchar == numofchars - 1)
  {
    bucketspec->specialsinbucket
      = (unsigned long)
        bcktab->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)];
    if (bucketspec->nonspecialsinbucket >= bucketspec->specialsinbucket)
    {
      bucketspec->nonspecialsinbucket -= bucketspec->specialsinbucket;
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
    rightchar = 0;
  } else
  {
    bucketspec->specialsinbucket = 0;
    rightchar++;
  }
  return rightchar;
}

#endif
