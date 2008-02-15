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
#include "seqpos-def.h"
#include "sfx-codespec.h"
#include "intcode-def.h"

typedef struct
{
  Seqpos left;
  unsigned long nonspecialsinbucket,
                specialsinbucket;
} Bucketboundaries;

/*@unused@*/ static inline unsigned int
              calcbucketboundaries(Bucketboundaries *bbound,
                                   const Seqpos *leftborder,
                                   const Seqpos *countspecialcodes,
                                   Codetype code,
                                   Codetype maxcode,
                                   Seqpos totalwidth,
                                   unsigned int rightchar,
                                   unsigned int numofchars)
{
  bbound->left = leftborder[code];
  if (code == maxcode)
  {
    assert(totalwidth >= bbound->left);
    bbound->nonspecialsinbucket = (unsigned long) (totalwidth - bbound->left);
  } else
  {
    if (leftborder[code+1] > 0)
    {
      bbound->nonspecialsinbucket
        = (unsigned long) (leftborder[code+1] - bbound->left);
    } else
    {
      bbound->nonspecialsinbucket = 0;
    }
  }
  assert(rightchar == code % numofchars);
  if (rightchar == numofchars - 1)
  {
    bbound->specialsinbucket
      = (unsigned long)
        countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)];
    if (bbound->nonspecialsinbucket >= bbound->specialsinbucket)
    {
      bbound->nonspecialsinbucket -= bbound->specialsinbucket;
    } else
    {
      bbound->nonspecialsinbucket = 0;
    }
    rightchar = 0;
  } else
  {
    bbound->specialsinbucket = 0;
    rightchar++;
  }
  return rightchar;
}

#endif
