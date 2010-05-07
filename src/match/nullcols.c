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

#include <stdlib.h>
#include "core/types_api.h"
#include "core/ma_api.h"

typedef struct
{
  unsigned long *endindex,
                *positions;
} Charatpos;

Charatpos *gt_newCharatpos(unsigned long patternlength,unsigned int alphasize)
{
  Charatpos *catpos;

  catpos = gt_malloc(sizeof (*catpos));
  catpos->endindex = gt_malloc(sizeof (unsigned long) * alphasize);
  catpos->positions = gt_malloc(sizeof (unsigned long) * patternlength);
  return catpos;
}

Charatpos *gt_reinitCharatpos(Charatpos *catpos,
                           const GtUchar *pattern,unsigned long patternlength,
                           unsigned int alphasize)
{
  const GtUchar *pptr;
  unsigned long partialsum, tmp;
  unsigned int idx;

  for (idx=0; idx<alphasize; idx++)
  {
    catpos->endindex[idx] = 0;
  }
  for (pptr=pattern; pptr<pattern+ patternlength; pptr++)
  {
    catpos->endindex[(int) *pptr]++;
  }
  partialsum = catpos->endindex[0];
  catpos->endindex[0] = 0;
  for (idx=1U; idx<alphasize; idx++)
  {
    tmp = catpos->endindex[idx];
    catpos->endindex[idx] = partialsum;
    partialsum += tmp;
  }
  for (pptr=pattern; pptr<pattern+patternlength; pptr++)
  {
    catpos->positions[catpos->endindex[*pptr]++]
      = (unsigned long) (pptr - pattern);
  }
  return catpos;
}

void gt_wrapCharatpos(Charatpos **catposptr)
{
  Charatpos *catpos = *catposptr;
  gt_free(catpos->endindex);
  gt_free(catpos->positions);
  gt_free(catpos);
  *catposptr = NULL;
}

void gt_maintainnullcols(const Charatpos *catpos,
                      unsigned long *front0,GtUchar cc,unsigned long depth)
{
  unsigned long idx;

  for (idx = (cc == 0) ? 0 : catpos->endindex[cc-1];
       idx < catpos->endindex[cc]; idx++)
  {
    unsigned long pos = catpos->positions[idx];
    if (front0[pos] == depth)
    {
      front0[pos]++;
    }
  }
}
