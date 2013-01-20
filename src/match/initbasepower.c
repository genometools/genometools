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

#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/codetype.h"
#include "core/unused_api.h"

unsigned int gt_maxbasepower(unsigned int numofchars)
{
  GtCodetype minfailure, thepower = (GtCodetype) 1;
  unsigned int i;

  minfailure = (~(GtCodetype) 0)/(GtCodetype) numofchars;
  for (i=0; thepower < minfailure; i++)
  {
    thepower *= numofchars;
  }
  return i;
}

GtCodetype *gt_initbasepower(unsigned int numofchars,unsigned int prefixlength)
{
  GtCodetype thepower = (GtCodetype) 1, *basepower;
  unsigned int i;

  basepower = gt_malloc(sizeof *basepower * (prefixlength+1));
  for (i=0; /* Nothing */; i++)
  {
    basepower[i] = thepower;
    if (i == prefixlength)
    {
      break;
    }
    gt_assert(thepower < ((~(GtCodetype) 0)/(GtCodetype) numofchars));
    thepower *= numofchars;
  }
  return basepower;
}

GtCodetype *gt_filllargestchartable(unsigned int numofchars,
                               unsigned int kmersize)
{
  GtCodetype code, *ptr, *filltable;

  filltable = gt_malloc(sizeof *filltable * kmersize);
  code = (GtCodetype) numofchars;
  for (ptr = filltable + kmersize - 1; ptr >= filltable; ptr--)
  {
    *ptr = code-1;
    code *= numofchars;
  }
  return filltable;
}

GtCodetype *gt_initfilltable(unsigned int numofchars,unsigned int prefixlength)
{
  unsigned int i;
  GtCodetype *filltable, *basepower;

  basepower = gt_initbasepower(numofchars,prefixlength);
  filltable = gt_malloc(sizeof *filltable * prefixlength);
  for (i=0; i<prefixlength; i++)
  {
    filltable[i] = basepower[prefixlength-i]-1;
  }
  gt_free(basepower);
  return filltable;
}

GtCodetype **gt_initmultimappower(unsigned int numofchars,unsigned int qvalue)
{
  int offset;
  unsigned int mapindex;
  GtCodetype thepower = (GtCodetype) 1, *mmptr, **multimappower;

  gt_array2dim_malloc(multimappower,(unsigned long) qvalue,numofchars);
  for (offset=(int) (qvalue - 1); offset>=0; offset--)
  {
    mmptr = multimappower[offset];
    mmptr[0] = 0;
    for (mapindex = 1U; mapindex < numofchars; mapindex++)
    {
      mmptr[mapindex] = mmptr[mapindex-1] + thepower;
    }
    thepower *= numofchars;
  }
  return multimappower;
}

void gt_multimappower_delete(GtCodetype **multimappower)
{
  if (multimappower != NULL)
  {
    gt_array2dim_delete(multimappower);
  }
}
