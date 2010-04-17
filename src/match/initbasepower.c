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

#include "core/assert_api.h"
#include "intcode-def.h"
#include "spacedef.h"

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
  GtCodetype thepower = (GtCodetype) 1, minfailure, *basepower;
  unsigned int i;

  ALLOCASSIGNSPACE(basepower,NULL,GtCodetype,prefixlength+1);
  minfailure = (~(GtCodetype) 0)/(GtCodetype) numofchars;
  for (i=0; /* Nothing */; i++)
  {
    basepower[i] = thepower;
    if (i == prefixlength)
    {
      break;
    }
    gt_assert(thepower < minfailure);
    thepower *= numofchars;
  }
  return basepower;
}

GtCodetype *gt_filllargestchartable(unsigned int numofchars,
                               unsigned int kmersize)
{
  GtCodetype code, *ptr, *filltable;

  filltable = gt_malloc(sizeof (GtCodetype) * kmersize);
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
  ALLOCASSIGNSPACE(filltable,NULL,GtCodetype,prefixlength);
  for (i=0; i<prefixlength; i++)
  {
    filltable[i] = basepower[prefixlength-i]-1;
  }
  FREESPACE(basepower);
  return filltable;
}

#define ARRAY2DIMMALLOC(ARRAY2DIM, ROWS, COLUMNS, TYPE)\
        {\
          unsigned int rownumber;\
          ALLOCASSIGNSPACE(ARRAY2DIM,NULL,TYPE *,ROWS);\
          ALLOCASSIGNSPACE((ARRAY2DIM)[0],NULL,TYPE,(ROWS) * (COLUMNS));\
          for (rownumber = 1U; rownumber < (ROWS); rownumber++)\
          {\
            (ARRAY2DIM)[rownumber] = (ARRAY2DIM)[rownumber-1] + (COLUMNS);\
          }\
        }

GtCodetype **gt_initmultimappower(unsigned int numofchars,unsigned int qvalue)
{
  int offset;
  unsigned int mapindex;
  GtCodetype thepower = (GtCodetype) 1, *mmptr, **multimappower;

  ARRAY2DIMMALLOC(multimappower,qvalue,numofchars,GtCodetype);
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

void gt_multimappowerfree(GtCodetype ***multimappower)
{
  if (*multimappower != NULL)
  {
    FREESPACE((*multimappower)[0]);
    FREESPACE((*multimappower));
    *multimappower = NULL;
  }
}
