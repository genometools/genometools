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

#include <assert.h>
#include "spacedef.h"
#include "intcode-def.h"

Codetype *initbasepower(unsigned int base,unsigned int len)
{
  unsigned int i;
  Codetype thepower = (Codetype) 1, minfailure, *basepower;

  ALLOCASSIGNSPACE(basepower,NULL,Codetype,len+1);
  minfailure = (~(Codetype) 0)/(Codetype) base;
  for (i=0; /* Nothing */; i++)
  {
    basepower[i] = thepower;
    if (i == len)
    {
      break;
    }
    if (thepower >= minfailure)
    {
      FREESPACE(basepower);
      fprintf(stderr,"overflow when computing %lu * %u",
              (unsigned long) thepower,base);
      exit(EXIT_FAILURE); /* programming error */
    }
    thepower *= base;
  }
  return basepower;
}

Codetype *initfilltable(const Codetype *basepower,unsigned int len)
{
  unsigned int i;
  Codetype *filltable;

  ALLOCASSIGNSPACE(filltable,NULL,Codetype,len);
  for (i=0; i<len; i++)
  {
    filltable[i] = basepower[len-i]-1;
  }
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

Codetype **initmultimappower(unsigned int numofchars,unsigned int qvalue)
{
  int offset;
  unsigned int mapindex;
  Codetype thepower = (Codetype) 1, *mmptr, **multimappower;

  ARRAY2DIMMALLOC(multimappower,qvalue,numofchars,unsigned int);
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

void multimappowerfree(Codetype ***multimappower)
{
  if (*multimappower != NULL)
  {
    FREESPACE((*multimappower)[0]);
    FREESPACE((*multimappower));
    *multimappower = NULL;
  }
}
