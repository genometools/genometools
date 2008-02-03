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

#include <limits.h>
#include <assert.h>
#include <math.h>
#include "spacedef.h"

unsigned int *initbasepower(unsigned int base,unsigned int len)
{
  unsigned int thepower = 1U, i, minfailure, *basepower;

  ALLOCASSIGNSPACE(basepower,NULL,unsigned int,len+1);
  minfailure = UINT_MAX/base;
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
      fprintf(stderr,"overflow when computing %u * %u",thepower,base);
      exit(EXIT_FAILURE); /* programming error */
    }
    thepower *= base;
  }
  return basepower;
}

unsigned int *initfilltable(const unsigned int *basepower,unsigned int len)
{
  unsigned int i, *filltable;

  ALLOCASSIGNSPACE(filltable,NULL,unsigned int,len);
  for (i=0; i<len; i++)
  {
    filltable[i] = basepower[len-i]-1;
  }
  return filltable;
}

unsigned int *initmappower(unsigned int numofchars,unsigned int qvalue)
{
  unsigned int *mappower,
               mapindex,
               thepower = (unsigned int) pow((double) numofchars,
                                             (double) (qvalue-1));

  assert(numofchars > 0);
  ALLOCASSIGNSPACE(mappower,NULL,unsigned int,numofchars);
  mappower[0] = 0;
  for (mapindex = 1U; mapindex < numofchars; mapindex++)
  {
    mappower[mapindex] = mappower[mapindex-1] + thepower;
  }
  return mappower;
}
