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

#include "spacedef.h"
#include "turnwheels.h"

#define MAXNUMOFWHEELS 18

 struct Turningwheel
{
  unsigned int wheelspace[MAXNUMOFWHEELS],
               numofwheels,
               asize,
               idx,
               minchanged;
};

Turningwheel *newTurningwheel(unsigned int numofwheels,
                              unsigned int asize)
{
  unsigned int i;
  Turningwheel *tw;

  ALLOCASSIGNSPACE(tw,NULL,Turningwheel,1);
  assert(numofwheels < (unsigned int) MAXNUMOFWHEELS);
  assert(numofwheels > 0);
  assert(asize > 0);
  for (i=0; i<numofwheels; i++)
  {
    tw->wheelspace[i] = 0;
  }
  tw->numofwheels = numofwheels;
  tw->asize = asize;
  tw->idx = numofwheels - 1;
  tw->minchanged = 0;
  return tw;
}

bool nextTurningwheel(Turningwheel *tw)
{
  while (true)
  {
    tw->wheelspace[tw->idx]++;
    tw->minchanged = tw->idx;
    if (tw->wheelspace[tw->idx] == tw->asize)
    {
      tw->wheelspace[tw->idx] = 0;
      if (tw->idx == 0)
      {
        return false;
      }
      tw->idx--;
    } else
    {
      tw->idx = tw->numofwheels-1;
      break;
    }
  }
  return true;
}

unsigned int minchangedTurningwheel(const Turningwheel *tw)
{
  return tw->minchanged;
}

void outputTurningwheel(const Turningwheel *tw)
{
  unsigned int i;

  for (i=0; i<tw->numofwheels; i++)
  {
    printf("%u",tw->wheelspace[i]);
  }
}

void freeTurningwheel(Turningwheel **tw)
{
  FREESPACE(*tw);
}
