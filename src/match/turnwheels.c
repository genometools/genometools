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

#include "core/ma.h"
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

Turningwheel *gt_turningwheel_new(unsigned int numofwheels,
                                  unsigned int asize)
{
  unsigned int i;
  Turningwheel *tw;

  tw = gt_malloc(sizeof (*tw));
  gt_assert(numofwheels < (unsigned int) MAXNUMOFWHEELS);
  gt_assert(numofwheels > 0);
  gt_assert(asize > 0);
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

size_t gt_turningwheel_size(void)
{
  return sizeof (Turningwheel);
}

bool gt_turningwheel_next(Turningwheel *tw)
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

unsigned int gt_turningwheel_minchanged(const Turningwheel *tw)
{
  return tw->minchanged;
}

void gt_turningwheel_output(const Turningwheel *tw)
{
  unsigned int i;

  for (i=0; i<tw->numofwheels; i++)
  {
    printf("%u",tw->wheelspace[i]);
  }
}

void gt_turningwheel_delete(Turningwheel *tw)
{
  if (tw != NULL)
  {
    gt_free(tw);
  }
}
