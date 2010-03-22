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

#include "core/divmodmul.h"

#include "spacedef.h"
#include "intcode-def.h"
#include "sfx-partssuf-def.h"

typedef struct
{
  GtCodetype nextcode;
  unsigned long widthofpart,
         suftaboffset,
         sumofwidth;
} Suftabpartcomponent;

 struct Suftabparts
{
  Suftabpartcomponent *components;
  unsigned int numofparts;
  unsigned long largestwidth;
};

static GtCodetype findfirstlarger(const unsigned long *leftborder,
                                GtCodetype numofallcodes,
                                unsigned long suftaboffset)
{
  GtCodetype left = 0, right = numofallcodes, mid, found = numofallcodes;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    if (suftaboffset == leftborder[mid])
    {
      return mid;
    }
    if (suftaboffset < leftborder[mid])
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

static void removeemptyparts(Suftabparts *suftabparts,
                             GtLogger *logger)
{
  if (suftabparts->numofparts > 0)
  {
    unsigned int destpart, srcpart;
    for (destpart = 0, srcpart = 0; srcpart < suftabparts->numofparts;
         srcpart++)
    {
      if (suftabparts->components[srcpart].widthofpart > 0)
      {
        if (destpart < srcpart)
        {
          suftabparts->components[destpart] = suftabparts->components[srcpart];
        }
        destpart++;
      }
    }
    if (destpart < srcpart)
    {
      suftabparts->numofparts -= (srcpart - destpart);
    }
    for (srcpart = 0; srcpart < suftabparts->numofparts; srcpart++)
    {
      gt_assert(suftabparts->components[srcpart].widthofpart > 0);
      gt_logger_log(logger,"widthofpart[%u]=%lu",
                  srcpart,
                  suftabparts->components[srcpart].widthofpart);
    }
  }
}

Suftabparts *newsuftabparts(unsigned int numofparts,
                            const unsigned long *leftborder,
                            GtCodetype numofallcodes,
                            unsigned long numofsuffixestoinsert,
                            unsigned long fullspecials,
                            GtLogger *logger)
{
  Suftabparts *suftabparts;

  ALLOCASSIGNSPACE(suftabparts,NULL,Suftabparts,(size_t) 1);
  gt_assert(suftabparts != NULL);
  if (numofsuffixestoinsert == 0)
  {
    suftabparts->numofparts = 0;
  } else
  {
    if (numofsuffixestoinsert < (unsigned long) numofparts)
    {
      suftabparts->numofparts = 1U;
    } else
    {
      suftabparts->numofparts = numofparts;
    }
  }
  if (suftabparts->numofparts == 0)
  {
    suftabparts->largestwidth = fullspecials/numofparts+1;
    suftabparts->components = NULL;
  } else
  {
    unsigned int part, remainder;
    unsigned long widthofsuftabpart,
           suftaboffset = 0,
           sumofwidth = 0;
    ALLOCASSIGNSPACE(suftabparts->components,NULL,Suftabpartcomponent,
                     numofparts);
    widthofsuftabpart = numofsuffixestoinsert/numofparts;
    remainder =
            (unsigned int) (numofsuffixestoinsert % (unsigned long) numofparts);
    suftabparts->largestwidth = 0;
    for (part=0; part < numofparts; part++)
    {
      if (remainder > 0)
      {
        suftaboffset += widthofsuftabpart + 1;
        remainder--;
      } else
      {
        suftaboffset += widthofsuftabpart;
      }
      if (part == numofparts - 1)
      {
        suftabparts->components[part].nextcode = numofallcodes;
      } else
      {
        suftabparts->components[part].nextcode = findfirstlarger(leftborder,
                                                                 numofallcodes,
                                                                 suftaboffset);
      }
      if (part == 0)
      {
        suftabparts->components[part].widthofpart
          = leftborder[suftabparts->components[part].nextcode];
        suftabparts->components[part].suftaboffset = 0;
      } else
      {
        suftabparts->components[part].widthofpart
          = leftborder[suftabparts->components[part].nextcode] -
            leftborder[suftabparts->components[part-1].nextcode];
        suftabparts->components[part].suftaboffset
          = leftborder[suftabparts->components[part-1].nextcode];
      }
      if (suftabparts->largestwidth <
         suftabparts->components[part].widthofpart)
      {
        suftabparts->largestwidth
          = suftabparts->components[part].widthofpart;
      }
      sumofwidth += suftabparts->components[part].widthofpart;
      suftabparts->components[part].sumofwidth = sumofwidth;
    }
    gt_assert(sumofwidth == numofsuffixestoinsert);
  }
  removeemptyparts(suftabparts,logger);
  return suftabparts;
}

GtCodetype stpgetcurrentmincode(unsigned int part,
                              const Suftabparts *suftabparts)
{
  if (part == 0)
  {
    return 0;
  }
  return suftabparts->components[part-1].nextcode + 1;
}

unsigned long stpgetcurrentsuftaboffset(unsigned int part,
                                   const Suftabparts *suftabparts)
{
  return suftabparts->components[part].suftaboffset;
}

GtCodetype stpgetcurrentmaxcode(unsigned int part,
                              const Suftabparts *suftabparts)
{
  if (part == suftabparts->numofparts - 1)
  {
    return suftabparts->components[part].nextcode - 1;
  }
  return suftabparts->components[part].nextcode;
}

unsigned long stpgetcurrentsumofwdith(unsigned int part,
                                 const Suftabparts *suftabparts)
{
  return suftabparts->components[part].sumofwidth;
}

unsigned long stpgetcurrentwidthofpart(unsigned int part,
                                const Suftabparts *suftabparts)
{
  return suftabparts->components[part].widthofpart;
}

unsigned long stpgetlargestwidth(const Suftabparts *suftabparts)
{
  return suftabparts->largestwidth;
}

unsigned int stpgetnumofparts(const Suftabparts *suftabparts)
{
  return suftabparts->numofparts;
}

void freesuftabparts(Suftabparts *suftabparts)
{
  if (suftabparts != NULL)
  {
    FREESPACE(suftabparts->components);
    FREESPACE(suftabparts);
  }
}
