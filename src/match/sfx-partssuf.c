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

#include <unistd.h>
#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/divmodmul.h"
#include "core/codetype.h"
#include "bcktab.h"
#include "sfx-partssuf.h"

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
  unsigned long largestwidth,
                numofsuffixestoinsert,
                pagesize;
};

#ifdef SKDEBUG
static void showrecord(const Suftabpartcomponent *component)
{
  printf("# part: width=%lu offset=%lu sumwidth=%lu nextcode=%lu\n",
          component->widthofpart,component->suftaboffset,
                                 component->sumofwidth,
                                 component->nextcode);
}

static void showallrecords(const Suftabparts *suftabparts)
{
  unsigned int idx;

  for (idx = 0; idx < suftabparts->numofparts; idx++)
  {
    showrecord(suftabparts->components + idx);
  }
}
#endif

static void removeemptyparts(Suftabparts *suftabparts,
                             GtLogger *logger)
{
#ifdef SKDEBUG
  printf("# before removial\n");
  showallrecords(suftabparts);
#endif
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
    if (suftabparts->components[suftabparts->numofparts-1].widthofpart == 0)
    {
      gt_assert(suftabparts->numofparts > 1U);
      destpart = suftabparts->numofparts-2;
      while (true)
      {
        if (suftabparts->components[destpart].widthofpart > 0)
        {
          suftabparts->components[destpart].nextcode
            = suftabparts->components[suftabparts->numofparts-1].nextcode;
          suftabparts->numofparts = destpart + 1;
          break;
        }
        if (destpart > 0)
        {
          destpart--;
        } else
        {
          gt_assert(false);
        }
      }
    } else
    {
      if (destpart < srcpart)
      {
        suftabparts->numofparts -= (srcpart - destpart);
        gt_assert(suftabparts->numofparts == destpart);
      }
    }
    for (srcpart = 0; srcpart < suftabparts->numofparts; srcpart++)
    {
      gt_assert(suftabparts->components[srcpart].widthofpart > 0);
      gt_logger_log(logger,"widthofpart[%u]=%lu",
                    srcpart,
                    suftabparts->components[srcpart].widthofpart);
    }
  }
#ifdef SKDEBUG
  printf("#after removal\n");
  showallrecords(suftabparts);
#endif
}

Suftabparts *gt_newsuftabparts(unsigned int numofparts,
                               const GtBcktab *bcktab,
                               unsigned long numofsuffixestoinsert,
                               unsigned long fullspecials,
                               GtLogger *logger)
{
  Suftabparts *suftabparts;

  suftabparts = gt_malloc(sizeof *suftabparts);
  suftabparts->numofsuffixestoinsert = numofsuffixestoinsert;
  suftabparts->pagesize = (unsigned long) sysconf((int) _SC_PAGESIZE);
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
    unsigned long widthofsuftabpart, suftaboffset = 0, sumofwidth = 0;

    suftabparts->components
      = gt_malloc(sizeof (*suftabparts->components) * numofparts);
    widthofsuftabpart = numofsuffixestoinsert/numofparts;
    remainder = (unsigned int) (numofsuffixestoinsert %
                                (unsigned long) numofparts);
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
        suftabparts->components[part].nextcode
          = gt_bcktab_numofallcodes(bcktab);
      } else
      {
        suftabparts->components[part].nextcode
          = gt_bcktab_findfirstlarger(bcktab,suftaboffset);
      }
      if (part == 0)
      {
        suftabparts->components[part].widthofpart
          = gt_bcktab_get(bcktab,suftabparts->components[part].nextcode);
        suftabparts->components[part].suftaboffset = 0;
      } else
      {
        suftabparts->components[part].widthofpart
          = gt_bcktab_get(bcktab,suftabparts->components[part].nextcode) -
            gt_bcktab_get(bcktab,suftabparts->components[part-1].nextcode);
        suftabparts->components[part].suftaboffset
          = gt_bcktab_get(bcktab,suftabparts->components[part-1].nextcode);
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
  return (part == 0) ? 0 : suftabparts->components[part-1].nextcode + 1;
}

unsigned long stpgetcurrentsuftaboffset(unsigned int part,
                                        const Suftabparts *suftabparts)
{
  return suftabparts->components[part].suftaboffset;
}

unsigned long stpgetcurrentleftborderoffset(unsigned int part,
                                            const Suftabparts *suftabparts)
{
  unsigned long mincode = stpgetcurrentmincode(part,suftabparts);

  return (mincode * sizeof (unsigned long) / suftabparts->pagesize)
         * suftabparts->pagesize;
}

GtCodetype stpgetcurrentmaxcode(unsigned int part,
                                const Suftabparts *suftabparts)
{
  return (part == suftabparts->numofparts - 1)
           ? suftabparts->components[part].nextcode - 1
           : suftabparts->components[part].nextcode;
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

unsigned long stpnumofsuffixestoinsert(const Suftabparts *suftabparts)
{
  return suftabparts->numofsuffixestoinsert;
}

void gt_freesuftabparts(Suftabparts *suftabparts)
{
  if (suftabparts != NULL)
  {
    gt_free(suftabparts->components);
    gt_free(suftabparts);
  }
}
