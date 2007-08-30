/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "seqpos-def.h"
#include "divmodmul.h"
#include "spacedef.h"
#include "intcode-def.h"
#include "sfx-partssuf-def.h"

typedef struct
{
  Codetype nextcode;
  Seqpos widthofpart;
  uint64_t suftaboffset,
           sumofwidth;
} Suftabpartcomponent;

 struct _Suftabparts
{
  Suftabpartcomponent *components;
  unsigned int numofparts;
  Seqpos largestwidth;
};

static Codetype findfirstlarger(const Seqpos *leftborder,
                                Codetype numofallcodes,
                                Seqpos suftaboffset)
{
  Codetype left = 0, right = numofallcodes, mid, found = numofallcodes;

  while (left+1 < right)
  {
    mid = DIV2(left+right);
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

static void removeemptyparts(Suftabparts *suftabparts)
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
      assert(suftabparts->components[srcpart].widthofpart > 0);
      /*
      XXX insert later
      printf("# widthofpart[%u]=" FormatSeqpos "\n",
              (unsigned int) srcpart,
              PRINTSeqposcast(suftabparts->components[srcpart].widthofpart));
      */
    }
  }
}

Suftabparts *newsuftabparts(unsigned int numofparts,
                            const Seqpos *leftborder,
                            Codetype numofallcodes,
                            Seqpos numofsuffixestoinsert,
                            Seqpos fullspecials,
                            Env *env)
{
  Suftabparts *suftabparts;

  env_error_check(env);
  ALLOCASSIGNSPACE(suftabparts,NULL,Suftabparts,(size_t) 1);
  assert(suftabparts != NULL);
  if (numofsuffixestoinsert == 0)
  {
    suftabparts->numofparts = 0;
  } else
  {
    if (numofsuffixestoinsert < (Seqpos) numofparts)
    {
      suftabparts->numofparts = (unsigned int) 1;
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
    uint64_t widthofsuftabpart,
             suftaboffset = 0,
             sumofwidth = 0;
    ALLOCASSIGNSPACE(suftabparts->components,NULL,Suftabpartcomponent,
                     numofparts);
    widthofsuftabpart = (uint64_t) numofsuffixestoinsert/numofparts;
    remainder = (unsigned int) (numofsuffixestoinsert % (Seqpos) numofparts);
    suftabparts->largestwidth = 0;
    /* XXX insert later
    printf("# numofsuffixestoinsert=" FormatSeqpos "\n",
           PRINTSeqposcast(numofsuffixestoinsert));
    */
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
          = (uint64_t) leftborder[suftabparts->components[part-1].nextcode];
      }
      if (suftabparts->largestwidth <
         suftabparts->components[part].widthofpart)
      {
        suftabparts->largestwidth
          = suftabparts->components[part].widthofpart;
      }
      sumofwidth += (uint64_t) suftabparts->components[part].widthofpart;
      suftabparts->components[part].sumofwidth = sumofwidth;
    }
    assert(sumofwidth == numofsuffixestoinsert);
  }
  removeemptyparts(suftabparts);
  return suftabparts;
}

Codetype stpgetcurrentmincode(unsigned int part,
                              const Suftabparts *suftabparts)
{
  if (part == 0)
  {
    return 0;
  }
  return suftabparts->components[part-1].nextcode + 1;
}

uint64_t stpgetcurrentsuftaboffset(unsigned int part,
                                   const Suftabparts *suftabparts)
{
  return suftabparts->components[part].suftaboffset;
}

Codetype stpgetcurrentmaxcode(unsigned int part,
                              const Suftabparts *suftabparts)
{
  if (part == suftabparts->numofparts - 1)
  {
    return suftabparts->components[part].nextcode - 1;
  }
  return suftabparts->components[part].nextcode;
}

uint64_t stpgetcurrentsumofwdith(unsigned int part,
                                 const Suftabparts *suftabparts)
{
  return suftabparts->components[part].sumofwidth;
}

Seqpos stpgetcurrentwidtofpart(unsigned int part,
                               const Suftabparts *suftabparts)
{
  return suftabparts->components[part].widthofpart;
}

Seqpos stpgetlargestwidth(const Suftabparts *suftabparts)
{
  return suftabparts->largestwidth;
}

unsigned int stpgetnumofparts(const Suftabparts *suftabparts)
{
  return suftabparts->numofparts;
}

void freesuftabparts(Suftabparts *suftabparts,Env *env)
{
  if (suftabparts != NULL)
  {
    FREESPACE(suftabparts->components);
    FREESPACE(suftabparts);
  }
}
