/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "types.h"
#include "divmodmul.h"
#include "spacedef.h"
#include "intcode-def.h"
#include "partssuf-def.h"

typedef struct
{
  Codetype nextcode;
  Seqpos widthofpart;
  Uint64 suftaboffset, 
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

Suftabparts *newsuftabparts(unsigned int numofparts,
                            const Seqpos *leftborder,
                            Codetype numofallcodes,
                            Uint64 numofsuffixestoinsert,
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
    if (numofsuffixestoinsert < (Uint64) numofparts)
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
    unsigned int part;
    Uint remainder,
         widthofsuftabpart,
         suftaboffset = 0;
    Uint64 sumofwidth = 0;
    ALLOCASSIGNSPACE(suftabparts->components,NULL,Suftabpartcomponent,
                     numofparts);
    CHECKUint64Cast(numofsuffixestoinsert/(Uint64) numofparts);
    widthofsuftabpart = (Uint) (numofsuffixestoinsert/(Uint64) numofparts);
    remainder = (Uint) (numofsuffixestoinsert % (Uint64) numofparts);
    suftabparts->largestwidth = 0;
    printf("# numofsuffixestoinsert=%lu\n",(Showuint) numofsuffixestoinsert);
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
          = (Uint64) leftborder[suftabparts->components[part-1].nextcode];
      }
      printf("widthofpart[%u]=%lu\n",part,
              (Showuint) suftabparts->components[part].widthofpart);
      if (suftabparts->largestwidth <
         suftabparts->components[part].widthofpart)
      {
        suftabparts->largestwidth
          = suftabparts->components[part].widthofpart;
      }
      sumofwidth += (Uint64) suftabparts->components[part].widthofpart;
      suftabparts->components[part].sumofwidth = sumofwidth;
    }
    assert(sumofwidth == numofsuffixestoinsert);
  }
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

Uint64 stpgetcurrentsuftaboffset(unsigned int part,
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

Uint64 stpgetcurrentsumofwdith(unsigned int part,
                               const Suftabparts *suftabparts)
{
  return suftabparts->components[part].sumofwidth;
}

Seqpos stpgetcurrentwidtofpart(unsigned int part,const Suftabparts *suftabparts)
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
