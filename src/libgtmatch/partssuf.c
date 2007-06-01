/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "types.h"
#include "divmodmul.h"
#include "spacedef.h"
#include "partssuf-def.h"

typedef struct
{
  Uint nextcode,
       widthofpart;
  Uint64 suftaboffset, sumofwidth;
} Suftabpartcomponent;

 struct _Suftabparts
{
  Suftabpartcomponent *components;
  unsigned int numofparts;
  Uint largestwidth;
};

static Uint findfirstlarger(const Uint *leftborder,
                            Uint numofallcodes,
                            Uint suftaboffset)
{
  Uint l= 0, r = numofallcodes, mid, found = numofallcodes;

  while (l+1 < r)
  {
    mid = DIV2(l+r);
    if (suftaboffset == leftborder[mid])
    {
      return mid;
    }
    if (suftaboffset < leftborder[mid])
    {
      found = mid;
      r = mid - 1;
    } else
    {
      l = mid + 1;
    }
  }
  return found;
}

Suftabparts *newsuftabparts(unsigned int numofparts,
                            const Uint *leftborder,
                            Uint numofallcodes,
                            Uint64 numofsuffixestoinsert,
                            Uint fullspecials,
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
                     (Uint) numofparts);
    CHECKIFFITS32BITS(numofsuffixestoinsert/(Uint64) numofparts);
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
      printf("widthofpart[%u]=%u\n",part,
              suftabparts->components[part].widthofpart);
      if (suftabparts->largestwidth <
         suftabparts->components[part].widthofpart)
      {
        suftabparts->largestwidth
          = suftabparts->components[part].widthofpart;
      }
      sumofwidth += (Uint64) suftabparts->components[part].widthofpart;
      suftabparts->components[part].sumofwidth = sumofwidth;
    }
    printf("sumofwidth = %u, numofsuffixestoinsert = %u\n",(Uint) sumofwidth,
                        (Uint) numofsuffixestoinsert);
    assert(sumofwidth == numofsuffixestoinsert);
  }
  return suftabparts;
}

Uint stpgetcurrentmincode(Uint part,const Suftabparts *suftabparts)
{
  if (part == 0)
  {
    return 0;
  }
  return suftabparts->components[part-1].nextcode + 1;
}

Uint64 stpgetcurrentsuftaboffset(Uint part,const Suftabparts *suftabparts)
{
  return suftabparts->components[part].suftaboffset;
}

Uint stpgetcurrentmaxcode(Uint part,const Suftabparts *suftabparts)
{
  if (part == suftabparts->numofparts - 1)
  {
    return suftabparts->components[part].nextcode - 1;
  }
  return suftabparts->components[part].nextcode;
}

Uint64 stpgetcurrentsumofwdith(Uint part,const Suftabparts *suftabparts)
{
  return suftabparts->components[part].sumofwidth;
}

Uint stpgetcurrentwidtofpart(Uint part,const Suftabparts *suftabparts)
{
  return suftabparts->components[part].widthofpart;
}

Uint stpgetlargestwidth(const Suftabparts *suftabparts)
{
  return suftabparts->largestwidth;
}

Uint stpgetnumofparts(const Suftabparts *suftabparts)
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
