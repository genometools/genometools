/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "types.h"
#include "encseq-def.h"

/* obsolete code */

typedef struct
{
  unsigned int prefixlength;
  Uint64 totallength;
  Uint countspecialmaxprefixlen0;
} CountCompletespecials;

static int lengthofspecialranges(/*@unused@*/ void *info,
                                 const PairUint64 *pair,/*@unused@*/ Env *env)
{
  unsigned int len = (Uint) (pair->uint1 - pair->uint0);
  CountCompletespecials *csp = (CountCompletespecials *) info;

  if (pair->uint0 == 0)
  {
    if (pair->uint1 == csp->totallength)
    {
      csp->countspecialmaxprefixlen0 += (Uint) (len+1);
      printf("add %u completespecials for (%u,%u)\n",
              len+1,(Uint) pair->uint0,(Uint) pair->uint1);
    } else
    {
      csp->countspecialmaxprefixlen0 += (Uint) len;
      printf("add %u completespecials for (%u,%u)\n",
              len,(Uint) pair->uint0,(Uint) pair->uint1);
    }
  } else
  {
    if (pair->uint1 == csp->totallength)
    {
      csp->countspecialmaxprefixlen0 += (Uint) len;
      printf("add %u completespecials for (%u,%u)\n",
              len,(Uint) pair->uint0,(Uint) pair->uint1);
    } else
    {
      if (len >= (unsigned int) 2)
      {
        csp->countspecialmaxprefixlen0 += (Uint) (len - 1);
        printf("add %u completespecials for (%u,%u)\n",
              len-1,(Uint) pair->uint0,(Uint) pair->uint1);
      }
    }
  }
  return 0;
}

Uint determinefullspecials(const Encodedsequence *encseq,
                           Uint64 totallength,
                           unsigned int prefixlength,
                           Env *env)
{
  CountCompletespecials csp;

  csp.countspecialmaxprefixlen0 = 0;
  csp.prefixlength = prefixlength;
  csp.totallength = totallength;
  (void) overallspecialranges(encseq,lengthofspecialranges,&csp,env);
  return csp.countspecialmaxprefixlen0;
}
