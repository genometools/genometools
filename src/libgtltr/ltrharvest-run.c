/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-maxpairs.pr"
#include "libgtmatch/sfx-map.pr"

#include "ltrharvest-opt.h"
#include "repeats.h"
#include "searchforLTRs.h"

static int runltrharvest(LTRharvestoptions *lo, Env *env)
{
  bool haserror = false;
  Suffixarray suffixarray;
  Seqpos totallength;

  INITARRAY (&lo->repeatinfo.repeats, Repeat);
  lo->repeatinfo.suffixarrayptr = &suffixarray;

  /* map suffix array */
  if(streamsuffixarray(&suffixarray,
                       &totallength,
		       SARR_LCPTAB | SARR_SUFTAB | SARR_ESQTAB,
		       lo->str_indexname,
		       false,
		       env) != 0)
  {
    haserror = true;
  }

  /* search for maximal repeats */ 
  if(!haserror && enumeratemaxpairs(&suffixarray,
                       (uint32_t)lo->minseedlength,
		       (void*)simpleexactselfmatchstore,
		       lo,
		       env) != 0)
  {
    haserror = true;
  }

  INITARRAY(&lo->arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if(searchforLTRs (&suffixarray, lo, env) != 0)
  {
     return -1; 
  }
  FREEARRAY(&lo->repeatinfo.repeats, Repeat);




  FREEARRAY(&lo->arrayLTRboundaries, LTRboundaries);
  /* free suffixarray */
  freesuffixarray(&suffixarray, env);

  return haserror ? -1 : 0;
}

int parseargsandcallltrharvest(int argc,const char *argv[],Env *env)
{
  LTRharvestoptions lo;
  int retval;
  bool haserr = false;

  lo.env = env;

  retval = ltrharvestoptions(&lo,argc,argv,env);
  if (retval == 0)
  {
    haserr = false;
    if (runltrharvest(&lo,env) < 0)
    {
      haserr = true;
    }
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }

  retval = wrapltrharvestoptions(&lo,env);
  if(retval != 0)
  {
    haserr = true;
  }

  return haserr ? -1 : 0;
}
