/*
  Copyright (C) 2007 by David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
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

#include "repeattypes.h"

#include "repeats.h"
#include "ltrharvest-opt.h"

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
		       simpleexactselfmatchoutput,
		       NULL,
		       env) != 0)
  {
    haserror = true;
  }

  FREEARRAY(&lo->repeatinfo.repeats, Repeat);
  /* free suffixarray */
  freesuffixarray(&suffixarray, env);

  return haserror ? -1 : 0;
}

int parseargsandcallltrharvest(int argc,const char *argv[],Env *env)
{
  LTRharvestoptions lo;
  int retval;
  bool haserr = false;
  //unsigned int spacepeak;

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
  wrapltrharvestoptions(&lo,env);
  return haserr ? -1 : 0;
}
