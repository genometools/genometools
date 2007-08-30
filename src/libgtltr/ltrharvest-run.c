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
#include "duplicates.h"
#include "outputstd.c"

static int runltrharvest(LTRharvestoptions *lo, Env *env)
{
  bool haserror = false;
  Suffixarray suffixarray;
  Seqpos totallength;
  
  env_error_check(env);

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

  /* test if motif is valid */
  if( testmotif(&lo->motif, suffixarray.alpha, env) != 0)
  {
    return -1; 
  }

  /* show defined option and values */
  if(lo->verbosemode)
  {
    showuserdefinedoptionsandvalues(lo);
  }

  /* init array for maximal repeats */
  INITARRAY (&lo->repeatinfo.repeats, Repeat);
  lo->repeatinfo.suffixarrayptr = &suffixarray;

  /* search for maximal repeats */ 
  if(!haserror && enumeratemaxpairs(&suffixarray,
                       (uint32_t)lo->minseedlength,
		       (void*)simpleexactselfmatchstore,
		       lo,
		       env) != 0)
  {
    haserror = true;
  }

  /* init array for candidate pairs */
  INITARRAY(&lo->arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if(searchforLTRs (&suffixarray, lo, env) != 0)
  {
     return -1; 
  }

  /* free array for maximal repeats */
  FREEARRAY(&lo->repeatinfo.repeats, Repeat);

  /* remove exact duplicates */
  removeduplicates(&lo->arrayLTRboundaries);

  /* remove overlapping predictions if desired */
  // function fehlt noch

  /* print multiple FASTA file of predictions */
  // fehlt noch

  /* print GFF3 format file of predictions */
  // fehlt noch

  /* print predictions to stdout */
  showinfoiffoundfullLTRs(lo, &suffixarray, env);

  /* free prediction array */
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
    printargsline(argv,argc);
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
