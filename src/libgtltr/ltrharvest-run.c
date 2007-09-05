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
#include "libgtmatch/pos2seqnum.pr"

#include "ltrharvest-opt.h"
#include "repeats.h"
#include "searchforLTRs.h"
#include "duplicates.h"
#include "outputstd.c"

static int runltrharvest(LTRharvestoptions *lo, Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  Seqpos *markpos = NULL;

  env_error_check(env);

  /* map suffix array */
  if(streamsuffixarray(&suffixarray,
                       &totallength,
		       SARR_LCPTAB | SARR_SUFTAB | SARR_ESQTAB,
		       lo->str_indexname,
		       false,
		       env) != 0)
  {
    return -1;
  }

  /* test if motif is valid and encode motif */
  if( testmotifandencodemotif(&lo->motif, suffixarray.alpha, env) != 0)
  {
    return -1; 
  }

  /* show defined option and values */
  if(lo->verbosemode)
  {
    showuserdefinedoptionsandvalues(lo);
  }

  // calculate markpos array for contig offset
  markpos = calculatemarkpositions(suffixarray.encseq, 
                                   suffixarray.numofdbsequences, 
                                   env);
  if(markpos == NULL)
  { 
    return -1;
  }

  /* init array for maximal repeats */
  INITARRAY (&lo->repeatinfo.repeats, Repeat);
  lo->repeatinfo.suffixarrayptr = &suffixarray;

  /* search for maximal repeats */ 
  if(enumeratemaxpairs(&suffixarray,
                       (uint32_t)lo->minseedlength,
		       (void*)simpleexactselfmatchstore,
		       lo,
		       env) != 0)
  {
    return -1;
  }

  /* init array for candidate pairs */
  INITARRAY(&lo->arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if(searchforLTRs (&suffixarray, lo, markpos, env) != 0)
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
  
  /* print inner region multiple FASTA file of predictions */
  // fehlt noch
  
  /* print GFF3 format file of predictions */
  /*if(lo->gff3output)
  {
    if(printgff3format(lo,
	  suffixarray,
	  markpos) != 0 )
    {
      return -1; 
    }
  }*/

  FREESPACE(markpos);

  /* print predictions to stdout */
  showinfoiffoundfullLTRs(lo, &suffixarray, env);

  /* free prediction array */
  FREEARRAY(&lo->arrayLTRboundaries, LTRboundaries);
  /* free suffixarray */
  freesuffixarray(&suffixarray, env);

  return 0;
}

int parseargsandcallltrharvest(int argc,const char *argv[],Env *env)
{
  LTRharvestoptions lo;
  int retval;
  bool haserr = false;

  lo.env = env; //for getting env in simpleexactselfmatchstore 

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
