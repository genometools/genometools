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
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/esa-maxpairs.pr"
#include "libgtmatch/sfx-map.pr"
#include "libgtmatch/pos2seqnum.pr"

#include "ltrharvest-opt.h"
#include "repeats.h"
#include "searchforLTRs.h"
#include "duplicates.h"
#include "outputstd.h"
#include "outputfasta.h"

static int runltrharvest(LTRharvestoptions *lo, Env *env)
{
  Sequentialsuffixarrayreader *ssar; /* suffix array */
  Seqpos *markpos = NULL;
  unsigned long numofdbsequences;

  env_error_check(env);
  
  ssar = newSequentialsuffixarrayreaderfromfile(lo->str_indexname,
		                  SARR_LCPTAB | SARR_SUFTAB | SARR_ESQTAB,
				  false,
				  env);
  if(ssar == NULL)
  {
    return -1;
  }

  /*
  if(streamsuffixarray(&suffixarray,
                       &totallength,
		       SARR_LCPTAB | SARR_SUFTAB | SARR_ESQTAB,
		       lo->str_indexname,
		       false,
		       env) != 0)
  {
    return -1;
  }
  */

  /* test if motif is valid and encode motif */
  if( testmotifandencodemotif(&lo->motif, 
                              alphabetSequentialsuffixarrayreader(ssar),
			      env) != 0)
  {
    return -1; 
  }

  /* show defined option and values */
  if(lo->verbosemode)
  {
    showuserdefinedoptionsandvalues(lo);
  }

  numofdbsequences = numofdbsequencesSequentialsuffixarrayreader(ssar);
  // calculate markpos array for contig offset
  if( numofdbsequences > 1)
  {
    markpos = encseq2markpositions( 
	encseqSequentialsuffixarrayreader(ssar),
	numofdbsequencesSequentialsuffixarrayreader(ssar), 
	env);
    if(markpos == NULL)
    { 
      return -1;
    }
  }

  /* init array for maximal repeats */
  INITARRAY (&lo->repeatinfo.repeats, Repeat);
  //lo->repeatinfo.suffixarrayptr = &suffixarray;
  lo->repeatinfo.ssarptr = ssar;

  /* search for maximal repeats */ 
  if(enumeratemaxpairs(ssar,
		       getnumofcharsAlphabet(
		         alphabetSequentialsuffixarrayreader(ssar)),
		       encseqSequentialsuffixarrayreader(ssar),
		       readmodeSequentialsuffixarrayreader(ssar),
                       (unsigned int)lo->minseedlength,
		       (void*)simpleexactselfmatchstore,
		       lo,
		       env) != 0)
  {
    return -1;
  }

  /* init array for candidate pairs */
  INITARRAY(&lo->arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if(searchforLTRs (ssar, lo, markpos, env) != 0)
  {
     return -1; 
  }

  /* free array for maximal repeats */
  FREEARRAY(&lo->repeatinfo.repeats, Repeat);

  /* remove exact duplicates */
  removeduplicates(&lo->arrayLTRboundaries);

  /* remove overlapping predictions if desired */
  if(lo->nooverlapallowed || lo->bestofoverlap)
  {
    removeoverlapswithlowersimilarity(&lo->arrayLTRboundaries,
                                      lo->nooverlapallowed);
  }

  /* print multiple FASTA file of predictions */
  if(lo->fastaoutput)
  {
    if (showpredictionsmultiplefasta(lo,
          markpos,
	  false,
	  60,
	  ssar,
	  true,
	  env) != 0)
    {
      return -1;
    }
  }

  /* print inner region multiple FASTA file of predictions */
  if(lo->fastaoutputinnerregion)
  {
    if (showpredictionsmultiplefasta(lo,
          markpos,
	  true,
	  60,
	  ssar,
	  true,
	  env) != 0)
    {
      return -1;
    }
  }

  /* print GFF3 format file of predictions */
  /*if(lo->gff3output)
  {
    if(printgff3format(lo,
	  ssar,
	  markpos) != 0 )
    {
      return -1; 
    }
  }*/

  FREESPACE(markpos);

  /* print predictions to stdout */
  showinfoiffoundfullLTRs(lo, ssar, env);

  /* free prediction array */
  FREEARRAY(&lo->arrayLTRboundaries, LTRboundaries);
  /* free suffixarray */
  freeSequentialsuffixarrayreader(&ssar, env);

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
