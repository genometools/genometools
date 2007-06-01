/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "env-errout.h"
#include "sarr-def.h"

#include "alphabet.pr"
#include "sfxmap.pr"
#include "suftaborder.pr"
#include "test-mappedstr.c"
#include "test-encseq.c"
#include "encodedseq.pr"

int main(int argc,const char *argv[])
{
  const char *indexname;
  Env *env;
  bool haserr = false;
  Suffixarray suffixarray;

  if (argc != 2)
  {
    fprintf(stderr,"Usage: %s <indexname>\n",argv[0]);
    return EXIT_FAILURE;
  }
  indexname = argv[1];
  env = env_new();
  if (mapsuffixarray(&suffixarray,true,true,indexname,env) != 0)
  {
    haserr = true;
  }
  if(!haserr)
  {
    if(testencodedsequence((const char **) suffixarray.filenametab,
                           suffixarray.numoffiles,
                           suffixarray.encseq,
                           getsymbolmapAlphabet(suffixarray.alpha),
                           env) != 0)
    {
      haserr = true;
    }
  }
  if(!haserr)
  {
    if(verifymappedstr(&suffixarray,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    checkentiresuftab(suffixarray.encseq,
                      getcharactersAlphabet(suffixarray.alpha),
                      suffixarray.suftab,
                      false, // specialsareequal,
                      true,  // specialsareequalatdepth0,
                      getencseqtotallength(suffixarray.encseq),
                      0,
                      env);
  }
  if (haserr)
  {
    ENVERROUT;
  } else
  {
    freesuffixarray(&suffixarray,env);
  }
  if (env_delete(env))
  {
    haserr = true;
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
