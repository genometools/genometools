/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include "types.h"
#include "spacedef.h"
#include "megabytes.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"

#include "alphabet.pr"
#include "makeprj.pr"
#include "sfx-opt.pr"
#include "measure-time.pr"
#include "outprjfile.pr"
#include "suffixerator.pr"
#include "opensfxfile.pr"

typedef struct
{
  FILE *outfpsuftab,
       *outfpbwttab;
} Outfileinfo;

static int suftab2file(void *info,
                       const Seqpos *suftab,
                       Seqpos widthofpart,
                       Env *env)
{
  Outfileinfo *outfileinfo = (Outfileinfo *) info;

  env_error_check(env);
  if (outfileinfo->outfpsuftab != NULL)
  {
    if (fwrite(suftab,
              sizeof (*suftab),
              (size_t) widthofpart,
              outfileinfo->outfpsuftab)
              != (size_t) widthofpart)
    {
      env_error_set(env,"cannot write %lu items of size %lu: errormsg=\"%s\"",
           (Showuint) widthofpart,
           (Showuint) sizeof (*suftab),
           strerror(errno));
      return -1;
    }
  }
  if (outfileinfo->outfpbwttab != NULL)
  {
    printf("# simulating output of %lu elements of the bwttab\n",
            (Showuint) widthofpart);
  }
  return 0;
}

static int outal1file(const Str *indexname,const Alphabet *alpha,Env *env)
{
  FILE *al1fp;
  bool haserr = false;

  env_error_check(env);
  al1fp = opensfxfile(indexname,".al1",env);
  if(al1fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    outputalphabet(al1fp,alpha);
    env_fa_xfclose(al1fp,env);
  }
  return haserr ? -1 : 0;
}

static int runsuffixerator(const Suffixeratoroptions *so,Env *env)
{
  unsigned char numofchars = 0;
  unsigned long numofsequences;
  Uint64 totallength;
  Alphabet *alpha;
  Specialcharinfo specialcharinfo;
  PairSeqpos *filelengthtab = NULL;
  Outfileinfo outfileinfo;
  bool haserr = false;
  Encodedsequence *encseq;
  Measuretime *mtime;

  env_error_check(env);
  inittheclock(&mtime,
               "determining sequence length and number of special symbols",
               env);
  outfileinfo.outfpsuftab = NULL;
  outfileinfo.outfpbwttab = NULL;
  alpha = assigninputalphabet(so->isdna,
                              so->isprotein,
                              so->str_smap,
                              strarray_get(so->filenametab,0),
                              env);
  if (alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (scanfastasequence(&numofsequences,
                          &totallength,
                          &specialcharinfo,
                          so->filenametab,
                          &filelengthtab,
                          getsymbolmapAlphabet(alpha), env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    assert(so->prefixlength > 0);
    if (outprjfile(so->str_indexname,
                   so->filenametab,
                   filelengthtab,
                   totallength,
                   numofsequences,
                   &specialcharinfo,
                   so->prefixlength,
                   env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    numofchars = (unsigned char) getnumofcharsAlphabet(alpha);
    if (outal1file(so->str_indexname,alpha,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if(so->outsuftab)
    {
      outfileinfo.outfpsuftab = opensfxfile(so->str_indexname,".suf",env);
      if(outfileinfo.outfpsuftab == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    if(so->outbwttab)
    {
      outfileinfo.outfpbwttab 
        = opensfxfile(so->str_indexname,".bwt",env);
      if(outfileinfo.outfpbwttab == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    deliverthetime(stdout,mtime,"computing sequence encoding",env);
    encseq = initencodedseq(true,
                            so->filenametab,
                            NULL,
                            (Uint64) totallength,
                            &specialcharinfo,
                            alpha,
                            str_length(so->str_sat) > 0
                                  ? str_get(so->str_sat)
                                  : NULL,
                            env);
    if (encseq == NULL)
    {
      haserr = true;
    } else
    {
      if (so->outtistab)
      {
        if (flushencseqfile(so->str_indexname,encseq,env) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (!haserr)
  {
    printf("# specialcharacters=%lu\n",
           (Showuint) specialcharinfo.specialcharacters);
    printf("# specialranges=%lu\n",
           (Showuint) specialcharinfo.specialranges);
    if (so->outsuftab || so->outbwttab)
    {
      if (suffixerator(suftab2file,
                       &outfileinfo,
                       totallength,
                       specialcharinfo.specialcharacters,
                       specialcharinfo.specialranges,
                       encseq,
                       numofchars,
                       so->prefixlength,
                       so->numofparts,
                       mtime,
                       env) != 0)
      {
        haserr = true;
      }
    }
  }
  env_fa_xfclose(outfileinfo.outfpsuftab,env);
  env_fa_xfclose(outfileinfo.outfpbwttab,env);
  FREESPACE(filelengthtab);
  freeAlphabet(&alpha,env);
  freeEncodedsequence(&encseq,env);
  deliverthetime(stdout,mtime,NULL,env);
  return haserr ? -1 : 0;
}

int parseargsandcallsuffixerator(int argc,const char *argv[],Env *env)
{
  Suffixeratoroptions so;
  int retval;
  bool haserr = false;

  retval = suffixeratoroptions(&so,argc,argv,env);
  if (retval == 0)
  {
    if (runsuffixerator(&so,env) < 0)
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
  wrapsfxoptions(&so,env);
  return haserr ? -1 : 0;
}
