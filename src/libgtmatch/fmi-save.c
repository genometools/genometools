/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <math.h>

#include "fmindex.h"

#include "fmi-mapspec.pr"
#include "opensfxfile.pr"

static int writefmascii (const Str *indexname,
                         const Fmindex *fm,
                         bool storeindexpos,
                         Env *env)
{
  FILE *fmafp;

  env_error_check(env);
  if ((fmafp = opensfxfile (indexname, FMASCIIFILESUFFIX,"wb",env)) == NULL)
  {
    return -1;
  }
  fprintf (fmafp, "bwtlength=" FormatSeqpos "\n",
           PRINTSeqposcast(fm->bwtlength));
  fprintf (fmafp, "longest=" FormatSeqpos "\n",
                   PRINTSeqposcast(fm->longestsuffixpos));
  fprintf (fmafp, "storeindexpos=%d\n", storeindexpos ? 1 : 0);
  fprintf (fmafp, "log2blocksize=%u\n", (unsigned int) fm->log2bsize);
  fprintf (fmafp, "log2markdist=%u\n", (unsigned int) fm->log2markdist);
  fprintf (fmafp, "specialcharacters=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.specialcharacters));
  fprintf (fmafp, "specialranges=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.specialranges));
  fprintf (fmafp, "lengthofspecialprefix=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.lengthofspecialprefix));
  fprintf (fmafp, "lengthofspecialsuffix=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.lengthofspecialsuffix));
  fprintf (fmafp, "suffixlength=%u\n", (unsigned int) fm->suffixlength);
  env_fa_xfclose(fmafp, env);
  return 0;
}

static int writefmdata (const Str *indexname,
                        Fmindex *fm, 
                        bool storeindexpos,
                        Env *env)
{
  FILE *fp;

  env_error_check(env);
  if ((fp = opensfxfile (indexname, FMDATAFILESUFFIX,"wb",env)) == NULL)
  {
    return -1;
  }
  if(flushfmindex2file(fp,fm,storeindexpos,env) != 0)
  {
    return -2;
  }
  env_fa_xfclose(fp, env);
  return 0;
}

int saveFmindex (const Str *indexname,Fmindex *fm,
                 bool storeindexpos,Env *env)
{
  env_error_check(env);
  if (writefmascii (indexname, fm, storeindexpos,env) != 0)
  {
    return -1;
  }
  if (writefmdata (indexname, fm, storeindexpos,env) != 0)
  {
    return -2;
  }
  return 0;
}
