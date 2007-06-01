/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include "libgtcore/env.h"
#include "types.h"
#include "spacedef.h"

#include "compfilenm.pr"

/*@null@*/ FILE *opensfxfile(const char *indexname,
                             const char *suffix,Env *env)
{
  char *tmpfilename;
  FILE *fp;

  tmpfilename = COMPOSEFILENAME(indexname,suffix);
  fp = env_fa_fopen(env,tmpfilename,"wb");
  if (fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s\n",tmpfilename,
                                                      strerror(errno));
  }
  FREESPACE(tmpfilename);
  return fp;
}

