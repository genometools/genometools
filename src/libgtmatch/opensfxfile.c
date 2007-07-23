/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include <string.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"

/*@null@*/ FILE *opensfxfile(const Str *indexname,
                             const char *suffix,
                             const char *mode,Env *env)
{
  Str *tmpfilename;
  FILE *fp;

  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,suffix,env);
  fp = env_fa_fopen(env,str_get(tmpfilename),mode);
  if (fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s\n",str_get(tmpfilename),
                                                      strerror(errno));
  }
  str_delete(tmpfilename,env);
  return fp;
}
