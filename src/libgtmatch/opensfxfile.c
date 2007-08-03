/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"

/*@null@*/ FILE *opensfxfile(const Str *indexname,
                             const char *suffix,
                             const char *mode,
                             Env *env)
{
  Str *tmpfilename;
  FILE *fp;

  env_error_check(env);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,suffix,env);
  fp = env_fa_fopen(env,str_get(tmpfilename),mode);
  if (fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s",str_get(tmpfilename),
                                                    strerror(errno));
  }
  str_delete(tmpfilename,env);
  return fp;
}

bool indexfilealreadyexists(const Str *indexname,const char *suffix,Env *env)
{
  struct stat statbuf;
  Str *tmpfilename;

  env_error_check(env);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,suffix,env);

  if(stat(str_get(tmpfilename),&statbuf) == 0)
  {
    str_delete(tmpfilename,env);
    return true;
  }
  str_delete(tmpfilename,env);
  return false;
}
