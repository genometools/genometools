/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"

typedef struct
{
  Str *filename;
  FILE *fpin;
} Fopeninfo;

static int callfileopen(void *info,const Str *path,Env *env)
{
  Fopeninfo *fileopeninfo = (Fopeninfo *) info;
  Str *tmpfilename;

  env_error_check(env);
  tmpfilename = str_clone(path,env);
  str_append_char(tmpfilename,'/',env);
  str_append_str(tmpfilename,fileopeninfo->filename,env); 
  fileopeninfo->fpin = env_fa_fopen(env,str_get(tmpfilename),"rb");
  if (fileopeninfo->fpin != NULL)
  {
    str_delete(tmpfilename,env);
    return 1;
  }
  str_delete(tmpfilename,env);
  return 0;
}

static int evalpathlist(const char *envarname,void *info,
                        int(*applypath)(void *,const Str *,Env *),
                        Env *env)
{
  const char *envptr = getenv(envarname);

  env_error_check(env);
  if (envptr != NULL)
  {
    unsigned long currentpos, startpos;
    Str *start = str_new(env);
    int ret;

    startpos = 0;
    for (currentpos=0; /* Nothing */ ; currentpos++)
    {
      if (envptr[currentpos] == '\0')
      {
        str_reset(start);
        str_append_cstr_nt(start,envptr + startpos,currentpos - startpos,env);
        if (applypath(info,start,env) < 0)
        {
          str_delete(start,env);
          return -1;
        }
        break;
      }
      if (envptr[currentpos] == ':')
      {
        str_reset(start);
        str_append_cstr_nt(start,envptr + startpos,currentpos - startpos,env);
        ret = applypath(info,start,env);
        startpos = currentpos + 1;
        if (ret < 0)
        {
          str_delete(start,env);
          return -2;
        }
        if (ret == 1) /* stop */
        {
          break;
        }
      }
    }
    str_delete(start,env);
  }
   return 0;
}

/*@null@*/ FILE *scanpathsforfile(const char *envstring,
                                  const Str *filename,
                                  Env *env)
{
  Fopeninfo fileopeninfo;

  env_error_check(env);
  fileopeninfo.fpin = env_fa_fopen(env,str_get(filename),"rb");
  if (fileopeninfo.fpin != NULL)
  {
    return fileopeninfo.fpin;
  }
  if (strchr(str_get(filename),'/') != NULL)
  {
    env_error_set(env,
                  "filename \"%s\" contains illegal symbol '/': the path list "
                  "specified by environment variable \"%s\" cannot be searched "
                  "for it",str_get(filename),envstring);
    return NULL;
  }
  fileopeninfo.filename = str_clone(filename,env);
  if (evalpathlist(envstring,(void *) &fileopeninfo,callfileopen,env) != 0)
  {
    str_delete(fileopeninfo.filename,env);
    return NULL;
  }
  if (fileopeninfo.fpin == NULL)
  {
    str_delete(fileopeninfo.filename,env);
    env_error_set(env,"cannot find file \"%s\"",str_get(filename));
    return NULL;
  }
  str_delete(fileopeninfo.filename,env);
  return fileopeninfo.fpin;
}
