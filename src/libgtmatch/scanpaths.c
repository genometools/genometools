/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include "spacedef.h"
#include "libgtcore/env.h"
#include "libgtcore/str.h"

#include "dstrdup.pr"

typedef struct
{
  Str *filename;
  FILE *fpin;
} Fopeninfo;

static int callfileopen(void *info,const char *path,Env *env)
{
  Fopeninfo *fileopeninfo = (Fopeninfo *) info;
  Str *tmpfilename;

  env_error_check(env);
  tmpfilename = str_new_cstr(path,env);
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
                        int(*applypath)(void *,const char *,Env *),
                        Env *env)
{
  char *envptr = getenv(envarname);

  env_error_check(env);
  if (envptr != NULL)
  {
    char *start, *ptr, *envstring;
    int ret;

    ASSIGNDYNAMICSTRDUP(envstring,envptr);
    start = envstring;
    for (ptr = envstring; /* Nothing */ ; ptr++)
    {
      assert(ptr != NULL);
      if (*ptr == '\0')
      {
        if (applypath(info,start,env) < 0)
        {
          FREESPACE(envstring);
          return -1;
        }
        break;
      }
      if (*ptr == ':')
      {
        *ptr = '\0';
        ret = applypath(info,start,env);
        start = ptr + 1;
        if (ret < 0)
        {
          FREESPACE(envstring);
          return -2;
        }
        if (ret == 1) /* stop */
        {
          break;
        }
      }
    }
    FREESPACE(envstring);
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
