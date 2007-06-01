#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include "types.h"
#include "spacedef.h"
#include "libgtcore/env.h"

#include "dstrdup.pr"
#include "compfilenm.pr"

typedef struct
{
  char *filename;
  FILE *fpin;
} Fopeninfo;

static int callfileopen(void *info,const char *path,Env *env)
{
  Fopeninfo *fileopeninfo = (Fopeninfo *) info;
  char *tmpfilename;

  env_error_check(env);
  tmpfilename = COMPOSEFILENAMEGENERIC(path,'/',fileopeninfo->filename);
  fileopeninfo->fpin = env_fa_fopen(env,tmpfilename,"rb");
  if (fileopeninfo->fpin != NULL)
  {
    FREESPACE(tmpfilename);
    return 1;
  }
  FREESPACE(tmpfilename);
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
                                  const char *filename,
                                  Env *env)
{
  Fopeninfo fileopeninfo;

  env_error_check(env);
  fileopeninfo.fpin = env_fa_fopen(env,filename,"rb");
  if (fileopeninfo.fpin != NULL)
  {
    return fileopeninfo.fpin;
  }
  if (strchr(filename,'/') != NULL)
  {
    env_error_set(env,
                  "filename \"%s\" contains illegal symbol '/': the path list "
                  "specified by environment variable \"%s\" cannot be searched "
                  "for it",filename,envstring);
    return NULL;
  }
  ASSIGNDYNAMICSTRDUP(fileopeninfo.filename,filename);
  if (evalpathlist(envstring,(void *) &fileopeninfo,callfileopen,env) != 0)
  {
    FREESPACE(fileopeninfo.filename);
    return NULL;
  }
  if (fileopeninfo.fpin == NULL)
  {
    FREESPACE(fileopeninfo.filename);
    env_error_set(env,"cannot find file \"%s\"",filename);
    return NULL;
  }
  FREESPACE(fileopeninfo.filename);
  return fileopeninfo.fpin;
}
