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

