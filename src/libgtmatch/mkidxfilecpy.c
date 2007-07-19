/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include <string.h>
#include <inttypes.h>
#include "libgtcore/str.h"
#include "libgtcore/env.h"

int makeindexfilecopy(const Str *destindex,
                      const Str *sourceindex,
                      const char *suffix,
                      uint64_t maxlength,
                      Env *env)
{
  Str *destfilename = NULL, *sourcefilename = NULL;
  FILE *fpdest = NULL, *fpsource = NULL;
  int cc;
  bool haserr = false;

  destfilename = str_clone(destindex,env);
  str_append_cstr(destfilename,suffix,env);
  fpdest = env_fa_fopen(env,str_get(destfilename),"wb");
  if (fpdest == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s",str_get(destfilename),
                                                    strerror(errno));
    haserr = true;
  }
  if(!haserr)
  {
    sourcefilename = str_clone(sourceindex,env);
    str_append_cstr(sourcefilename,suffix,env);
    fpsource = env_fa_fopen(env,str_get(sourcefilename),"rb");
    if (fpsource == NULL)
    {
      env_error_set(env,"cannot open file \"%s\": %s",str_get(sourcefilename),
                                                      strerror(errno));
      haserr = true;
    }
  }
  if(!haserr)
  {
    if(maxlength == 0)
    {
      while((cc = fgetc(fpsource)) != EOF)
      {
        (void) putc(cc,fpdest);
      }
    } else
    {
      uint64_t pos;
  
      for(pos = 0; pos < maxlength; pos++)
      {
        if((cc = fgetc(fpsource)) == EOF)
        {
          break;
        }
        (void) putc(cc,fpdest);
      }
    }
  }
  env_fa_xfclose(fpdest,env);
  env_fa_xfclose(fpsource,env);
  str_delete(destfilename,env);
  str_delete(sourcefilename,env);
  return haserr ? -1 : 0;
}
