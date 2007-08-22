/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <string.h>
#include <inttypes.h>
#include "libgtcore/str.h"
#include "libgtcore/env.h"

#include "opensfxfile.pr"

int makeindexfilecopy(const Str *destindex,
                      const Str *sourceindex,
                      const char *suffix,
                      uint64_t maxlength,
                      Env *env)
{
  FILE *fpdest = NULL, *fpsource = NULL;
  int cc;
  bool haserr = false;

  env_error_check(env);
  fpdest = opensfxfile(destindex,suffix,"wb",env);
  if (fpdest == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    fpsource = opensfxfile(sourceindex,suffix,"rb",env);
    if (fpsource == NULL)
    {
      haserr = true;
    }
  }
  printf("# cp %s%s %s%s\n",
           str_get(sourceindex),suffix,str_get(destindex),suffix);
  if (!haserr)
  {
    if (maxlength == 0)
    {
      while ((cc = fgetc(fpsource)) != EOF)
      {
        (void) putc(cc,fpdest);
      }
    } else
    {
      uint64_t pos;

      for (pos = 0; pos < maxlength; pos++)
      {
        if ((cc = fgetc(fpsource)) == EOF)
        {
          break;
        }
        (void) putc(cc,fpdest);
      }
    }
  }
  env_fa_xfclose(fpdest,env);
  env_fa_xfclose(fpsource,env);
  return haserr ? -1 : 0;
}
