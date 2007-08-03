/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <string.h>
#include "libgtcore/env.h"
#include "readmode-def.h"

static char *readmodes[] = {"fwd", 
                            "rev", 
                            "cpl", 
                            "rcl"};

const char *showreadmode(Readmode readmode)
{
  return readmodes[(int) readmode];
}

int parsereadmode(const char *dirargstring,Env *env)
{
  size_t i;

  env_error_check(env);
  for(i=0; i<sizeof(readmodes)/sizeof(readmodes[0]); i++)
  {
    if(strcmp(dirargstring,readmodes[i]) == 0)
    {
      return (int) i;
    }
  }
  env_error_set(env,"argument to option -dir must be fwd or rev or cpl or rcl");
  return -1;
}
