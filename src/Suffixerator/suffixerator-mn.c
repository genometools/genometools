#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "env-errout.h"

#include "runsuffixerator.pr"

int main(int argc,const char *argv[])
{
  Env *env;
  bool haserr = false;

  env = env_new();
  if (parseargsandcallsuffixerator(argc,argv,env) != 0)
  {
    ENVERROUT;
    haserr = true;
  }
  if (env_delete(env))
  {
    haserr = true;
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
