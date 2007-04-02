/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/env.h>
#include <libgtcore/tooldriver.h>

int tooldriver(int(*tool)(int argc, const char **argv, Env*),
               int argc, char *argv[])
{
  Env *env;
  int rval;
  env = env_new();
  rval = tool(argc, (const char**) argv, env);
  if (env_error_is_set(env)) {
    fprintf(stderr, "error: %s\n", env_error_get(env));
    assert(rval);
  }
  if (env_delete(env))
    return 2; /* programmer error */
  return rval;
}
