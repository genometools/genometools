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
  int had_err;
  env = env_new();
  had_err = tool(argc, (const char**) argv, env);
  if (env_error_is_set(env)) {
    fprintf(stderr, "error: %s\n", env_error_get(env));
    assert(had_err);
  }
  if (env_delete(env))
    return 2; /* programmer error */
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
