/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "gtr.h"

/* The GenomeTools (gt) genome analysis system */
int main(int argc, char *argv[])
{
  Env *env;
  GTR *gtr;
  int rval;
  env = env_new();
  env_error_set_progname(env, argv[0]);
  gtr = gtr_new(env);
  gtr_register_components(gtr);
  switch (gtr_parse(gtr, &rval, argc, (const char**) argv, env)) {
    case OPTIONPARSER_OK:
      argc -= rval;
      argv += rval;
      rval = gtr_run(gtr, argc, (const char**) argv, env);
      break;
    case OPTIONPARSER_ERROR:
      rval = 1; /* user error */
      break;
    case OPTIONPARSER_REQUESTS_EXIT:
      rval = 0; /* everything went fine */
  }
  if (env_error_is_set(env)) {
    fprintf(stderr, "%s: error: %s\n", env_error_get_progname(env),
            env_error_get(env));
    assert(rval);
  }
  gtr_delete(gtr, env);
  if (env_delete(env))
    return 2; /* programmer error */
  return rval;
}
