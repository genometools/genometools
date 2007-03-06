/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgt/mutate.h>

Seq* mutate(const char *description, const char *orig_seq, unsigned long len,
            Alpha *alpha, unsigned int rate, Env *env)
{
  env_error_check(env);
  assert(description && orig_seq && alpha);
  assert(rate >= 0 && rate <= 100);

  /* XXX */
  assert(0);
  return NULL;
}
