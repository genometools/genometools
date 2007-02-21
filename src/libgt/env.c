/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "env.h"
#include "xansi.h"

struct Env {
  MA *ma;
  Error *error;
};

Env* env_new(void)
{
  Env *e = xmalloc(sizeof (Env));
  e->ma = ma_new();
  e->error = error_new();
  return e;
}

MA* env_ma(const Env *e)
{
  assert(e && e->ma);
  return e->ma;
}

Error* env_error(const Env *e)
{
  assert(e && e->error);
  return e->error;
}

void env_delete(Env *e)
{
  if (!e) return;
  ma_delete(e->ma);
  error_delete(e->error);
  free(e);
}

void env_error_set(Env *e, const char *format, ...)
{
  va_list ap;
  assert(e && format);
  va_start(ap, format);
  error_vset(env_err(e), format, ap);
  va_end(ap);
}
