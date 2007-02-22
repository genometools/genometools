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
  Log *log;
};

Env* env_new(void)
{
  Env *e = xcalloc(1, sizeof (Env));
  e->ma = ma_new();
  e->error = error_new(e->ma);
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

Log* env_log(const Env *e)
{
  assert(e);
  return e->log;
}

void env_set_log(Env *e, Log *log)
{
  assert(e);
  e->log = log;
}

int env_delete(Env *e)
{
  int rval;
  assert(e);
  log_delete(e->log, e->ma);
  error_delete(e->error, e->ma);
  rval = ma_check_space_leak(e->ma);
  ma_delete(e->ma);
  free(e);
  return rval;
}

void env_ma_free(void *ptr, Env *env)
{
  assert(env);
  if (!ptr) return;
  ma_free(ptr, env_ma(env));
}

void env_error_set(Env *e, const char *format, ...)
{
  va_list ap;
  assert(e && format);
  va_start(ap, format);
  error_vset(env_error(e), format, ap);
  va_end(ap);
}

void env_log_log(Env *e, const char *format, ...)
{
  va_list ap;
  assert(e && format);
  if (!env_log(e)) return;
  va_start(ap, format);
  log_vlog(env_log(e), format, ap);
  va_end(ap);
}
