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
  Env *env = xcalloc(1, sizeof (Env));
  env->ma = ma_new();
  ma_init(env->ma, env);
  env->error = error_new(env->ma);
  return env;
}

MA* env_ma(const Env *env)
{
  assert(env && env->ma);
  return env->ma;
}

Error* env_error(const Env *env)
{
  assert(env && env->error);
  return env->error;
}

Log* env_log(const Env *env)
{
  assert(env);
  return env->log;
}

void env_set_log(Env *env, Log *log)
{
  assert(env);
  env->log = log;
}

int env_delete(Env *env)
{
  int rval;
  assert(env);
  log_delete(env->log, env->ma);
  error_delete(env->error, env->ma);
  rval = ma_check_space_leak(env->ma, env);
  ma_clean(env->ma, env);
  ma_delete(env->ma);
  free(env);
  return rval;
}

void env_ma_free_func(void *ptr, Env *env)
{
  assert(env);
  if (!ptr) return;
  ma_free(ptr, env_ma(env));
}

void env_error_set(Env *env, const char *format, ...)
{
  va_list ap;
  assert(env && format);
  va_start(ap, format);
  error_vset(env_error(env), format, ap);
  va_end(ap);
}

void env_log_log(Env *env, const char *format, ...)
{
  va_list ap;
  assert(env && format);
  if (!env_log(env)) return;
  va_start(ap, format);
  log_vlog(env_log(env), format, ap);
  va_end(ap);
}
