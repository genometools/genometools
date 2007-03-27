/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/array.h>
#include <libgtcore/cstr.h>
#include <libgtcore/strarray.h>

struct StrArray {
  Array *cstrings;
};

StrArray* strarray_new(Env *env)
{
  StrArray *sa;
  env_error_check(env);
  sa = env_ma_malloc(env, sizeof (StrArray));
  sa->cstrings = array_new(sizeof (char*), env);
  return sa;
}

void strarray_add_cstr(StrArray *sa, const char *cstr, Env *env)
{
  char *cstr_copy;
  env_error_check(env);
  assert(sa && cstr);
  cstr_copy = cstr_dup(cstr, env);
  array_add(sa->cstrings, cstr_copy, env);
}

const char* strarray_get(const StrArray *sa, unsigned long strnum)
{
  assert(sa && strnum < array_size(sa->cstrings));
  return *(char**) array_get(sa->cstrings, strnum);
}

unsigned long strarray_size(const StrArray *sa)
{
  assert(sa);
  return array_size(sa->cstrings);
}

void strarray_delete(StrArray *sa, Env *env)
{
  unsigned long i;
  if (!sa) return;
  for (i = 0; i < array_size(sa->cstrings); i++)
    env_ma_free(*(char**) array_get(sa->cstrings, i), env);
  array_delete(sa->cstrings, env);
  env_ma_free(sa, env);
}
