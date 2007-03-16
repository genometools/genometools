/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <libgtcore/cstr.h>
#include <libgtcore/xansi.h>

void* cstr_dup(const char *cstr, Env *env)
{
  size_t size;
  char *copy;
  assert(cstr);
  size = strlen(cstr) + 1;
  copy = env_ma_malloc(env, size);
  memcpy(copy, cstr, size);
  return copy;
}

void cstr_show(const char *cstr, unsigned long length, FILE *fp)
{
  unsigned long i;
  assert(cstr && fp);
  for (i = 0; i < length; i++)
    xfputc(cstr[i], fp);
}

char** cstr_array_prefix_first(const char **cstr_array, const char *p, Env *env)
{
  unsigned long i, a_len, f_len;
  char **a;
  assert(cstr_array && p);
  a_len = cstr_array_size(cstr_array);
  a = env_ma_malloc(env, sizeof (char*) * (a_len + 1));
  f_len = strlen(p) + strlen(cstr_array[0]) + 2; /* blank + '\0' */
  a[0] = env_ma_malloc(env, sizeof (char) * f_len);
  (void) snprintf(a[0], f_len, "%s %s", p, cstr_array[0]);
  for (i = 1; i < a_len; i++)
    a[i] = cstr_dup(cstr_array[i], env);
  a[a_len] = NULL;
  return a;
}

char** cstr_array_preprend(const char **cstr_array, const char *p, Env *env)
{
  unsigned long i, a_len;
  char **a;
  assert(cstr_array && p);
  a_len = cstr_array_size(cstr_array);
  a = env_ma_malloc(env, sizeof (char*) * (a_len + 2));
  a[0] = cstr_dup(p, env);
  for (i = 0; i < a_len; i++)
    a[i+1] = cstr_dup(cstr_array[i], env);
  a[a_len+1] = NULL;
  return a;
}

void cstr_array_show(char **cstr_array, FILE *fp)
{
  unsigned long i = 0;
  while (cstr_array[i]) {
    if (i)
      xfputc(' ', fp);
    xfputs(cstr_array[i], fp);
    i++;
  }
  xfputc('\n', fp);
}

unsigned long cstr_array_size(const char **cstr_array)
{
  unsigned long i = 0;
  while (cstr_array[i])
    i++;
  return i;
}

void cstr_array_delete(char **cstr_array, Env *env)
{
  unsigned long i = 0;
  if (!cstr_array) return;
  while (cstr_array[i])
    env_ma_free(cstr_array[i++], env);
  env_ma_free(cstr_array, env);
}
