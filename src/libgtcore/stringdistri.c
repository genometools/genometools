/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/stringdistri.h"

struct StringDistri {
  Hashtable *hashdist;
  unsigned long num_of_occurrences;
};

StringDistri* stringdistri_new(Env *env)
{
  StringDistri *sd;
  env_error_check(env);
  sd = ma_malloc(sizeof *sd);
  sd->hashdist = hashtable_new(HASH_STRING, env_ma_free_func, env_ma_free_func,
                               env);
  sd->num_of_occurrences = 0;
  return sd;
}

void stringdistri_add(StringDistri *d, const char *key, Env *env)
{
  unsigned long *valueptr;
  assert(d && key);
  valueptr = hashtable_get(d->hashdist, (void*) key);
  if (!valueptr) {
    valueptr = ma_malloc(sizeof *valueptr);
    *valueptr = 1;
    hashtable_add(d->hashdist, cstr_dup(key, env), valueptr, env);
  }
  else
    (*valueptr)++;
  d->num_of_occurrences++;
}

typedef struct {
  StringDistriIterFunc func;
  void *data;
  unsigned long num_of_occurrences;
} ForeachInfo;

static int foreach_iterfunc(void *key, void *value, void *data, Env *env)
{
  unsigned long occurrences;
  ForeachInfo *info;
  env_error_check(env);
  assert(key && value && data);
  occurrences = *(unsigned long*) value;
  info = (ForeachInfo*) data;
  info->func(key, occurrences, (double) occurrences / info->num_of_occurrences,
             info->data);
  return 0;
}

void stringdistri_foreach(const StringDistri *d, StringDistriIterFunc func,
                        void *data, Env *env)
{
  ForeachInfo info;
  int rval;
  env_error_check(env);
  assert(d);
  if (d->hashdist) {
    info.func = func;
    info.data = data;
    info.num_of_occurrences = d->num_of_occurrences;
    rval = hashtable_foreach_ao(d->hashdist, foreach_iterfunc, &info, env);
    assert(!rval); /* foreach_iterfunc() is sane */
  }
}

void stringdistri_delete(StringDistri *d, Env *env)
{
  if (!d) return;
  hashtable_delete(d->hashdist, env);
  ma_free(d);
}
