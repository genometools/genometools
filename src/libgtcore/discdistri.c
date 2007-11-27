/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include <stdio.h>
#include "libgtcore/discdistri.h"
#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"

struct DiscDistri {
  Hashtable *hashdist;
  unsigned long long num_of_occurrences;
};

DiscDistri* discdistri_new(Env *env)
{
  return ma_calloc(1, sizeof (DiscDistri));
}

void discdistri_add(DiscDistri *d, unsigned long key, Env *env)
{
  discdistri_add_multi(d, key, 1, env);
}

void discdistri_add_multi(DiscDistri *d, unsigned long key,
                          unsigned long long occurrences, Env *env)
{
  unsigned long long *valueptr;
  assert(d);

  if (!d->hashdist)
    d->hashdist = hashtable_new(HASH_DIRECT, NULL, ma_free_func, env);

  valueptr = hashtable_get(d->hashdist, (void*) key);
  if (!valueptr) {
    valueptr = ma_malloc(sizeof *valueptr);
    *valueptr = occurrences;
    hashtable_add(d->hashdist, (void*) key, valueptr, env);
  }
  else
    (*valueptr) += occurrences;

  d->num_of_occurrences += occurrences;
}

unsigned long long discdistri_get(const DiscDistri *d, unsigned long key)
{
  unsigned long long *valueptr;
  assert(d);
  if (!d->hashdist || !(valueptr = hashtable_get(d->hashdist, (void*) key)))
    return 0;
  return *valueptr;
}

void discdistri_show(const DiscDistri *d, Env *env)
{
  assert(d);
  discdistri_show_generic(d, NULL, env);
}

typedef struct {
  double cumulative_probability;
  unsigned long long num_of_occurrences;
  GenFile *genfile;
} ShowValueInfo;

static int showvalue(void *key, void *value, void *data, Env *env)
{
  unsigned long long occurrences;
  double probability;
  ShowValueInfo *info;

  env_error_check(env);
  assert(key && value && data);

  occurrences = *(unsigned long long*) value;
  assert(occurrences);
  info = (ShowValueInfo*) data;

  probability = (double) occurrences / info->num_of_occurrences;
  info->cumulative_probability += probability;
  genfile_xprintf(info->genfile, "%lu: %llu (prob=%.4f,cumulative=%.4f)\n",
                  (unsigned long) key, occurrences, probability,
                  info->cumulative_probability);

  return 0;
}

void discdistri_show_generic(const DiscDistri *d, GenFile *genfile, Env *env)
{
  ShowValueInfo showvalueinfo;
  int rval;

  env_error_check(env);
  assert(d);

  if (d->hashdist) {
    showvalueinfo.cumulative_probability = 0.0;
    showvalueinfo.num_of_occurrences = d->num_of_occurrences;
    showvalueinfo.genfile = genfile;
    rval = hashtable_foreach_no(d->hashdist, showvalue, &showvalueinfo, env);
    assert(!rval); /* showvalue() is sane */
  }
}

typedef struct {
  DiscDistriIterFunc func;
  void *data;
} ForeachInfo;

static int foreach_iterfunc(void *key, void *value, void *data, Env *env)
{
  ForeachInfo *info;
  env_error_check(env);
  assert(key && value && data);
  info = (ForeachInfo*) data;
  info->func((unsigned long) key, *(unsigned long long*) value, info->data);
  return 0;
}

void discdistri_foreach(const DiscDistri *d, DiscDistriIterFunc func,
                        void *data, Env *env)
{
  ForeachInfo info;
  int rval;
  env_error_check(env);
  assert(d);
  if (d->hashdist) {
    info.func = func;
    info.data = data;
    rval = hashtable_foreach_no(d->hashdist, foreach_iterfunc, &info, env);
    assert(!rval); /* foreach_iterfunc() is sane */
  }
}

int discdistri_unit_test(Env *env)
{
  DiscDistri *d;
  int had_err = 0;

  env_error_check(env);

  d = discdistri_new(env);

  ensure(had_err, discdistri_get(d, 0) == 0);
  ensure(had_err, discdistri_get(d, 100) == 0);
  if (!had_err) {
    discdistri_add(d, 0, env);
    discdistri_add_multi(d, 100, 256, env);
  }
  ensure(had_err, discdistri_get(d, 0) == 1);
  ensure(had_err, discdistri_get(d, 100) == 256);

  discdistri_delete(d, env);

  return had_err;
}

void discdistri_delete(DiscDistri *d, Env *env)
{
  if (!d) return;
  hashtable_delete(d->hashdist, env);
  ma_free(d);
}
