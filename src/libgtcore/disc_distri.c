/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/disc_distri.h"
#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"

struct DiscDistri {
  Hashtable *hashdist;
  unsigned long long num_of_occurrences;
};

DiscDistri* disc_distri_new(void)
{
  return ma_calloc(1, sizeof (DiscDistri));
}

void disc_distri_add(DiscDistri *d, unsigned long key)
{
  disc_distri_add_multi(d, key, 1);
}

void disc_distri_add_multi(DiscDistri *d, unsigned long key,
                          unsigned long long occurrences)
{
  unsigned long long *valueptr;
  assert(d);

  if (!d->hashdist)
    d->hashdist = hashtable_new(HASH_DIRECT, NULL, ma_free_func);

  valueptr = hashtable_get(d->hashdist, (void*) key);
  if (!valueptr) {
    valueptr = ma_malloc(sizeof *valueptr);
    *valueptr = occurrences;
    hashtable_add(d->hashdist, (void*) key, valueptr);
  }
  else
    (*valueptr) += occurrences;

  d->num_of_occurrences += occurrences;
}

unsigned long long disc_distri_get(const DiscDistri *d, unsigned long key)
{
  unsigned long long *valueptr;
  assert(d);
  if (!d->hashdist || !(valueptr = hashtable_get(d->hashdist, (void*) key)))
    return 0;
  return *valueptr;
}

void disc_distri_show(const DiscDistri *d)
{
  assert(d);
  disc_distri_show_generic(d, NULL);
}

typedef struct {
  double cumulative_probability;
  unsigned long long num_of_occurrences;
  GenFile *genfile;
} ShowValueInfo;

static int showvalue(void *key, void *value, void *data, UNUSED Error *err)
{
  unsigned long long occurrences;
  double probability;
  ShowValueInfo *info;

  error_check(err);
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

void disc_distri_show_generic(const DiscDistri *d, GenFile *genfile)
{
  ShowValueInfo showvalueinfo;
  int rval;

  assert(d);

  if (d->hashdist) {
    showvalueinfo.cumulative_probability = 0.0;
    showvalueinfo.num_of_occurrences = d->num_of_occurrences;
    showvalueinfo.genfile = genfile;
    rval = hashtable_foreach_no(d->hashdist, showvalue, &showvalueinfo, NULL);
    assert(!rval); /* showvalue() is sane */
  }
}

typedef struct {
  DiscDistriIterFunc func;
  void *data;
} ForeachInfo;

static int foreach_iterfunc(void *key, void *value, void *data,
                            UNUSED Error *err)
{
  ForeachInfo *info;
  error_check(err);
  assert(value && data);
  info = (ForeachInfo*) data;
  info->func((unsigned long) key, *(unsigned long long*) value, info->data);
  return 0;
}

void disc_distri_foreach(const DiscDistri *d, DiscDistriIterFunc func,
                        void *data)
{
  ForeachInfo info;
  int rval;
  assert(d);
  if (d->hashdist) {
    info.func = func;
    info.data = data;
    rval = hashtable_foreach_no(d->hashdist, foreach_iterfunc, &info, NULL);
    assert(!rval); /* foreach_iterfunc() is sane */
  }
}

int disc_distri_unit_test(Error *err)
{
  DiscDistri *d;
  int had_err = 0;

  error_check(err);

  d = disc_distri_new();

  ensure(had_err, disc_distri_get(d, 0) == 0);
  ensure(had_err, disc_distri_get(d, 100) == 0);
  if (!had_err) {
    disc_distri_add(d, 0);
    disc_distri_add_multi(d, 100, 256);
  }
  ensure(had_err, disc_distri_get(d, 0) == 1);
  ensure(had_err, disc_distri_get(d, 100) == 256);

  disc_distri_delete(d);

  return had_err;
}

void disc_distri_delete(DiscDistri *d)
{
  if (!d) return;
  hashtable_delete(d->hashdist);
  ma_free(d);
}
