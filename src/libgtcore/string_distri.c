/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/string_distri.h"
#include "libgtcore/unused.h"

struct StringDistri {
  Hashtable *hashdist;
  unsigned long num_of_occurrences;
};

StringDistri* string_distri_new(void)
{
  StringDistri *sd;
  sd = ma_malloc(sizeof *sd);
  sd->hashdist = hashtable_new(HASH_STRING, ma_free_func, ma_free_func);
  sd->num_of_occurrences = 0;
  return sd;
}

void string_distri_add(StringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  assert(sd && key);
  valueptr = hashtable_get(sd->hashdist, (void*) key);
  if (!valueptr) {
    valueptr = ma_malloc(sizeof *valueptr);
    *valueptr = 1;
    hashtable_add(sd->hashdist, cstr_dup(key), valueptr);
  }
  else
    (*valueptr)++;
  sd->num_of_occurrences++;
}

void string_distri_sub(StringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  assert(sd && key && string_distri_get(sd, key) && sd->num_of_occurrences);
  valueptr = hashtable_get(sd->hashdist, (void*) key);
  (*valueptr)--;
  if (!(*valueptr))
    hashtable_remove(sd->hashdist, key);
  sd->num_of_occurrences--;
}

unsigned long string_distri_get(const StringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  assert(sd && key);
  if ((valueptr = hashtable_get(sd->hashdist, (void*) key)))
    return *valueptr;
  else
    return 0;
}

typedef struct {
  StringDistriIterFunc func;
  void *data;
  unsigned long num_of_occurrences;
} ForeachInfo;

static int foreach_iterfunc(void *key, void *value, void *data,
                            UNUSED Error *err)
{
  unsigned long occurrences;
  ForeachInfo *info;
  error_check(err);
  assert(key && value && data);
  occurrences = *(unsigned long*) value;
  info = (ForeachInfo*) data;
  info->func(key, occurrences, (double) occurrences / info->num_of_occurrences,
             info->data);
  return 0;
}

void string_distri_foreach(const StringDistri *sd, StringDistriIterFunc func,
                          void *data)
{
  ForeachInfo info;
  int rval;
  assert(sd);
  if (sd->hashdist) {
    info.func = func;
    info.data = data;
    info.num_of_occurrences = sd->num_of_occurrences;
    rval = hashtable_foreach_ao(sd->hashdist, foreach_iterfunc, &info, NULL);
    assert(!rval); /* foreach_iterfunc() is sane */
  }
}

void string_distri_delete(StringDistri *sd)
{
  if (!sd) return;
  hashtable_delete(sd->hashdist);
  ma_free(sd);
}
