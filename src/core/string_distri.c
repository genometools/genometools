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
#include "core/cstr.h"
#include "core/hashmap-generic.h"
#include "core/ma.h"
#include "core/string_distri.h"
#include "core/unused_api.h"

struct StringDistri {
  Hashtable *hashdist;
  unsigned long num_of_occurrences;
};

DECLARE_HASHMAP(char *, cstr, unsigned long, ul, static, inline)
DEFINE_HASHMAP(char *, cstr, unsigned long, ul, ht_cstr_elem_hash,
               ht_cstr_elem_cmp, gt_free, NULL_DESTRUCTOR, static,
               inline)

StringDistri* string_distri_new(void)
{
  StringDistri *sd;
  sd = gt_malloc(sizeof *sd);
  sd->hashdist = cstr_ul_hashmap_new();
  sd->num_of_occurrences = 0;
  return sd;
}

void string_distri_add(StringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  assert(sd && key);
  valueptr = cstr_ul_hashmap_get(sd->hashdist, key);
  if (!valueptr) {
    cstr_ul_hashmap_add(sd->hashdist, gt_cstr_dup(key), 1);
  }
  else
    (*valueptr)++;
  sd->num_of_occurrences++;
}

void string_distri_sub(StringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  assert(sd && key && string_distri_get(sd, key) && sd->num_of_occurrences);
  valueptr = cstr_ul_hashmap_get(sd->hashdist, key);
  (*valueptr)--;
  if (!(*valueptr))
    cstr_ul_hashmap_remove(sd->hashdist, key);
  sd->num_of_occurrences--;
}

unsigned long string_distri_get(const StringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  assert(sd && key);
  if ((valueptr = cstr_ul_hashmap_get(sd->hashdist, key)))
    return *valueptr;
  else
    return 0;
}

typedef struct {
  StringDistriIterFunc func;
  void *data;
  unsigned long num_of_occurrences;
} ForeachInfo;

static enum iterator_op
foreach_iterfunc(char *key, unsigned long occurrences, void *data,
                 GT_UNUSED GtError *err)
{
  ForeachInfo *info;
  gt_error_check(err);
  assert(key && data);
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
    rval = cstr_ul_hashmap_foreach_in_default_order(
      sd->hashdist, foreach_iterfunc, &info, NULL);
    assert(!rval); /* foreach_iterfunc() is sane */
  }
}

void string_distri_delete(StringDistri *sd)
{
  if (!sd) return;
  hashtable_delete(sd->hashdist);
  gt_free(sd);
}
