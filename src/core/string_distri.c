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

#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/hashmap-generic.h"
#include "core/ma.h"
#include "core/string_distri.h"
#include "core/unused_api.h"

struct GtStringDistri {
  GtHashtable *hashdist;
  unsigned long num_of_occurrences;
};

DECLARE_HASHMAP(char *, cstr, unsigned long, ul, static, inline)
DEFINE_HASHMAP(char *, cstr, unsigned long, ul, gt_ht_cstr_elem_hash,
               gt_ht_cstr_elem_cmp, gt_free, NULL_DESTRUCTOR, static,
               inline)

GtStringDistri* gt_string_distri_new(void)
{
  GtStringDistri *sd;
  sd = gt_malloc(sizeof *sd);
  sd->hashdist = cstr_ul_gt_hashmap_new();
  sd->num_of_occurrences = 0;
  return sd;
}

void gt_string_distri_add(GtStringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  gt_assert(sd && key);
  valueptr = cstr_ul_gt_hashmap_get(sd->hashdist, key);
  if (!valueptr) {
    cstr_ul_gt_hashmap_add(sd->hashdist, gt_cstr_dup(key), 1);
  }
  else
    (*valueptr)++;
  sd->num_of_occurrences++;
}

void gt_string_distri_sub(GtStringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  gt_assert(sd && key && gt_string_distri_get(sd, key) &&
            sd->num_of_occurrences);
  valueptr = cstr_ul_gt_hashmap_get(sd->hashdist, key);
  (*valueptr)--;
  if (!(*valueptr))
    cstr_ul_gt_hashmap_remove(sd->hashdist, key);
  sd->num_of_occurrences--;
}

unsigned long gt_string_distri_get(const GtStringDistri *sd, const char *key)
{
  unsigned long *valueptr;
  gt_assert(sd && key);
  if ((valueptr = cstr_ul_gt_hashmap_get(sd->hashdist, key)))
    return *valueptr;
  else
    return 0;
}

double gt_string_distri_get_prob(const GtStringDistri *sd, const char *key)
{
  unsigned long occ;
  gt_assert(sd && key);
  if ((occ = gt_string_distri_get(sd, key)))
    return (double) occ / sd->num_of_occurrences;
  return 0.0;
}

typedef struct {
  GtStringDistriIterFunc func;
  void *data;
  unsigned long num_of_occurrences;
} StringDistriForeachInfo;

static enum iterator_op
string_distri_foreach_iterfunc(char *key, unsigned long occurrences, void *data,
                               GT_UNUSED GtError *err)
{
  StringDistriForeachInfo *info;
  gt_error_check(err);
  gt_assert(key && data);
  info = (StringDistriForeachInfo*) data;
  info->func(key, occurrences, (double) occurrences / info->num_of_occurrences,
             info->data);
  return 0;
}

void gt_string_distri_foreach(const GtStringDistri *sd,
                              GtStringDistriIterFunc func, void *data)
{
  StringDistriForeachInfo info;
  GT_UNUSED int rval;
  gt_assert(sd);
  if (sd->hashdist) {
    info.func = func;
    info.data = data;
    info.num_of_occurrences = sd->num_of_occurrences;
    rval = cstr_ul_gt_hashmap_foreach_in_default_order(
      sd->hashdist, string_distri_foreach_iterfunc, &info, NULL);
    gt_assert(!rval); /* string_distri_foreach_iterfunc() is sane */
  }
}

void gt_string_distri_delete(GtStringDistri *sd)
{
  if (!sd) return;
  gt_hashtable_delete(sd->hashdist);
  gt_free(sd);
}
