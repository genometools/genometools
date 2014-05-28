/*
  Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2008 Thomas Jahns <Thomas.Jahns@gmx.net>
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

#include <stdio.h>

#include "core/assert_api.h"
#include "core/compat.h"
#include "core/ensure.h"
#include "core/hashmap-generic.h"
#include "core/ma.h"
#include "core/unused_api.h"

#include "core/disc_distri_api.h"
#include "core/disc_distri.h"

struct GtDiscDistri {
  GtHashtable *hashdist;
  GtUint64 num_of_occurrences;
};

GtDiscDistri* gt_disc_distri_new(void)
{
  return gt_calloc((size_t) 1, sizeof (GtDiscDistri));
}

void gt_disc_distri_add(GtDiscDistri *d, GtUword key)
{
  gt_disc_distri_add_multi(d, key, (GtUint64) 1);
}

DECLARE_HASHMAP(GtUword, ul, GtUint64, ull, static, inline)
DEFINE_HASHMAP(GtUword, ul, GtUint64, ull, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

void gt_disc_distri_add_multi(GtDiscDistri *d, GtUword key,
                              GtUint64 occurrences)
{
  GtUint64 *valueptr;
  gt_assert(d);

  if (!d->hashdist)
    d->hashdist = ul_ull_gt_hashmap_new();

  valueptr = ul_ull_gt_hashmap_get(d->hashdist, key);
  if (!valueptr) {
    ul_ull_gt_hashmap_add(d->hashdist, key, occurrences);
  }
  else
    (*valueptr) += occurrences;

  d->num_of_occurrences += occurrences;
}

GtUint64 gt_disc_distri_get(const GtDiscDistri *d, GtUword key)
{
  GtUint64 *valueptr;
  gt_assert(d);
  if (!d->hashdist || !(valueptr = ul_ull_gt_hashmap_get(d->hashdist, key)))
    return 0;
  return *valueptr;
}

typedef struct {
  double cumulative_probability;
  GtUint64 num_of_occurrences;
  GtFile *outfp;
} GtShowValueInfo;

static enum iterator_op
showvalue(GtUword key, GtUint64 occurrences,
          void *data, GT_UNUSED GtError *err)
{
  double probability;
  GtShowValueInfo *info;

  gt_error_check(err);
  gt_assert(data && occurrences);
  info = (GtShowValueInfo*) data;

  probability = (double) ((double) occurrences / info->num_of_occurrences);
  info->cumulative_probability += probability;
  gt_file_xprintf(info->outfp, GT_WU": "GT_LLU" (prob=%.4f,cumulative=%.4f)\n",
                  key, occurrences, probability, info->cumulative_probability);
  return CONTINUE_ITERATION;
}

void gt_disc_distri_show(const GtDiscDistri *d, GtFile *outfp)
{
  GtShowValueInfo showvalueinfo;
  GT_UNUSED int rval;

  gt_assert(d);

  if (d->hashdist != NULL) {
    showvalueinfo.cumulative_probability = 0.0;
    showvalueinfo.num_of_occurrences = d->num_of_occurrences;
    showvalueinfo.outfp = outfp;
    rval = ul_ull_gt_hashmap_foreach_in_default_order(d->hashdist, showvalue,
                                                      &showvalueinfo, NULL);
    gt_assert(!rval); /* showvalue() is sane */
  }
}

typedef struct {
  GtDiscDistriIterFunc func;
  void *data;
} DiscDistriForeachInfo;

static enum iterator_op
disc_distri_foreach_iterfunc(GtUword key, GtUint64 occurrences,
                             void *data, GT_UNUSED GtError *err)
{
  DiscDistriForeachInfo *info;
  gt_error_check(err);
  gt_assert(data);
  info = (DiscDistriForeachInfo*) data;
  info->func(key, occurrences, info->data);
  return CONTINUE_ITERATION;
}

static
void gt_disc_distri_foreach_generic(const GtDiscDistri *d,
                                    GtDiscDistriIterFunc func,
                                    void *data,
                                    ul_ull_gt_hashmap_KeyCmp cmp)
{
  DiscDistriForeachInfo info;
  GT_UNUSED int rval;
  gt_assert(d);
  if (d->hashdist != NULL) {
    info.func = func;
    info.data = data;
    if (cmp != NULL)
      rval = ul_ull_gt_hashmap_foreach_ordered(d->hashdist,
                                               disc_distri_foreach_iterfunc,
                                               &info, cmp, NULL);
    else
      rval = ul_ull_gt_hashmap_foreach_in_default_order(d->hashdist,
                                               disc_distri_foreach_iterfunc,
                                               &info, NULL);
    gt_assert(!rval); /* disc_distri_foreach_iterfunc() is sane */
  }
}

void gt_disc_distri_foreach(const GtDiscDistri *d, GtDiscDistriIterFunc func,
                            void *data)
{
  gt_disc_distri_foreach_generic(d,func,data,NULL);
}

static int
rev_key_cmp(const GtUword a, const GtUword b)
{
  return -gt_ht_ul_elem_cmp(&a,&b);
}

void gt_disc_distri_foreach_in_reverse_order(const GtDiscDistri *d,
                                             GtDiscDistriIterFunc func,
                                             void *data)
{
  gt_disc_distri_foreach_generic(d,func,data, rev_key_cmp);
}

#define DISC_DISTRI_FOREACHTESTSIZE 3

/* data for foreach unit test */
struct ForeachTesterData
{
  int counter;
  int expkeys[DISC_DISTRI_FOREACHTESTSIZE];
  int expvalues[DISC_DISTRI_FOREACHTESTSIZE];
  int *had_err;
  GtError *err;
};

/* helper function for unit test of foreach */
static void foreachtester(GtUword key,
                          GtUint64 value, void *data)
{
  struct ForeachTesterData *tdata = data;
  int had_err = *(tdata->had_err);
  GtError *err = tdata->err;
  tdata->counter++;
  gt_ensure(tdata->counter < DISC_DISTRI_FOREACHTESTSIZE);
  gt_ensure((GtUword) tdata->expkeys[tdata->counter] == key);
  gt_ensure((GtUint64) tdata->expvalues[tdata->counter] == value);
  *(tdata->had_err) = had_err;
}

int gt_disc_distri_unit_test(GtError *err)
{
  GtDiscDistri *d;
  int had_err = 0;
  struct ForeachTesterData tdata;

  gt_error_check(err);

  d = gt_disc_distri_new();

  gt_ensure(gt_disc_distri_get(d, 0UL) == 0);
  gt_ensure(gt_disc_distri_get(d, 100UL) == 0);
  if (!had_err) {
    gt_disc_distri_add(d, 0);
    gt_disc_distri_add_multi(d, 100UL, 256ULL);
  }
  gt_ensure(gt_disc_distri_get(d, 0UL) == 1ULL);
  gt_ensure(gt_disc_distri_get(d, 100UL) == 256ULL);

  /* test foreach and foreach_in_reverse_order: */
  gt_disc_distri_add(d, 2UL);
  if (!had_err) {
    tdata.counter = -1;
    tdata.expkeys[0] = 0;
    tdata.expvalues[0] = 1;
    tdata.expkeys[1] = 2;
    tdata.expvalues[1] = 1;
    tdata.expkeys[2] = 100;
    tdata.expvalues[2] = 256;
    tdata.had_err = &had_err;
    tdata.err = err;
    gt_disc_distri_foreach(d, foreachtester, &tdata);
  }
  if (!had_err) {
    tdata.counter = -1;
    tdata.expkeys[0] = 100;
    tdata.expvalues[0] = 256;
    tdata.expkeys[1] = 2;
    tdata.expvalues[1] = 1;
    tdata.expkeys[2] = 0;
    tdata.expvalues[2] = 1;
    tdata.had_err = &had_err;
    tdata.err = err;
    gt_disc_distri_foreach_in_reverse_order(d, foreachtester, &tdata);
  }

  gt_disc_distri_delete(d);

  return had_err;
}

void gt_disc_distri_delete(GtDiscDistri *d)
{
  if (!d) return;
  gt_hashtable_delete(d->hashdist);
  gt_free(d);
}
