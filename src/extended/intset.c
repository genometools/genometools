/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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
#include "extended/intset_rep.h"

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/types_api.h"
#include "core/unused_api.h"

GtIntset *gt_intset_create(const GtIntsetClass *intset_c)
{
  GtIntset *intset;
  gt_assert(intset_c && intset_c->size);
  intset = gt_calloc((size_t) 1, intset_c->size);
  intset->c_class = intset_c;
  intset->members = gt_calloc((size_t) 1, sizeof (GtIntsetMembers));
  return intset;
}

GtIntset *gt_intset_ref(GtIntset *intset)
{
  gt_assert(intset);
  intset->members->refcount++;
  return intset;
}

void *gt_intset_cast(GT_UNUSED const GtIntsetClass *intset_c,
                     GtIntset *intset)
{
  gt_assert(intset_c && intset &&
            intset->c_class == intset_c);
  return intset;
}

void gt_intset_add(GtIntset *intset, GtUword elem)
{
  gt_assert(intset != NULL);
  gt_assert(intset->c_class != NULL);
  if (intset->c_class->add_func != NULL)
    intset->c_class->add_func(intset, elem);
}

GtUword gt_intset_get(GtIntset *intset, GtUword idx)
{
  gt_assert(intset != NULL);
  gt_assert(intset->c_class != NULL);
  if (intset->c_class->get_func != NULL)
    return intset->c_class->get_func(intset, idx);
  return intset->members->num_of_elems;
}

bool gt_intset_is_member(GtIntset *intset, GtUword elem)
{
  gt_assert(intset != NULL);
  gt_assert(intset->c_class != NULL);
  if (intset->c_class->is_member_func != NULL)
    return intset->c_class->is_member_func(intset, elem);
  return false;
}

GtUword gt_intset_get_idx_smaller_geq(GtIntset *intset, GtUword pos)
{
  gt_assert(intset != NULL);
  gt_assert(intset->c_class != NULL);
  if (intset->c_class->idx_sm_geq_func != NULL)
    return intset->c_class->idx_sm_geq_func(intset, pos);
  return GT_UWORD_MAX;
}

const GtIntsetClass *gt_intset_class_new(size_t size,
                                         GtIntsetAddFunc add_func,
                                         GtIntsetGetFunc get_func,
                                         GtIntsetIdxSmGeqFunc idx_sm_geq_func,
                                         GtIntsetIsMemberFunc is_member_func,
                                         GtIntsetDeleteFunc delete_func)
{
  GtIntsetClass *intset_c = gt_class_alloc(sizeof (*intset_c));
  intset_c->size = size;
  intset_c->add_func = add_func;
  intset_c->get_func = get_func;
  intset_c->is_member_func = is_member_func;
  intset_c->idx_sm_geq_func = idx_sm_geq_func;
  intset_c->delete_func = delete_func;
  return intset_c;
}

void gt_intset_delete(GtIntset *intset)
{
  if (intset != NULL) {
    if (intset->members->refcount) {
      intset->members->refcount--;
      return;
    }
    gt_assert(intset->c_class);
    if (intset->c_class->delete_func != NULL)
      intset->c_class->delete_func(intset);
    gt_free(intset->members->sectionstart);
    gt_free(intset->members);
    gt_free(intset);
  }
}

int gt_intset_unit_test_notinset(GtIntset *intset, GtUword start,
                                 GtUword end, GtError *err)
{
  int had_err = 0;
  GtUword test;
  for (test = start; test <= end; ++test) {
    gt_ensure(gt_intset_is_member(intset, test) == false);
  }
  return had_err;
}

int gt_intset_unit_test_check_seqnum(GtIntset *intset, GtUword start,
                                     GtUword end, GtUword num, GtError *err)
{
  int had_err = 0;
  GtUword test;
  for (test = start; test <= end; ++test) {
    gt_ensure(gt_intset_get_idx_smaller_geq(intset, test) == num);
  }
  return had_err;
}
