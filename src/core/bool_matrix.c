/*
  Copyright (c) 2005-2007, 2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007       Center for Bioinformatics, University of Hamburg

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

#include "core/array.h"
#include "core/bool_matrix.h"
#include "core/dyn_bittab.h"
#include "core/ma.h"
#include "core/undef_api.h"

struct GtBoolMatrix {
  GtArray *dyn_bittabs;
};

GtBoolMatrix* gt_bool_matrix_new(void)
{
  GtBoolMatrix *bm;
  bm = gt_malloc(sizeof (GtBoolMatrix));
  bm->dyn_bittabs = gt_array_new(sizeof (GtDynBittab*));
  return bm;
}

bool gt_bool_matrix_get(GtBoolMatrix *bm, unsigned long firstdim,
                        unsigned long seconddim)
{
  GtDynBittab *bt;
  gt_assert(bm);
  if (firstdim < gt_array_size(bm->dyn_bittabs)) {
    bt = *(GtDynBittab**) gt_array_get(bm->dyn_bittabs, firstdim);
    if (bt && gt_dyn_bittab_bit_is_set(bt, seconddim))
      return true;
  }
  return false;
}

void gt_bool_matrix_set(GtBoolMatrix *bm, unsigned long firstdim,
                        unsigned long seconddim, bool b)
{
  GtDynBittab *bt;
  unsigned long i, elems_to_add;

  gt_assert(bm);
  /* make sure first dimension is large enough */
  if (firstdim >= gt_array_size(bm->dyn_bittabs)) {
    elems_to_add = firstdim - gt_array_size(bm->dyn_bittabs) + 1;
    for (i = 0; i < elems_to_add; i++) {
      bt = gt_dyn_bittab_new();
      gt_array_add(bm->dyn_bittabs, bt);
    }
  }

  /* set bit, if necessary */
  bt = *(GtDynBittab**) gt_array_get(bm->dyn_bittabs, firstdim);
  gt_assert(bt);
  if (b)
    gt_dyn_bittab_set_bit(bt, seconddim);
  else
    gt_dyn_bittab_unset_bit(bt, seconddim);

  /* value has been set */
  gt_assert(gt_bool_matrix_get(bm, firstdim, seconddim) == b);
}

unsigned long gt_bool_matrix_get_first_column(const GtBoolMatrix *bm,
                                              unsigned long firstdim)
{
  GtDynBittab *bt;
  gt_assert(bm);
  if (firstdim < gt_array_size(bm->dyn_bittabs)) {
    if ((bt = *(GtDynBittab**) gt_array_get(bm->dyn_bittabs, firstdim)))
      return gt_dyn_bittab_get_first_bitnum(bt);
  }
  return GT_UNDEF_ULONG;
}

unsigned long gt_bool_matrix_get_last_column(const GtBoolMatrix *bm,
                                             unsigned long firstdim)
{
  GtDynBittab *bt;
  gt_assert(bm);
  if (firstdim < gt_array_size(bm->dyn_bittabs)) {
    if ((bt = *(GtDynBittab**) gt_array_get(bm->dyn_bittabs, firstdim)))
      return gt_dyn_bittab_get_last_bitnum(bt);
  }
  return GT_UNDEF_ULONG;
}

unsigned long gt_bool_matrix_get_next_column(const GtBoolMatrix *bm,
                                             unsigned long firstdim,
                                             unsigned long i)
{
  GtDynBittab *bt;
  gt_assert(bm);
  if (firstdim < gt_array_size(bm->dyn_bittabs)) {
    if ((bt = *(GtDynBittab**) gt_array_get(bm->dyn_bittabs, firstdim)))
      return gt_dyn_bittab_get_next_bitnum(bt, i);
  }
  return GT_UNDEF_ULONG;
}

void gt_bool_matrix_delete(GtBoolMatrix *bm)
{
  unsigned long i;
  if (!bm) return;
  for (i = 0; i < gt_array_size(bm->dyn_bittabs); i++)
    gt_dyn_bittab_delete(*(GtDynBittab**) gt_array_get(bm->dyn_bittabs, i));
  gt_array_delete(bm->dyn_bittabs);
  gt_free(bm);
}
