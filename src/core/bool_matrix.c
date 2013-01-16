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

struct GtBoolMatrix {
  GtArray *dyn_bittabs;
};

GtBoolMatrix* gt_bool_matrix_new(void)
{
  GtBoolMatrix *boolmatrix;
  boolmatrix = gt_malloc(sizeof (GtBoolMatrix));
  boolmatrix->dyn_bittabs = gt_array_new(sizeof (GtDynBittab*));
  return boolmatrix;
}

bool gt_bool_matrix_get(GtBoolMatrix *boolmatrix, unsigned long firstdim,
                       unsigned long seconddim)
{
  GtDynBittab *bt;
  gt_assert(boolmatrix);
  if (firstdim < gt_array_size(boolmatrix->dyn_bittabs)) {
    bt = *(GtDynBittab**) gt_array_get(boolmatrix->dyn_bittabs, firstdim);
    if (bt && gt_dyn_bittab_bit_is_set(bt, seconddim))
      return true;
  }
  return false;
}

void gt_bool_matrix_set(GtBoolMatrix *boolmatrix, unsigned long firstdim,
                       unsigned long seconddim, bool b)
{
  GtDynBittab *bt;
  unsigned long i, elems_to_add;

  gt_assert(boolmatrix);
  /* make sure first dimension is large enough */
  if (firstdim <= gt_array_size(boolmatrix->dyn_bittabs)) {
    elems_to_add = firstdim - gt_array_size(boolmatrix->dyn_bittabs) + 1;
    for (i = 0; i < elems_to_add; i++) {
      bt = gt_dyn_bittab_new();
      gt_array_add(boolmatrix->dyn_bittabs, bt);
    }
  }

  /* set bit, if necessary */
  bt = *(GtDynBittab**) gt_array_get(boolmatrix->dyn_bittabs, firstdim);
  gt_assert(bt);
  if (b)
    gt_dyn_bittab_set_bit(bt, seconddim);
  else
    gt_dyn_bittab_unset_bit(bt, seconddim);

  /* value has been set */
  gt_assert(gt_bool_matrix_get(boolmatrix, firstdim, seconddim) == b);
}

void gt_bool_matrix_delete(GtBoolMatrix *boolmatrix)
{
  unsigned long i;
  if (!boolmatrix) return;
  for (i = 0; i < gt_array_size(boolmatrix->dyn_bittabs); i++) {
    gt_dyn_bittab_delete(*(GtDynBittab**)
                        gt_array_get(boolmatrix->dyn_bittabs, i));
  }
  gt_array_delete(boolmatrix->dyn_bittabs);
  gt_free(boolmatrix);
}
