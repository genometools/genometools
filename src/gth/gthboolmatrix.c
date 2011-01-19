/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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
#include "core/dynbittab.h"
#include "core/ma.h"
#include "gth/gthboolmatrix.h"

struct GthBoolMatrix {
  GtArray *dynbittabs;
};

GthBoolMatrix* gthboolmatrix_new(void)
{
  GthBoolMatrix *boolmatrix;
  boolmatrix = gt_malloc(sizeof (GthBoolMatrix));
  boolmatrix->dynbittabs = gt_array_new(sizeof (GtDynBittab*));
  return boolmatrix;
}

bool gthboolmatrix_get(GthBoolMatrix *boolmatrix, unsigned long firstdim,
                       unsigned long seconddim)
{
  GtDynBittab *bt;
  gt_assert(boolmatrix);
  if (firstdim < gt_array_size(boolmatrix->dynbittabs)) {
    bt = *(GtDynBittab**) gt_array_get(boolmatrix->dynbittabs, firstdim);
    if (bt && gt_dynbittab_bit_is_set(bt, seconddim))
      return true;
  }
  return false;
}

void gthboolmatrix_set(GthBoolMatrix *boolmatrix, unsigned long firstdim,
                       unsigned long seconddim, bool b)
{
  GtDynBittab *bt;
  unsigned long i, elems_to_add;

  gt_assert(boolmatrix);
  /* make sure first dimension is large enough */
  if (firstdim <= gt_array_size(boolmatrix->dynbittabs)) {
    elems_to_add = firstdim - gt_array_size(boolmatrix->dynbittabs) + 1;
    for (i = 0; i < elems_to_add; i++) {
      bt = gt_dynbittab_new();
      gt_array_add(boolmatrix->dynbittabs, bt);
    }
  }

  /* set bit, if necessary */
  bt = *(GtDynBittab**) gt_array_get(boolmatrix->dynbittabs, firstdim);
  gt_assert(bt);
  if (b)
    gt_dynbittab_set_bit(bt, seconddim);
  else
    gt_dynbittab_unset_bit(bt, seconddim);

  /* value has been set */
  gt_assert(gthboolmatrix_get(boolmatrix, firstdim, seconddim) == b);
}

void gthboolmatrix_delete(GthBoolMatrix *boolmatrix)
{
  unsigned long i;
  if (!boolmatrix) return;
  for (i = 0; i < gt_array_size(boolmatrix->dynbittabs); i++) {
    gt_dynbittab_delete(*(GtDynBittab**)
                        gt_array_get(boolmatrix->dynbittabs, i));
  }
  gt_array_delete(boolmatrix->dynbittabs);
  gt_free(boolmatrix);
}
