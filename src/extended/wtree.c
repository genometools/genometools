/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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
#include "extended/wtree_rep.h"

#include "core/assert_api.h"
#include "core/ma.h"
#include "core/unused_api.h"

#include <limits.h>

GtWtree *gt_wtree_create(const GtWtreeClass *wtree_c)
{
  GtWtree *wtree;
  gt_assert(wtree_c && wtree_c->size);
  wtree = gt_calloc((size_t) 1, wtree_c->size);
  wtree->c_class = wtree_c;
  wtree->members = gt_calloc((size_t) 1, sizeof (GtWtreeMembers));
  return wtree;
}

GtWtree *gt_wtree_ref(GtWtree *wtree)
{
  gt_assert(wtree);
  wtree->members->refcount++;
  return wtree;
}

void *gt_wtree_cast(GT_UNUSED const GtWtreeClass *wtree_c,
                    GtWtree *wtree)
{
  gt_assert(wtree_c && wtree &&
            wtree->c_class == wtree_c);
  return wtree;
}

GtWtreeSymbol gt_wtree_access(GtWtree *wtree,
                              GtUword pos)
{
  gt_assert(wtree != NULL);
  gt_assert(wtree->c_class != NULL);
  if (wtree->c_class->access_func != NULL)
    return wtree->c_class->access_func(wtree, pos);
  return (GtWtreeSymbol) ULONG_MAX;
}

GtUword gt_wtree_rank(GtWtree *wtree,
                      GtUword pos,
                      GtWtreeSymbol symbol)
{
  gt_assert(wtree != NULL);
  gt_assert(wtree->c_class != NULL);
  if (wtree->c_class->rank_func != NULL)
    return wtree->c_class->rank_func(wtree, pos, symbol);
  return ULONG_MAX;
}

GtUword gt_wtree_select(GtWtree *wtree,
                        GtUword i,
                        GtWtreeSymbol symbol)
{
  gt_assert(wtree != NULL);
  gt_assert(wtree->c_class != NULL);
  if (wtree->c_class->select_func != NULL)
    return wtree->c_class->select_func(wtree, i, symbol);
  return ULONG_MAX;
}

GtUword gt_wtree_length(GtWtree *wtree)
{
  return wtree->members->length;
}

GtUword gt_wtree_num_of_symbols(GtWtree *wtree)
{
  return wtree->members->num_of_symbols;
}

const GtWtreeClass* gt_wtree_class_new(size_t size,
                                       GtWtreeAccessFunc access_func,
                                       GtWtreeRankFunc rank_func,
                                       GtWtreeSelectFunc select_func,
                                       GtWtreeDeleteFunc delete_func)
{
  GtWtreeClass *wtree_c = gt_malloc(sizeof (*wtree_c));
  wtree_c->size = size;
  wtree_c->access_func = access_func;
  wtree_c->rank_func = rank_func;
  wtree_c->select_func = select_func;
  wtree_c-> delete_func = delete_func;
  return wtree_c;
}

void gt_wtree_delete(GtWtree *wtree)
{
  if (wtree != NULL) {
    if (wtree->members->refcount) {
      wtree->members->refcount--;
      return;
    }
    gt_assert(wtree->c_class);
    if (wtree->c_class->delete_func != NULL)
      wtree->c_class->delete_func(wtree);
    gt_free(wtree->members);
    gt_free(wtree);
  }
}
