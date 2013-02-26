/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/class_alloc_lock.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "annotationsketch/line_breaker_bases.h"
#include "annotationsketch/line_breaker_rep.h"

struct GtLineBreakerBases {
  const GtLineBreaker parent_instance;
  GtHashmap *itrees;
};

#define gt_line_breaker_bases_cast(LB)\
        gt_line_breaker_cast(gt_line_breaker_bases_class(), LB)

int gt_line_breaker_bases_is_line_occupied(GtLineBreaker* lb, bool *result,
                                           GtLine *line, GtBlock *block,
                                           GT_UNUSED GtError *err)
{
  GtLineBreakerBases *lbb;
  GtRange r;
  GtIntervalTree *t;
  gt_assert(lb && block && line);
  r = gt_block_get_range(block);
  lbb = gt_line_breaker_bases_cast(lb);
  if (!(t = gt_hashmap_get(lbb->itrees, line)))
    *result = false;
  else
    *result = (gt_interval_tree_find_first_overlapping(t, r.start, r.end));
  return 0; /* gt_line_breaker_bases_is_line_occupied() is sane */
}

int gt_line_breaker_bases_register_block(GtLineBreaker *lb,
                                         GtLine *line,
                                         GtBlock *block,
                                         GT_UNUSED GtError *err)
{
  GtLineBreakerBases *lbb;
  GtIntervalTree *t;
  GtIntervalTreeNode *new_node;
  GtRange *rng;
  gt_assert(lb && block && line);
  lbb = gt_line_breaker_bases_cast(lb);
  rng = gt_block_get_range_ptr(block);
  new_node = gt_interval_tree_node_new(rng, rng->start, rng->end);
  if (!(t = gt_hashmap_get(lbb->itrees, line)))
  {
    t = gt_interval_tree_new(NULL);
    gt_hashmap_add(lbb->itrees, line, t);
  }
  gt_interval_tree_insert(t, new_node);
  return 0; /* gt_line_breaker_bases_register_block() is sane */
}

void gt_line_breaker_bases_delete(GtLineBreaker *lb)
{
  GtLineBreakerBases *lbb;
  if (!lb) return;
  lbb = gt_line_breaker_bases_cast(lb);
  gt_hashmap_delete(lbb->itrees);
}

const GtLineBreakerClass* gt_line_breaker_bases_class(void)
{
  static const GtLineBreakerClass *lbc = NULL;
  gt_class_alloc_lock_enter();
  if (!lbc)
  {
    lbc = gt_line_breaker_class_new(sizeof (GtLineBreakerBases),
                                    gt_line_breaker_bases_is_line_occupied,
                                    gt_line_breaker_bases_register_block,
                                    gt_line_breaker_bases_delete);
  }
  gt_class_alloc_lock_leave();
  return lbc;
}

GtLineBreaker* gt_line_breaker_bases_new()
{
  GtLineBreakerBases *lbb;
  GtLineBreaker *lb;
  lb = gt_line_breaker_create(gt_line_breaker_bases_class());
  lbb = gt_line_breaker_bases_cast(lb);
  lbb->itrees = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                               (GtFree) gt_interval_tree_delete);
  return lb;
}
