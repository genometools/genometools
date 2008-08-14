/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

#include "libgtcore/interval_tree.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtview/line_breaker_bases.h"
#include "libgtview/line_breaker_rep.h"

struct LineBreakerBases {
  const LineBreaker parent_instance;
  Hashtable *itrees;
};

#define line_breaker_bases_cast(LB)\
        line_breaker_cast(line_breaker_bases_class(), LB)

bool line_breaker_bases_is_line_occupied(LineBreaker* lb,
                                         Line *line,
                                         Block *block)
{
  LineBreakerBases *lbb;
  Range r;
  IntervalTree *t;
  assert(lb && block && line);
  r = block_get_range(block);
  lbb = line_breaker_bases_cast(lb);
  if (!(t = hashtable_get(lbb->itrees, line)))
    return false;
  else
    return (interval_tree_find_first_overlapping(t, r.start, r.end));
}

void line_breaker_bases_register_block(LineBreaker *lb,
                                       Line *line,
                                       Block *block)
{
  LineBreakerBases *lbb;
  IntervalTree *t;
  IntervalTreeNode *new_node;
  Range *rng;
  assert(lb && block && line);
  lbb = line_breaker_bases_cast(lb);
  rng = block_get_range_ptr(block);
  new_node = interval_tree_node_new(rng, rng->start, rng->end);
  if (!(t = hashtable_get(lbb->itrees, line)))
  {
    t = interval_tree_new(NULL);
    hashtable_add(lbb->itrees, line, t);
  }
  interval_tree_insert(t, new_node);
}

void line_breaker_bases_delete(LineBreaker *lb)
{
  LineBreakerBases *lbb;
  if (!lb) return;
  lbb = line_breaker_bases_cast(lb);
  hashtable_delete(lbb->itrees);
}

const LineBreakerClass* line_breaker_bases_class(void)
{
  static const LineBreakerClass line_breaker_class =
    { sizeof (LineBreakerBases),
      line_breaker_bases_is_line_occupied,
      line_breaker_bases_register_block,
      line_breaker_bases_delete };
  return &line_breaker_class;
}

LineBreaker* line_breaker_bases_new()
{
  LineBreakerBases *lbb;
  LineBreaker *lb;
  lb = line_breaker_create(line_breaker_bases_class());
  lbb = line_breaker_bases_cast(lb);
  lbb->itrees = hashtable_new(HASH_DIRECT, NULL,
                              (FreeFunc) interval_tree_delete);
  return lb;
}
