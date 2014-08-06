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

#ifndef WTREE_H
#define WTREE_H

#include "core/types_api.h"

/* Abstract class GtWtree,
   Used to store sequences or Permutations.
   Each implementation must implement a function for mapping <GtWtreeSymbol> to
   the implementation-specific alphabet symbol. */
typedef struct GtWtree GtWtree;

/* Type used by GtWtree to represent symbols */
typedef GtUword GtWtreeSymbol;

/* Increases the reference count of the <GtWtree>. */
GtWtree*      gt_wtree_ref(GtWtree *wtree);

/* Returns the symbol at <pos> from <wtree>, note that <pos> < length of
   <wtree>.
   Returns ULONG_MAX if function is not implemented. */
GtWtreeSymbol gt_wtree_access(GtWtree *wtree,
                              GtUword pos);

/* Returns the number of symbols <symbol> in the prefix of <wtree> upto position
   <pos>. Note that <pos> < length of <wtree>.
   Returns ULONG_MAX if function is not implemented. */
GtUword       gt_wtree_rank(GtWtree *wtree,
                            GtUword pos,
                            GtWtreeSymbol symbol);

/* Returns the position of the <i>th <symbol> in wtree, returns ULONG_MAX if the
   <wtree> contains less than <i> symbols. Note that 0 < <i> <= length of
   <wtree>.
   Returns ULONG_MAX if function is not implemented. */
GtUword       gt_wtree_select(GtWtree *wtree,
                              GtUword i,
                              GtWtreeSymbol symbol);

/* Returns the length of the sequence/permutation encoded by <wtree>. */
GtUword       gt_wtree_length(GtWtree *wtree);

/* Returns the number of distinct symbols recognized by <wtree>. Corresponding
   to alphabet or set-size. */
GtUword       gt_wtree_num_of_symbols(GtWtree *wtree);

/* Frees all memory associated with <wtree>. */
void          gt_wtree_delete(GtWtree *wtree);

#endif
