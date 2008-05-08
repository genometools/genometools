/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef UNION_FIND_H
#define UNION_FIND_H

#include "libgtcore/error.h"

/*
  This class implements the union-find data structure (with union-by-rank
  heuristic and path compression).
  For a description see for example page 446 to page 449 of the book

  T.H. Cormen, C.E. Leiserson and R.L. Rivest. Introduction to Algorithms.
  MIT Press: Cambridge, MA, 1990.
*/

typedef struct UnionFind UnionFind;

/* Create a new union-find data structures representing <num_of_elems> many
   elements (numbered from 0 up to <num_of_elems> - 1) contained in disjoined
   sets. */
UnionFind*    union_find_new(unsigned long num_of_elems);

/* Delete the given union-find data structure. */
void          union_find_delete(UnionFind*);

/* Find the representative set for the given <elem>. */
unsigned long union_find_find(UnionFind*, unsigned long elem);

/* Union the set containing <elem_a> with the set containing <elem_b>. */
void          union_find_union(UnionFind*, unsigned long elem_a,
                                           unsigned long elem_b);

int           union_find_unit_test(Error*);

#endif
