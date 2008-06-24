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

#include <assert.h>
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtext/union_find.h"

#define UNION_FIND_TEST_SIZE  1024

typedef struct {
  unsigned long parent,
                rank;
} UnionFindElement;

struct UnionFind {
  UnionFindElement *elems;
  unsigned long num_of_elems;
};

UnionFind* union_find_new(unsigned long num_of_elems)
{
  UnionFind *uf;
  unsigned long i;
  assert(num_of_elems);
  uf = ma_malloc(sizeof *uf);
  uf->elems = ma_calloc(sizeof (UnionFindElement), num_of_elems);
  for (i = 0; i < num_of_elems; i++)
    uf->elems[i].parent = i;
  uf->num_of_elems = num_of_elems;
  return uf;
}

void union_find_delete(UnionFind *uf)
{
  if (!uf) return;
  ma_free(uf->elems);
  ma_free(uf);
}

unsigned long union_find_find(UnionFind *uf, unsigned long elem)
{
  assert(uf && elem < uf->num_of_elems);
  if (elem != uf->elems[elem].parent) /* path compression */
    uf->elems[elem].parent = union_find_find(uf, uf->elems[elem].parent);
  return uf->elems[elem].parent;
}

void union_find_union(UnionFind *uf, unsigned long elem_a, unsigned long elem_b)
{
  unsigned long x, y;
  assert(uf && elem_a < uf->num_of_elems && elem_b < uf->num_of_elems);
  x = union_find_find(uf, elem_a);
  y = union_find_find(uf, elem_b);
  /* union-by-rank heuristic */
  if (uf->elems[x].rank > uf->elems[y].rank)
    uf->elems[y].parent = x;
  else {
    uf->elems[x].parent = y;
    if (uf->elems[x].rank == uf->elems[y].rank)
      uf->elems[y].rank++;
  }
}

int union_find_unit_test(Error *err)
{
  unsigned long i;
  UnionFind *uf;
  int had_err = 0;
  error_check(err);

  /* one element */
  uf = union_find_new(1);
  ensure(had_err, union_find_find(uf, 0) == 0);
  union_find_delete(uf);

  /* two elements */
  if (!had_err) {
    uf = union_find_new(2);
    ensure(had_err, union_find_find(uf, 0) != union_find_find(uf, 1));
    union_find_union(uf, 0, 1);
    ensure(had_err, union_find_find(uf, 0) == union_find_find(uf, 1));
    union_find_delete(uf);
  }

  /* many elements */
  if (!had_err) {
    uf = union_find_new(UNION_FIND_TEST_SIZE);
    for (i = 1; !had_err && i < UNION_FIND_TEST_SIZE; i++)
      ensure(had_err, union_find_find(uf, 0) != union_find_find(uf, i));
    for (i = 1; !had_err && i < UNION_FIND_TEST_SIZE; i++)
      union_find_union(uf, 0, i);
    for (i = 1; !had_err && i < UNION_FIND_TEST_SIZE; i++)
      ensure(had_err, union_find_find(uf, 0) == union_find_find(uf, i));
    union_find_delete(uf);
  }

  return had_err;
}
