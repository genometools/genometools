/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "extended/union_find.h"

#define UNION_FIND_TEST_SIZE  1024

typedef struct {
  GtUword parent/*, rank*/;
} GtUnionFindElement;

struct GtUnionFind {
  GtUnionFindElement *elems;
  GtUword num_of_elems;
  GtUword allocated;
};

GtUnionFind* gt_union_find_new(GtUword num_of_elems)
{
  GtUnionFind *uf;
  GtUword i;
  gt_assert(num_of_elems);
  uf = gt_malloc(sizeof *uf);
  uf->elems = gt_calloc(sizeof (GtUnionFindElement), num_of_elems);
  uf->allocated = num_of_elems;
  for (i = 0; i < num_of_elems; i++)
    uf->elems[i].parent = i;
  uf->num_of_elems = num_of_elems;
  return uf;
}

void gt_union_find_reset(GtUnionFind *uf, GtUword num_of_elems)
{
  GtUword i;
  gt_assert(num_of_elems);
  if (num_of_elems > uf->allocated)
  {
    uf->elems = gt_realloc(uf->elems, sizeof (GtUnionFindElement) *
        num_of_elems);
    uf->allocated = num_of_elems;
  }
  for (i = 0; i < num_of_elems; i++)
    uf->elems[i].parent = i;
  uf->num_of_elems = num_of_elems;
}

void gt_union_find_delete(GtUnionFind *uf)
{
  if (!uf) return;
  gt_free(uf->elems);
  gt_free(uf);
}

GtUword gt_union_find_find(GtUnionFind *uf, GtUword elem)
{
  gt_assert(uf && elem < uf->num_of_elems);
  if (elem != uf->elems[elem].parent) /* path compression */
    uf->elems[elem].parent = gt_union_find_find(uf, uf->elems[elem].parent);
  return uf->elems[elem].parent;
}

void gt_union_find_union(GtUnionFind *uf, GtUword elem_a,
                         GtUword elem_b)
{
  GtUword x, y;
  gt_assert(uf && elem_a < uf->num_of_elems && elem_b < uf->num_of_elems);
  x = gt_union_find_find(uf, elem_a);
  y = gt_union_find_find(uf, elem_b);
  uf->elems[y].parent = x;
  /* union-by-rank heuristic
  if (uf->elems[x].rank > uf->elems[y].rank)
    uf->elems[y].parent = x;
  else {
    uf->elems[x].parent = y;
    if (uf->elems[x].rank == uf->elems[y].rank)
      uf->elems[y].rank++;
  }
  */
}

int gt_union_find_unit_test(GtError *err)
{
  GtUword i;
  GtUnionFind *uf;
  int had_err = 0;
  gt_error_check(err);

  /* one element */
  uf = gt_union_find_new(1);
  gt_ensure(gt_union_find_find(uf, 0) == 0);
  gt_union_find_delete(uf);

  /* two elements */
  if (!had_err) {
    uf = gt_union_find_new(2);
    gt_ensure(gt_union_find_find(uf, 0) != gt_union_find_find(uf, 1));
    gt_union_find_union(uf, 0, 1);
    gt_ensure(gt_union_find_find(uf, 0) == gt_union_find_find(uf, 1));
    gt_union_find_delete(uf);
  }

  /* many elements */
  if (!had_err) {
    uf = gt_union_find_new(UNION_FIND_TEST_SIZE);
    for (i = 1; !had_err && i < UNION_FIND_TEST_SIZE; i++)
      gt_ensure(gt_union_find_find(uf, 0) != gt_union_find_find(uf,
                                                                         i));
    for (i = 1; !had_err && i < UNION_FIND_TEST_SIZE; i++)
      gt_union_find_union(uf, 0, i);
    for (i = 1; !had_err && i < UNION_FIND_TEST_SIZE; i++)
      gt_ensure(gt_union_find_find(uf, 0) == gt_union_find_find(uf,
                                                                         i));
    gt_union_find_delete(uf);
  }

  return had_err;
}
