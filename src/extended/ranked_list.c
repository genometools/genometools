/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "extended/ranked_list.h"
#include "extended/rbtree.h"

struct GtRankedList
{
  unsigned long currentsize,           /* current size of the ranked list */
                maxsize;               /* maximal size of the ranked list */
  GtRBTree *root;                      /* root of tree */
  GtCompareWithData comparefunction;
  void *worstelement,                  /* reference to worst key */
       *compareinfo;                   /* info needed by compare function */
};

GtRankedList* gt_ranked_list_new(unsigned long maxsize,
                                 GtCompareWithData comparefunction,
                                 void *compareinfo)
{
  GtRankedList *ranked_list;

  ranked_list = gt_malloc(sizeof (*ranked_list));
  ranked_list->currentsize = 0;
  ranked_list->maxsize = maxsize;
  ranked_list->comparefunction = comparefunction;
  ranked_list->compareinfo = compareinfo;
  ranked_list->root = gt_rbtree_new(comparefunction, NULL, NULL);
  ranked_list->worstelement = NULL;
  return ranked_list;
}

void gt_ranked_list_insert(GtRankedList *ranked_list, void *elem)
{
  bool nodecreated = false;

  if (ranked_list->currentsize < ranked_list->maxsize)
  {
    if (ranked_list->currentsize == 0 ||
        ranked_list->comparefunction(elem,ranked_list->worstelement,
                                     ranked_list->compareinfo) < 0)
    {
      ranked_list->worstelement = elem;
    }
    (void) gt_rbtree_search(ranked_list->root, elem, &nodecreated);
    if (nodecreated)
    {
      ranked_list->currentsize++;
    }
  } else
  {
/*
  new element is not as bad as worst element, so insert it and
  and delete the worst element
*/
    if (ranked_list->comparefunction(ranked_list->worstelement, elem,
                                     ranked_list->compareinfo) < 0)
    {
      (void) gt_rbtree_search(ranked_list->root, elem, &nodecreated);
      if (nodecreated)
      {
        if (gt_rbtree_erase(ranked_list->root, ranked_list->worstelement) != 0)
        {
          fprintf(stderr,"%s: deletion failed\n",__func__);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        ranked_list->worstelement = gt_rbtree_minimum_key(ranked_list->root);
      }
    }
  }
}

void* gt_ranked_list_last(const GtRankedList *ranked_list)
{
  return gt_rbtree_minimum_key(ranked_list->root);
}

void* gt_ranked_list_first(const GtRankedList *ranked_list)
{
  return gt_rbtree_maximum_key(ranked_list->root);
}

unsigned long gt_ranked_list_size(const GtRankedList *ranked_list)
{
  return ranked_list->currentsize;
}

GtRankedListIter* gt_ranked_list_iter_new_from_first(GtRankedList *ranked_list)
{
  return gt_rbtree_iter_new_from_first(ranked_list->root);
}

GtRankedListIter *gt_ranked_list_iter_new_from_last(GtRankedList *ranked_list)
{
  return gt_rbtree_iter_new_from_last(ranked_list->root);
}

void *gt_ranked_list_iter_next(GtRankedListIter *ranked_list_iter)
{
  return gt_rbtree_iter_next(ranked_list_iter);
}

void *gt_ranked_list_iter_prev(GtRankedListIter *ranked_list_iter)
{
  return gt_rbtree_iter_prev(ranked_list_iter);
}

void gt_ranked_list_delete(GtRankedList *ranked_list)
{
  gt_rbtree_delete(ranked_list->root);
  gt_free(ranked_list);
}

static int gt_ranked_list_cmp_numbers(const void *n1, const void *n2,
                                      GT_UNUSED void *info)
{
  int l1 = *(int*) n1;
  int l2 = *(int*) n2;
  if (l1 == l2)
    return 0;
  else if (l1 < l2)
    return -1;
  else
    return 1;
}

int gt_ranked_list_unit_test(GtError *err)
{
  int had_err = 0;
  GtRankedList *rl;
  int values[8] = {-3, 4, 1, 545, 24, 33, 22, 42},
      i;
  gt_error_check(err);

  rl = gt_ranked_list_new(5UL, gt_ranked_list_cmp_numbers, NULL);
  gt_ensure(had_err, rl != NULL);
  gt_ensure(had_err, gt_ranked_list_size(rl) == 0);

  for (i = 0; i < 8; i++) {
    gt_ranked_list_insert(rl, values+i);
    if (i < 5)
      gt_ensure(had_err, gt_ranked_list_size(rl) == (unsigned long) i + 1UL);
    else
      gt_ensure(had_err, gt_ranked_list_size(rl) == 5UL);
  }
  gt_ensure(had_err, (*(int*) gt_ranked_list_first(rl)) == 545);
  gt_ensure(had_err, (*(int*) gt_ranked_list_last(rl)) == 22);

  gt_ranked_list_delete(rl);

  return had_err;
}
