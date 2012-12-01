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

#include "core/ma_api.h"
#include "extended/rbtree.h"
#include "extended/ranked_list.h"

struct GtRankedlist
{
  unsigned long currentsize,           /* current size of the ranked list */
                maxsize;               /* maximal size of the ranked list */
  GtRBTree *root;                     /* root of tree */
  GtRankedlistCompareFunc comparefunction;
  void *worstelement,                  /* reference to worst key */
       *compareinfo;                   /* info needed by compare function */
};

GtRankedlist *gt_ranked_list_new(unsigned long maxsize,
                                 GtRankedlistCompareFunc comparefunction,
                                 void *compareinfo)
{
  GtRankedlist *ranked_list;

  ranked_list = gt_malloc(sizeof (*ranked_list));
  ranked_list->currentsize = 0;
  ranked_list->maxsize = maxsize;
  ranked_list->comparefunction = comparefunction;
  ranked_list->compareinfo = compareinfo;
  ranked_list->root = gt_rbtree_new(comparefunction, NULL, NULL);
  ranked_list->worstelement = NULL;
  return ranked_list;
}

void gt_ranked_list_insert(GtRankedlist *ranked_list,void *elemin)
{
  bool nodecreated = false;

  if (ranked_list->currentsize < ranked_list->maxsize)
  {
    if (ranked_list->currentsize == 0 ||
        ranked_list->comparefunction(elemin,ranked_list->worstelement,
                                     ranked_list->compareinfo) < 0)
    {
      ranked_list->worstelement = elemin;
    }
    (void) gt_rbtree_search(ranked_list->root, elemin, &nodecreated);
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
    if (ranked_list->comparefunction(ranked_list->worstelement,elemin,
                                     ranked_list->compareinfo) < 0)
    {
      (void) gt_rbtree_search(ranked_list->root, elemin, &nodecreated);
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

void *gt_ranked_list_minimum_key(const GtRankedlist *ranked_list)
{
  return gt_rbtree_minimum_key(ranked_list->root);
}

void *gt_ranked_list_maximum_key(const GtRankedlist *ranked_list)
{
  return gt_rbtree_maximum_key(ranked_list->root);
}

unsigned long gt_ranked_list_currentsize(const GtRankedlist *ranked_list)
{
  return ranked_list->currentsize;
}

GtRankedListIter *gt_ranked_list_iter_new_from_first(GtRankedlist *ranked_list)
{
  return gt_rbtree_iter_new_from_first(ranked_list->root);
}

GtRankedListIter *gt_ranked_list_iter_new_from_last(GtRankedlist *ranked_list)
{
  return gt_rbtree_iter_new_from_last(ranked_list->root);
}

void *gt_ranked_list_iter_next(GtRankedListIter *trav)
{
  return gt_rbtree_iter_next(trav);
}

void *gt_ranked_list_iter_prev(GtRankedListIter *trav)
{
  return gt_rbtree_iter_prev(trav);
}

void gt_ranked_list_delete(GtRankedlist *ranked_list)
{
  gt_rbtree_delete(ranked_list->root);
  gt_free(ranked_list);
}
