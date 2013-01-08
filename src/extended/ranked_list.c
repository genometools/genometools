/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/array_api.h"
#include "core/dlist_api.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/yarandom.h"
#include "extended/ranked_list.h"
#include "extended/rbtree.h"

struct GtRankedList
{
  unsigned long currentsize,
                maxsize;
  GtRBTree *root;
  GtCompareWithData comparefunction;
  GtFree free_func;
  void *worstelement,
       *compareinfo;
  GtDlist *list;
};

struct GtRankedListIter
{
  GtDlistelem *current_elem;
};

GtRankedList* gt_ranked_list_new(unsigned long maxsize,
                                 GtCompareWithData comparefunction,
                                 GtFree free_func,
                                 void *compareinfo)
{
  GtRankedList *ranked_list;
  gt_assert(maxsize > 0 && comparefunction != NULL);

  ranked_list = gt_malloc(sizeof (*ranked_list));
  ranked_list->currentsize = 0;
  ranked_list->maxsize = maxsize;
  ranked_list->free_func = free_func;
  ranked_list->comparefunction = comparefunction;
  ranked_list->compareinfo = compareinfo;
  /* ranked_list->root = gt_rbtree_new(comparefunction, free_func, NULL); */
  ranked_list->worstelement = NULL;
  ranked_list->list = gt_dlist_new_with_data(comparefunction, compareinfo);
  return ranked_list;
}

void gt_ranked_list_insert(GtRankedList *ranked_list, void *elem)
{
  /* bool nodecreated = false; */
  if (ranked_list->currentsize < ranked_list->maxsize)
  {
    if (ranked_list->currentsize == 0 ||
        ranked_list->comparefunction(elem,ranked_list->worstelement,
                                     ranked_list->compareinfo) < 0)
    {
      ranked_list->worstelement = elem;
    }
    /* (void) gt_rbtree_search(ranked_list->root, elem, &nodecreated); */
    gt_dlist_add(ranked_list->list, elem);
    /* if (nodecreated)
    { */
    ranked_list->currentsize++;
    /* } */
  } else
  {
/* new element is not as bad as worst element, so insert it and
  and delete the worst element */
    if (ranked_list->comparefunction(ranked_list->worstelement, elem,
                                     ranked_list->compareinfo) <= 0)
    {
      GtDlistelem *oldelem = gt_dlist_first(ranked_list->list);
      if (ranked_list->free_func != NULL)
        ranked_list->free_func(gt_dlistelem_get_data(oldelem));
      gt_dlist_remove(ranked_list->list, oldelem);
      gt_dlist_add(ranked_list->list, elem);
      ranked_list->worstelement =
                       gt_dlistelem_get_data(gt_dlist_first(ranked_list->list));
      /* (void) gt_rbtree_search(ranked_list->root, elem, &nodecreated);
      if (nodecreated)
      {
        if (gt_rbtree_erase(ranked_list->root, ranked_list->worstelement) != 0)
        {
          fprintf(stderr,"%s: deletion failed\n",__func__);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        ranked_list->worstelement = gt_rbtree_minimum_key(ranked_list->root);
      } */
    } else if (ranked_list->free_func != NULL)
    {
      ranked_list->free_func(elem);
    }
  }
}

void* gt_ranked_list_last(const GtRankedList *ranked_list)
{
  /* return gt_rbtree_minimum_key(ranked_list->root); */
  GtDlistelem *elem = gt_dlist_first(ranked_list->list);
  if (elem != NULL)
    return gt_dlistelem_get_data(elem);
  else
    return NULL;
}

void* gt_ranked_list_first(const GtRankedList *ranked_list)
{
  /* return gt_rbtree_maximum_key(ranked_list->root); */
  GtDlistelem *elem = gt_dlist_last(ranked_list->list);
  if (elem != NULL)
    return gt_dlistelem_get_data(elem);
  else
    return NULL;
}

unsigned long gt_ranked_list_size(const GtRankedList *ranked_list)
{
  /* return ranked_list->currentsize; */
  return gt_dlist_size(ranked_list->list);
}

GtRankedListIter* gt_ranked_list_iter_new_from_first(GtRankedList *ranked_list)
{
  /* return gt_rbtree_iter_new_from_last(ranked_list->root); */
  GtRankedListIter *iter = gt_calloc((size_t) 1, sizeof (*iter));
  iter->current_elem = gt_dlist_last(ranked_list->list);
  return iter;
}

GtRankedListIter *gt_ranked_list_iter_new_from_last(GtRankedList *ranked_list)
{
  /* return gt_rbtree_iter_new_from_first(ranked_list->root); */
  GtRankedListIter *iter = gt_calloc((size_t) 1, sizeof (*iter));
  iter->current_elem = gt_dlist_first(ranked_list->list);
  return iter;
}

void *gt_ranked_list_iter_next(GtRankedListIter *ranked_list_iter)
{
  /* return gt_rbtree_iter_prev(ranked_list_iter); */
  void *data = NULL;
  if (ranked_list_iter->current_elem != NULL) {
    data = gt_dlistelem_get_data(ranked_list_iter->current_elem);
    ranked_list_iter->current_elem =
                          gt_dlistelem_previous(ranked_list_iter->current_elem);
  }
  return data;
}

void *gt_ranked_list_iter_prev(GtRankedListIter *ranked_list_iter)
{
 /* return gt_rbtree_iter_next(ranked_list_iter); */
  void *data = NULL;
  if (ranked_list_iter->current_elem != NULL) {
    data =  gt_dlistelem_get_data(ranked_list_iter->current_elem);
    ranked_list_iter->current_elem =
                              gt_dlistelem_next(ranked_list_iter->current_elem);
  }
  return data;
}

void gt_ranked_list_delete(GtRankedList *ranked_list)
{
  if (ranked_list == NULL) return;
  /* gt_rbtree_delete(ranked_list->root); */
  if (ranked_list->free_func != NULL) {
    GtDlistelem *dlistelem;
    for (dlistelem = gt_dlist_first(ranked_list->list); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      void *data = gt_dlistelem_get_data(dlistelem);
      ranked_list->free_func(data);
    }
  }

  gt_dlist_delete(ranked_list->list);

  gt_free(ranked_list);
}

void gt_ranked_list_iter_delete(GtRankedListIter *ranked_list_iter)
{
  if (ranked_list_iter == NULL) return;
  gt_free(ranked_list_iter);
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

typedef struct {
  unsigned long id, score;
} GtRankedListTestStruct;

static int gt_ranked_list_cmp_teststructs(const void *n1, const void *n2,
                                          GT_UNUSED void *info)
{
  GtRankedListTestStruct *l1 = (GtRankedListTestStruct*) n1;
  GtRankedListTestStruct *l2 = (GtRankedListTestStruct*) n2;
  if (l1->score == l2->score)
    return 0;
  else if (l1->score < l2->score)
    return -1;
  else
    return 1;
}

int gt_ranked_list_unit_test(GtError *err)
{
  int had_err = 0;
  GtRankedList *rl;
  GtRankedListIter *iter;
  GtArray *arr;
  const unsigned long nof_best = 30UL, nof_tests = 100UL;
  GtRankedListTestStruct *mystr;
  int values[8] = {-3, 4, 1, 545, 24, 33, 22, 42},
      i, j;
  gt_error_check(err);

  rl = gt_ranked_list_new(5UL, gt_ranked_list_cmp_numbers, NULL, NULL);
  gt_ensure(had_err, rl != NULL);
  gt_ensure(had_err, gt_ranked_list_size(rl) == 0);

  iter = gt_ranked_list_iter_new_from_first(rl);
  mystr = gt_ranked_list_iter_next(iter);
  gt_ensure(had_err, mystr == NULL);
  mystr = gt_ranked_list_iter_next(iter);
  gt_ensure(had_err, mystr == NULL);
  gt_ranked_list_iter_delete(iter);

  iter = gt_ranked_list_iter_new_from_last(rl);
  mystr = gt_ranked_list_iter_prev(iter);
  gt_ensure(had_err, mystr == NULL);
  mystr = gt_ranked_list_iter_prev(iter);
  gt_ensure(had_err, mystr == NULL);
  gt_ranked_list_iter_delete(iter);

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

  for (j = 0; (unsigned long) j < nof_tests; j++) {
    rl = gt_ranked_list_new(30UL, gt_ranked_list_cmp_teststructs, gt_free_func,
                            NULL);
    arr = gt_array_new(sizeof (GtRankedListTestStruct));
    for (i = 0; i < 200; i++) {
      GtRankedListTestStruct newstr,
                             *ptr;
      newstr.id = (unsigned long) i;
      newstr.score = (unsigned long) (random() % (2*nof_best));
      gt_array_add(arr, newstr);
      ptr = gt_malloc(sizeof (*ptr));
      ptr->id = newstr.id;
      ptr->score = newstr.score;
      gt_ranked_list_insert(rl, ptr);
        if ((unsigned long) i < nof_best)
        gt_ensure(had_err, gt_ranked_list_size(rl) == (unsigned long) i + 1UL);
      else
        gt_ensure(had_err, gt_ranked_list_size(rl) == nof_best);
    }
    gt_array_sort_stable_with_data(arr, gt_ranked_list_cmp_teststructs, NULL);
    gt_array_reverse(arr);

    gt_ensure(had_err, gt_ranked_list_size(rl) == nof_best);
    iter = gt_ranked_list_iter_new_from_first(rl);

    i = 0;
    for (mystr = gt_ranked_list_iter_next(iter);
         mystr != NULL;
         mystr = gt_ranked_list_iter_next(iter)) {
      GtRankedListTestStruct *str = (GtRankedListTestStruct*)
                                         gt_array_get(arr, (unsigned long) i++);
      gt_ensure(had_err, mystr != NULL);
      gt_ensure(had_err, mystr->id == str->id);
      gt_ensure(had_err, mystr->score == str->score);
      /* printf("id: %lu/%lu, score %lu/%lu\n", mystr->id, str->id,
                                                mystr->score, str->score); */
    }
    gt_ranked_list_iter_delete(iter);

    gt_array_delete(arr);
    gt_ranked_list_delete(rl);
  }
  return had_err;
}
