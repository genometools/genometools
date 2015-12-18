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

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/error_api.h"
#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "extended/priority_queue.h"

#define GT_HEAP_PARENT(X)  GT_DIV2(X)
#define GT_HEAP_LEFT(X)    GT_MULT2(X)
#define GT_MINPQSIZE       16

struct GtPriorityQueue
{
  GtCompare cmpfun;
  GtUword capacity, numofelements;
  void *minelement, **elements;
};

GtPriorityQueue *gt_priority_queue_new(GtCompare cmpfun,
                                       GtUword maxnumofelements)
{
  GtPriorityQueue *pq = gt_malloc(sizeof *pq);
  pq->elements = gt_calloc((maxnumofelements + 1), sizeof (*pq->elements));
  pq->minelement = NULL;
  pq->cmpfun = cmpfun;
  pq->capacity = maxnumofelements;
  pq->numofelements = 0;
  return pq;
}

bool gt_priority_queue_is_empty(const GtPriorityQueue *pq)
{
  gt_assert(pq != NULL);
  return pq->numofelements == 0 ? true : false;
}

bool gt_priority_queue_is_full(const GtPriorityQueue *pq)
{
  gt_assert(pq != NULL);
  return pq->numofelements == pq->capacity ? true : false;
}

void gt_priority_queue_add(GtPriorityQueue *pq, void *value)
{
  gt_assert(pq != NULL && !gt_priority_queue_is_full(pq));
  if (pq->capacity < (GtUword) GT_MINPQSIZE)
  {
    void **ptr;

    /* store elements in reverse order, i.e.\ with the minimum element
       at the last index */
    /* move elements to the right until an element larger or equal than
       the key is found. */
    for (ptr = pq->elements + pq->numofelements; ptr > pq->elements; ptr--)
    {
      if (*(ptr-1) == NULL || pq->cmpfun(*(ptr-1),value) < 0)
      {
        *ptr = *(ptr-1);
      } else
      {
        break;
      }
    }
    *ptr = value;
    pq->numofelements++;
  } else
  {
    GtUword idx = ++pq->numofelements;

    while (true)
    {
      const GtUword newidx = GT_HEAP_PARENT(idx);

      if (newidx == 0 || pq->cmpfun(pq->elements[newidx],value) <= 0)
      {
        break;
      }
      gt_assert(idx > 0);
      pq->elements[idx] = pq->elements[newidx];
      idx = newidx;
    }
    gt_assert(idx > 0);
    pq->elements[idx] = value;
  }
}

void *gt_priority_queue_extract_min(GtPriorityQueue *pq)
{
  gt_assert(pq != NULL && !gt_priority_queue_is_empty(pq));
  if (pq->capacity < (GtUword) GT_MINPQSIZE)
  {
    gt_assert(pq->numofelements > 0);
    pq->minelement = pq->elements[--pq->numofelements];
  } else
  {
    GtUword idx, child;
    void *lastelement;

    pq->minelement = pq->elements[1];
    lastelement = pq->elements[pq->numofelements--];
    for (idx = 1UL; GT_HEAP_LEFT(idx) <= pq->numofelements; idx = child)
    {
      child = GT_HEAP_LEFT(idx);
      gt_assert(child > 0);
      if (child != pq->numofelements &&
          pq->cmpfun(pq->elements[child + 1],pq->elements[child]) < 0)
      {
        child++;
      }
      if (pq->cmpfun(lastelement,pq->elements[child]) > 0)
      {
        pq->elements[idx] = pq->elements[child];
      } else
      {
        break;
      }
    }
    gt_assert(idx > 0);
    pq->elements[idx] = lastelement;
  }
  return pq->minelement;
}

const void *gt_priority_queue_find_min(const GtPriorityQueue *pq)
{
  gt_assert(pq != NULL && !gt_priority_queue_is_empty(pq));
  return *(pq->elements +
         (pq->capacity < (GtUword) GT_MINPQSIZE ? pq->numofelements-1
                                                      : 1UL));
}

void gt_priority_queue_delete(GtPriorityQueue *pq)
{
  if (pq != NULL)
  {
    gt_free(pq->elements);
    gt_free(pq);
  }
}

static int gt_priority_queue_cmpulong(const void *aptr,const void *bptr)
{
  GtUword a, b;

  gt_assert(aptr != NULL && bptr != NULL);
  a = *((const GtUword *) aptr);
  b = *((const GtUword *) bptr);
  return a < b ? -1 : (a > b ? 1 : 0);
}

static void gt_priority_sort(GtUword *numbers,GtUword len)
{
  GtUword j, previousvalue = ULONG_MAX;
  GtPriorityQueue *pq = gt_priority_queue_new(gt_priority_queue_cmpulong, len);

  for (j = 0; j < len; j++)
  {
    gt_priority_queue_add(pq, numbers + j);
  }
  gt_assert(gt_priority_queue_is_full(pq));
  for (j = 0; j < len; j++)
  {
    void *elem = gt_priority_queue_extract_min(pq);

    if (previousvalue != ULONG_MAX)
    {
      gt_assert(previousvalue <= *((GtUword *) elem));
    }
    previousvalue = *((GtUword *) elem);
  }
  gt_assert(gt_priority_queue_is_empty(pq));
  gt_priority_queue_delete(pq);
}

int gt_priority_queue_unit_test(GtError *err)
{
  int had_err = 0;
  GtUword idx,
          maxsize = 10000UL,
          trials = 1000UL,
          *numbers = gt_malloc(sizeof *numbers * maxsize),
          *numbers_copy = gt_malloc(sizeof *numbers_copy * maxsize);
  GtUword arr[] = {76UL, 132UL, 136UL, 538UL, 545UL, 401UL};
  GtPriorityQueue *tmp = NULL;
  gt_error_check(err);

  tmp = gt_priority_queue_new(gt_priority_queue_cmpulong, 42);
  gt_ensure(gt_priority_queue_is_empty(tmp));
  gt_priority_queue_delete(tmp);

  tmp = gt_priority_queue_new(gt_priority_queue_cmpulong, 0);
  gt_ensure(gt_priority_queue_is_empty(tmp));
  gt_priority_queue_delete(tmp);

  if (!had_err) {
    gt_priority_sort(arr, (GtUword) sizeof arr/sizeof arr[0]);
    for (idx = 0; !had_err && idx < trials; idx++)
    {
      GtUword j,
              size = gt_rand_max(maxsize),
              maximal_value = 1 + gt_rand_max(1000UL);
      GtPriorityQueue *pq = gt_priority_queue_new(gt_priority_queue_cmpulong,
                                                  size);
      gt_ensure(gt_priority_queue_is_empty(pq));

      for (j = 0; !had_err && j < size; j++)
      {
        numbers_copy[j] = numbers[j] = gt_rand_max(maximal_value);
        gt_priority_queue_add(pq, numbers_copy + j);
      }
      gt_ensure(gt_priority_queue_is_full(pq));
      qsort(numbers,(size_t) size, sizeof *numbers, gt_priority_queue_cmpulong);

      for (j = 0; !had_err && j < size; j++)
      {
        GtUword this_min, elem;
        this_min = *((GtUword*) gt_priority_queue_find_min(pq));
        elem = *((GtUword *) gt_priority_queue_extract_min(pq));
        gt_ensure(elem == this_min);
        gt_ensure(elem == numbers[j]);
      }
      gt_ensure(gt_priority_queue_is_empty(pq));
      gt_priority_queue_delete(pq);
    }
  }
  gt_free(numbers);
  gt_free(numbers_copy);

  return had_err;
}
