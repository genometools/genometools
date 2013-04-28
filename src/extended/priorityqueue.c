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

#include <stdbool.h>
#include <stdlib.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/error_api.h"
#include "core/mathsupport.h"
#include "core/ensure.h"
#include "priorityqueue.h"

#define GT_MINPQSIZE 16

struct GtPriorityQueue
{
  unsigned long capacity, numofelements;
  GtPriorityQueueElementType minelement, *elements;
};

GtPriorityQueue *gt_priority_queue_new(unsigned long maxnumofelements)
{
  GtPriorityQueue *pq = gt_malloc(sizeof *pq);
  pq->elements = gt_malloc(sizeof (*pq->elements) * (maxnumofelements + 1));
  pq->capacity = maxnumofelements;
  pq->numofelements = 0;
  pq->elements[0].sortkey = 0;
  pq->elements[0].value = 0;
  return pq;
}

static void gt_priority_queue_checkorder(const GtPriorityQueue *pq)
{
  GtPriorityQueueElementType *ptr;

  gt_assert(pq != NULL);
  for (ptr = pq->elements; ptr < pq->elements + pq->numofelements - 1; ptr++)
  {
    gt_assert(ptr->sortkey >= (ptr+1)->sortkey);
  }
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

void gt_priority_queue_add(GtPriorityQueue *pq, unsigned long sortkey,
                          unsigned long value)
{
  gt_assert(pq != NULL && !gt_priority_queue_is_full(pq));
  if (pq->capacity < (unsigned long) GT_MINPQSIZE)
  {
    GtPriorityQueueElementType *ptr;

    /* store elements in reverse order, i.e.\ with the minimum element
       at the last index */
    /* move elements to the right until an element larger or equal than
       the key is found. */
    for (ptr = pq->elements + pq->numofelements; ptr > pq->elements; ptr--)
    {
      if ((ptr-1)->sortkey < sortkey)
      {
        *ptr = *(ptr-1);
      } else
      {
        break;
      }
    }
    ptr->sortkey = sortkey;
    ptr->value = value;
    pq->numofelements++;
  } else
  {
    unsigned long idx = ++pq->numofelements;

    while (true)
    {
      const unsigned long newidx = GT_DIV2(idx);

      if (newidx == 0 || pq->elements[newidx].sortkey <= sortkey)
      {
        break;
      }
      gt_assert(idx > 0);
      pq->elements[idx] = pq->elements[newidx];
      idx = newidx;
    }
    gt_assert(idx > 0);
    pq->elements[idx].sortkey = sortkey;
    pq->elements[idx].value = value;
  }
}

GtPriorityQueueElementType *gt_priority_queue_delete_min(GtPriorityQueue *pq)
{
  gt_assert(pq != NULL && !gt_priority_queue_is_empty(pq));
  if (pq->capacity < (unsigned long) GT_MINPQSIZE)
  {
    gt_assert(pq->numofelements > 0);
    pq->minelement = pq->elements[--pq->numofelements];
    gt_priority_queue_checkorder(pq);
  } else
  {
    unsigned long idx, child;
    GtPriorityQueueElementType lastelement;

    pq->minelement = pq->elements[1];
    lastelement = pq->elements[pq->numofelements--];
    for (idx = 1UL; GT_MULT2(idx) <= pq->numofelements; idx = child)
    {
      child = GT_MULT2(idx);
      gt_assert(child > 0);
      if (child != pq->numofelements &&
          pq->elements[child + 1].sortkey < pq->elements[child].sortkey)
      {
        child++;
      }
      if (lastelement.sortkey > pq->elements[child].sortkey)
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
  return &pq->minelement;
}

const GtPriorityQueueElementType *gt_priority_queue_find_min(const
                                                            GtPriorityQueue *pq)
{
  gt_assert(pq != NULL && !gt_priority_queue_is_empty(pq));
  return pq->elements +
         (pq->capacity < (unsigned long) GT_MINPQSIZE ? pq->numofelements-1
                                                      : 1UL);
}

void gt_priority_queue_delete(GtPriorityQueue *pq)
{
  if (pq != NULL)
  {
    gt_free(pq->elements);
    gt_free(pq);
  }
}

static int cmpUlong(const void *aptr,const void *bptr)
{
  unsigned long a = *((const unsigned long *) aptr);
  unsigned long b = *((const unsigned long *) bptr);

  return a < b ? -1 : (a > b ? 1 : 0);
}

int gt_priority_queue_unit_test(GtError *err)
{
  int had_err = 0;
  unsigned long idx, maxsize = 1000UL,
                trials = 100UL,
                *numbers = gt_malloc(sizeof *numbers * maxsize);

  gt_error_check (err);
  for (idx = 0; idx < trials; idx++)
  {
    unsigned long j,
                  size = gt_rand_max(maxsize),
                  maximal_value = 1 + gt_rand_max(1000UL);
    GtPriorityQueue *pq = gt_priority_queue_new(size);
    GtPriorityQueueElementType *elem;

    for (j = 0; j< size; j++)
    {
      numbers[j] = gt_rand_max(maximal_value);
      gt_priority_queue_add(pq, numbers[j], 0);
    }
    gt_ensure(had_err,gt_priority_queue_is_full(pq));
    qsort(numbers,(size_t) size,sizeof *numbers,cmpUlong);
    for (j = 0; j < size; j++)
    {
      elem = gt_priority_queue_delete_min(pq);
      gt_ensure(had_err,elem->sortkey == numbers[j]);
    }
    gt_ensure(had_err,gt_priority_queue_is_empty(pq));
    gt_priority_queue_delete(pq);
  }
  gt_free(numbers);
  return had_err;
}
