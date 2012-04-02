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
#include "priorityqueue.h"

#define GT_MINPQSIZE 10

struct GtPriorityQueue
{
  unsigned long capacity, numofelements;
  GtPriorityQueueElementType minelement, *elements;
};

GtPriorityQueue *gt_priorityqueue_new(unsigned long maxnumofelements)
{
  GtPriorityQueue *pq;

  pq = gt_malloc(sizeof (*pq));
  pq->elements = gt_malloc(sizeof (*pq->elements) * (maxnumofelements + 1));
  pq->capacity = maxnumofelements;
  pq->numofelements = 0;
  pq->elements[0].sortkey = 0;
  pq->elements[0].value = 0;
  return pq;
}

static void gt_priorityqueue_checkorder(const GtPriorityQueue *pq)
{
  GtPriorityQueueElementType *ptr;

  gt_assert(pq != NULL);
  for (ptr = pq->elements; ptr < pq->elements + pq->numofelements - 1; ptr++)
  {
    gt_assert(ptr->sortkey >= (ptr+1)->sortkey);
  }
}

bool gt_priorityqueue_is_empty(const GtPriorityQueue *pq)
{
  gt_assert(pq != NULL);
  return pq->numofelements == 0 ? true : false;
}

bool gt_priorityqueue_is_full(const GtPriorityQueue *pq)
{
  gt_assert(pq != NULL);
  return pq->numofelements == pq->capacity ? true : false;
}

void gt_priorityqueue_add(GtPriorityQueue *pq, unsigned long sortkey,
                          unsigned long value)
{
  gt_assert(pq != NULL && !gt_priorityqueue_is_full(pq));
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
    unsigned long idx;

    for (idx = ++pq->numofelements;
         idx/2 > 0 && pq->elements[idx/2].sortkey > sortkey;
         idx /= 2)
    {
      gt_assert(idx > 0);
      pq->elements[idx] = pq->elements[idx/2];
    }
    gt_assert(idx > 0);
    pq->elements[idx].sortkey = sortkey;
    pq->elements[idx].value = value;
  }
}

GtPriorityQueueElementType *gt_priorityqueue_delete_min(GtPriorityQueue *pq)
{
  gt_assert(pq != NULL && !gt_priorityqueue_is_empty(pq));
  if (pq->capacity < (unsigned long) GT_MINPQSIZE)
  {
    pq->minelement = pq->elements[--pq->numofelements];
    gt_priorityqueue_checkorder(pq);
  } else
  {
    unsigned long idx, child;
    GtPriorityQueueElementType lastelement;

    pq->minelement = pq->elements[1];
    lastelement = pq->elements[pq->numofelements--];
    for (idx = 1UL; idx * 2 <= pq->numofelements; idx = child)
    {
      child = idx * 2;
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

const GtPriorityQueueElementType *gt_priorityqueue_find_min(const
                                                            GtPriorityQueue *pq)
{
  gt_assert(pq != NULL && !gt_priorityqueue_is_empty(pq));
  return pq->elements +
         (pq->capacity < (unsigned long) GT_MINPQSIZE ? pq->numofelements-1
                                                      : 1UL);
}

void gt_priorityqueue_delete(GtPriorityQueue *pq)
{
  if (pq != NULL)
  {
    gt_free(pq->elements);
    gt_free(pq);
  }
}
