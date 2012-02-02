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
  GtPQelementtype *elements;
};

GtPriorityQueue *priorityqueue_new(unsigned long maxnumofelements)
{
  GtPriorityQueue *pq;

  pq = gt_malloc(sizeof (*pq));
  pq->elements = gt_malloc(sizeof (*pq->elements) * (maxnumofelements + 1));
  pq->capacity = maxnumofelements;
  pq->numofelements = 0;
  pq->elements[0] = 0;
  return pq;
}

static void priorityqueue_checkorder(const GtPriorityQueue *pq)
{
  GtPQelementtype *ptr;

  for (ptr = pq->elements; ptr < pq->elements + pq->numofelements - 1; ptr++)
  {
    gt_assert(*ptr >= *(ptr+1));
  }
}

void priorityqueue_add(GtPriorityQueue *pq, GtPQelementtype value)
{
  gt_assert(!priorityqueue_is_full(pq));
  if (pq->capacity < (unsigned long) GT_MINPQSIZE)
  {
    GtPQelementtype *ptr, tmp;

    pq->elements[pq->numofelements++] = value;
    for (ptr = pq->elements + pq->numofelements - 1; ptr > pq->elements; ptr--)
    {
      if (*(ptr-1) < *ptr)
      {
        tmp = *ptr;
        *ptr = *(ptr-1);
        *(ptr-1) = tmp;
      }
    }
    priorityqueue_checkorder(pq);
  } else
  {
    unsigned long idx;

    for (idx = ++pq->numofelements; idx/2 > 0 && pq->elements[idx/2] > value;
         idx /= 2)
    {
      gt_assert(idx > 0);
      pq->elements[idx] = pq->elements[idx/2];
    }
    gt_assert(idx > 0);
    pq->elements[idx] = value;
  }
}

GtPQelementtype priorityqueue_delete_min(GtPriorityQueue *pq)
{
  GtPQelementtype minelement;

  gt_assert(!priorityqueue_is_empty(pq));
  if (pq->capacity < (unsigned long) GT_MINPQSIZE)
  {
    minelement = pq->elements[--pq->numofelements];
    priorityqueue_checkorder(pq);
  } else
  {
    unsigned long idx, child;
    GtPQelementtype lastelement;

    minelement = pq->elements[1];
    lastelement = pq->elements[pq->numofelements--];
    for (idx = 1UL; idx * 2 <= pq->numofelements; idx = child)
    {
      child = idx * 2;
      gt_assert(child > 0);
      if (child != pq->numofelements &&
          pq->elements[child + 1] < pq->elements[child])
      {
        child++;
      }
      if (lastelement > pq->elements[child])
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
  return minelement;
}

GtPQelementtype priorityqueue_find_min(const GtPriorityQueue *pq)
{
  gt_assert(!priorityqueue_is_empty(pq));
  return pq->elements[pq->capacity < (unsigned long) GT_MINPQSIZE ? 0 : 1];
}

bool priorityqueue_is_empty(const GtPriorityQueue *pq)
{
  return pq->numofelements == 0 ? true : false;
}

bool priorityqueue_is_full(const GtPriorityQueue *pq)
{
  return pq->numofelements == pq->capacity ? true : false;
}

void priorityqueue_delete(GtPriorityQueue *pq)
{
  gt_free(pq->elements);
  gt_free(pq);
}
