/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef QUEUE_INLINE_H
#define QUEUE_INLINE_H

#include <stdio.h>
#include <stdlib.h>
#include "core/assert_api.h"
#include "divmodmul.h"
#include "spacedef.h"

typedef struct
{
 Inl_Queueelem *queuespace;      /* the space to store the queue elements */
 unsigned long enqueueindex, /* entry into which element is to be enqued */
               dequeueindex, /* last element of queue */
               queuesize,    /* size of the queue */
               noofelements; /* no ofelements between enqueueindex+1
                                and dequeindex */
} Inl_Queue;

typedef int (*Inl_Queueprocessor)(Inl_Queueelem *,void *info);

/*
  The following function delivers an empty queue with a reservoir of
  \texttt{size} elements to be stored in the queue. The
  reservoir can, if necessary, be enlarged.
*/

static Inl_Queue *gt_inl_queue_new(unsigned long queuesize)
{
  Inl_Queue *q;

  q = gt_malloc(sizeof(*q));
  gt_assert(queuesize > 0);
  ALLOCASSIGNSPACE(q->queuespace,NULL,Inl_Queueelem,queuesize);
  q->noofelements = 0;
  q->queuesize = queuesize;
  q->dequeueindex = q->enqueueindex = queuesize - 1;
  return q;
}

/*
  The following function frees the space required for the queue.
*/

void gt_inl_queue_delete(Inl_Queue *q)
{
  FREESPACE(q->queuespace);
  FREESPACE(q);
}

/*
  The following function returns true iff the queue is empty.
*/

bool gt_inl_queue_isempty(const Inl_Queue *q)
{
  return (q->noofelements == 0) ? true : false;
}

/*
  The following function resizes the queue by doubling the
  space reservoir.
*/

static void doublequeuesize(Inl_Queue *q)
{
  unsigned long i, j, newsize;

  newsize = MULT2(q->queuesize); /* double the size */
  ALLOCASSIGNSPACE(q->queuespace,q->queuespace,Inl_Queueelem,newsize);
  j=q->queuesize;
  gt_assert(q->enqueueindex >= q->dequeueindex);
  if (q->enqueueindex >= q->dequeueindex)
  {
    for (i=q->enqueueindex+1; i<=q->queuesize-1; i++)
    {
      q->queuespace[j++] = q->queuespace[i];
    }
    for (i=0; i<=q->dequeueindex; i++)
    {
      q->queuespace[j++] = q->queuespace[i];
    }
  } else
  {
    gt_assert(q->dequeueindex - q->enqueueindex == q->queuesize);
    printf("queuesize=%lu,enqueue=%lu,dequeue=%lu\n",
           q->queuesize,
           q->dequeueindex,
           q->enqueueindex);
    gt_assert(q->dequeueindex == q->queuesize-1);
    gt_assert(q->enqueueindex == 0);
    for (i=q->enqueueindex+1; i<=q->dequeueindex; i++)
    {
      q->queuespace[j++] = q->queuespace[i];
    }
  }
  q->dequeueindex = newsize-1;
  q->enqueueindex = q->queuesize-1;
  q->queuesize = newsize;
}

static void constextendqueuesize(Inl_Queue *q)
{
  const unsigned long addconst = q->dequeueindex+1;
  unsigned long i, j, newsize = q->queuesize + addconst; /* double the size */
  ALLOCASSIGNSPACE(q->queuespace,q->queuespace,Inl_Queueelem,newsize);
  gt_assert(q->enqueueindex >= q->dequeueindex);
  j=q->queuesize;
  for (i=0; i<=q->dequeueindex; i++)
  {
    q->queuespace[j++] = q->queuespace[i];
  }
  q->dequeueindex = newsize-1;
  q->enqueueindex = q->queuesize-1;
  printf("from queue of size %lu to queue of size %lu\n",q->queuesize,newsize);
  q->queuesize = newsize;
}

/*
  The following function adds an element \texttt{elem} to the end of
  the queue.
*/

void gt_inl_queue_add(Inl_Queue *q,Inl_Queueelem elem)
{
  if (q->noofelements == q->queuesize)
  {
    bool useconst = true;

    if (useconst)
    {
      constextendqueuesize(q);
    } else
    {
      doublequeuesize(q);
    }
  }
  q->noofelements++;
  q->queuespace[q->enqueueindex] = elem;
  if (q->enqueueindex > 0)
  {
    q->enqueueindex--;
  } else
  {
    q->enqueueindex = q->queuesize - 1; /* dequeuindex < queuesize-1 */
  }
}

/*
  The following function removes the element \texttt{elem} from the
  start of the queue.
*/

Inl_Queueelem gt_inl_queue_get(Inl_Queue *q)
{
  Inl_Queueelem value;

  gt_assert(q->noofelements > 0);
  q->noofelements--;
  value = q->queuespace[q->dequeueindex];
  if (q->dequeueindex > 0)
  {
    q->dequeueindex--;
  } else
  {
    q->dequeueindex = q->queuesize - 1;  /* != enqueueindex, since at least
                                            one elem */
  }
  return value;
}

Inl_Queueelem *gt_inl_queue_head(const Inl_Queue *q)
{
  gt_assert(q->noofelements > 0);
  return &q->queuespace[q->dequeueindex];
}

Inl_Queueelem *gt_inl_queue_tail(const Inl_Queue *q)
{
  gt_assert(q->noofelements > 0);
  if (q->enqueueindex == q->queuesize-1)
  {
    return &q->queuespace[0];
  }
  return &q->queuespace[q->enqueueindex+1];
}

void gt_inl_queue_deletehead(Inl_Queue *q)
{
  gt_assert(q->noofelements > 0);
  q->noofelements--;
  if (q->dequeueindex > 0)
  {
    q->dequeueindex--;
  } else
  {
    q->dequeueindex = q->queuesize - 1; /* != enqueueindex, since at least
                                           one elem */
  }
}

int gt_inl_queue_iterate(const Inl_Queue *q,
                         Inl_Queueprocessor queueprocessor,void *info)
{
  unsigned long idx;

  if (q->noofelements > 0)
  {
    if (q->enqueueindex < q->dequeueindex)
    {
      for (idx=q->enqueueindex+1; idx<=q->dequeueindex; idx++)
      {
        if (queueprocessor(q->queuespace + idx,info) != 0)
        {
          return -1;
        }
      }
    } else
    {
      for (idx=q->enqueueindex+1; idx<=q->queuesize-1; idx++)
      {
        if (queueprocessor(q->queuespace + idx,info) != 0)
        {
          return -1;
        }
      }
      for (idx=0; idx<=q->dequeueindex; idx++)
      {
        if (queueprocessor(q->queuespace + idx,info) != 0)
        {
          return -1;
        }
      }
    }
  }
  return 0;
}
#endif
