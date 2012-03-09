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
#include "core/divmodmul.h"

typedef struct
{
 GtInl_Queueelem *queuespace;  /* the space to store the queue elements */
 unsigned long enqueueindex, /* entry into which element is to be enqued */
               dequeueindex, /* last element of queue */
               queuesize,    /* size of the queue */
               noofelements; /* no ofelements between enqueueindex+1
                                and dequeindex */
} GtInl_Queue;

typedef int (*GtInl_Queueprocessor)(GtInl_Queueelem *,void *info);

/*
  The following function delivers an empty queue with a reservoir of
  \texttt{size} elements to be stored in the queue. The
  reservoir can, if necessary, be enlarged.
*/

static inline GtInl_Queue *gt_inl_queue_new(unsigned long queuesize)
{
  GtInl_Queue *q;

  q = gt_malloc(sizeof (*q));
  gt_assert(queuesize > 0);
  q->queuespace = gt_malloc(sizeof (*q->queuespace) * queuesize);
  q->noofelements = 0;
  q->queuesize = queuesize;
  q->dequeueindex = q->enqueueindex = queuesize - 1;
  return q;
}

/*
  The following function frees the space required for the queue.
*/

static inline void gt_inl_queue_delete(GtInl_Queue *q)
{
  if (q != NULL)
  {
    gt_free(q->queuespace);
    gt_free(q);
  }
}

/*
  The following function returns true iff the queue is empty.
*/

static inline bool gt_inl_queue_isempty(const GtInl_Queue *q)
{
  gt_assert(q != NULL);
  return (q->noofelements == 0) ? true : false;
}

/*
  The following function resizes the queue by doubling the
  space reservoir.
*/

static inline void extendqueuesize(GtInl_Queue *q,bool doublesize)
{
  unsigned long addconst, idx, newsize;

  gt_assert(q != NULL);
  if (doublesize)
  {
    addconst = q->queuesize;
  } else
  {
    addconst = MIN(1024UL,q->queuesize);
  }
  newsize = q->queuesize + addconst;
  q->queuespace = gt_realloc(q->queuespace,sizeof (*q->queuespace) * newsize);
  gt_assert(q->enqueueindex == q->dequeueindex);
  gt_assert(addconst > 0);
  for (idx=q->queuesize-1; idx>q->enqueueindex; idx--)
  {
    q->queuespace[idx+addconst] = q->queuespace[idx];
  }
  q->enqueueindex += addconst;
  /*
  printf("from queue of size %lu to queue of size %lu\n",q->queuesize,newsize);
  printf("now enqueindex=%lu,dequeuindex=%lu\n",
         q->enqueueindex,q->dequeueindex);
  */
  q->queuesize = newsize;
}

/*
  The following function adds an element \texttt{elem} to the end of
  the queue.
*/

static inline void gt_inl_queue_add(GtInl_Queue *q, GtInl_Queueelem elem,
                                    bool doublesize)
{
  gt_assert(q != NULL);
  if (q->noofelements == q->queuesize)
  {
    extendqueuesize(q,doublesize);
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

static inline GtInl_Queueelem gt_inl_queue_get(GtInl_Queue *q)
{
  GtInl_Queueelem value;

  gt_assert(q != NULL && q->noofelements > 0);
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

/*@unused@*/ static inline GtInl_Queueelem *gt_inl_queue_head(
                                                     const GtInl_Queue *q)
{
  gt_assert(q != NULL && q->noofelements > 0);
  return q->queuespace + q->dequeueindex;
}

/*@unused@*/ static inline GtInl_Queueelem *gt_inl_queue_tail(
                                                      const GtInl_Queue *q)
{
  gt_assert(q != NULL && q->noofelements > 0);
  if (q->enqueueindex == q->queuesize-1)
  {
    return q->queuespace;
  }
  return q->queuespace + q->enqueueindex + 1;
}

/*@unused@*/ static inline void gt_inl_queue_deletehead(GtInl_Queue *q)
{
  gt_assert(q != NULL && q->noofelements > 0);
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

/*@unused@*/ static inline int gt_inl_queue_iterate(const GtInl_Queue *q,
                                           GtInl_Queueprocessor queueprocessor,
                                           void *info)
{
  gt_assert(q != NULL);
  if (q->noofelements > 0)
  {
    unsigned long idx;

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
