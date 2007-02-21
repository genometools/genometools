/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "array.h"
#include "queue.h"
#include "xansi.h"

struct Queue {
  Array *queue;
  unsigned long current_index;
};

Queue* queue_new(size_t size_of_elem)
{
  Queue *q = xmalloc(sizeof (Queue));

  assert(size_of_elem);

  q->queue = array_new(size_of_elem);
  q->current_index = 0;

  return q;
}

void* queue_get(Queue *q)
{
  void *r;

  assert(q && q->current_index < array_size(q->queue));

  r = array_get(q->queue, q->current_index);
  q->current_index++;
  /* reset */
  if (q->current_index == array_size(q->queue)) {
    array_set_size(q->queue, 0);
    q->current_index = 0;
  }

  return r;
}

void* queue_get_elem(Queue *q, unsigned long idx)
{
  assert(q && q->current_index + idx < array_size(q->queue));
  return array_get(q->queue, q->current_index + idx);
}

void queue_add_elem(Queue *q, void *elem, size_t size_of_elem)
{
  assert(q && elem && size_of_elem);
  array_add_elem(q->queue, elem, size_of_elem);
}

unsigned long queue_size(const Queue *q)
{
  assert(q && q->current_index <= array_size(q->queue));
  return array_size(q->queue) - q->current_index;
}

void queue_delete(Queue *q)
{
  if (!q) return;
  array_delete(q->queue);
  free(q);
}
