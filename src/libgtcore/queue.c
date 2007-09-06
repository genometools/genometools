/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <libgtcore/array.h>
#include <libgtcore/queue.h>
#include <libgtcore/xansi.h>

struct Queue {
  Array *queue;
  unsigned long current_index;
};

Queue* queue_new(size_t size_of_elem, Env *env)
{
  Queue *q = env_ma_malloc(env, sizeof (Queue));
  assert(size_of_elem);
  q->queue = array_new(size_of_elem, env);
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
    array_reset(q->queue);
    q->current_index = 0;
  }
  return r;
}

void* queue_get_elem(Queue *q, unsigned long idx)
{
  assert(q && q->current_index + idx < array_size(q->queue));
  return array_get(q->queue, q->current_index + idx);
}

void queue_add_elem(Queue *q, void *elem, size_t size_of_elem, Env *env)
{
  assert(q && elem && size_of_elem);
  array_add_elem(q->queue, elem, size_of_elem, env);
}

unsigned long queue_size(const Queue *q)
{
  assert(q && q->current_index <= array_size(q->queue));
  return array_size(q->queue) - q->current_index;
}

void queue_delete(Queue *q, Env *env)
{
  if (!q) return;
  array_delete(q->queue, env);
  env_ma_free(q, env);
}
