/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtcore/queue.h"

typedef struct QueueElem
{
  void *contents;
  struct QueueElem *previous, *next;
} QueueElem;

struct Queue
{
  QueueElem *head,
            *tail;
  unsigned long num_of_elements;
};

Queue* queue_new(void)
{
  return ma_calloc(1, sizeof (Queue));
}

void queue_add(Queue *q, void *contents)
{
  QueueElem *newqueueelem;

  assert(q);

  newqueueelem = ma_malloc(sizeof (QueueElem));
  newqueueelem->contents = contents;
  newqueueelem->previous = NULL;
  newqueueelem->next = q->tail;
  if (q->num_of_elements == 0)
    q->head = newqueueelem;
  else
    q->tail->previous = newqueueelem;
  q->tail = newqueueelem;
  q->num_of_elements++;
}

void* queue_get(Queue *q)
{
  QueueElem *oldheadptr;
  void *contents;

  assert(q && q->num_of_elements);

  oldheadptr = q->head;
  q->head = q->head->previous;
  if (q->head == NULL)
    q->tail = NULL;
  else
    q->head->next = NULL;
  contents = oldheadptr->contents;
  ma_free(oldheadptr);
  q->num_of_elements--;

  return contents;
}

void* queue_head(Queue *q)
{
  assert(q && q->num_of_elements);
  return q->head->contents;
}

int queue_iterate(Queue *q, QueueProcessor queueprocessor, void *info, Error *e)
{
  QueueElem *current;
  error_check(e);
  assert(q && queueprocessor);
  if (q->num_of_elements) {
    for (current = q->head; current; current = current->previous) {
      if (queueprocessor(current->contents, info, e))
        return -1;
    }
  }
  return 0;
}

unsigned long queue_size(const Queue *q)
{
  assert(q);
  return q->num_of_elements;
}

static void queue_wrap(Queue *q, bool freecontents)
{
  if (!q) return;
  while (q->num_of_elements) {
    if (freecontents)
      ma_free(q->head->contents);
    (void) queue_get(q);
  }
  ma_free(q);
}

void queue_delete_with_contents(Queue *q)
{
  queue_wrap(q, true);
}

void queue_delete(Queue *q)
{
  queue_wrap(q, false);
}
