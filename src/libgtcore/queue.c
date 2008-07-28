/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "libgtcore/dynalloc.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/queue.h"
#include "libgtcore/unused.h"

struct Queue
{
  void **contents;
  unsigned long front, /* f */
                back;  /* b */
  size_t size,
         allocated;
};

/*
  empty queue:
  ---------------------------------------
  f
  b

  filled queue (no wraparound):
  ---#############-----------------------
     f            b

  filled queue (before wraparound):
  ---####################################
     f                                   b

  full queue (no wraparound):
  #######################################
  f                                      b

  filled queue (wraparound):
  ###-------#############################
     b      f

  full queue (wraparound):
  #######################################
            f
            b
*/

Queue* queue_new(void)
{
  return ma_calloc(1, sizeof (Queue));
}

static void queue_wrap(Queue *q, bool free_contents)
{
  assert(q);
  if (free_contents) {
    while (queue_size(q))
      ma_free(queue_get(q));
  }
  ma_free(q->contents);
  ma_free(q);
}

void queue_delete(Queue *q)
{
  if (!q) return;
  queue_wrap(q, false);
}

void queue_delete_with_contents(Queue *q)
{
  if (!q) return;
  queue_wrap(q, true);
}

static void check_space(Queue *q)
{
  if (!q->allocated) { /* empty queue without allocated memory */
    q->contents = dynalloc(q->contents, &q->allocated, sizeof (void*));
    q->size = q->allocated / sizeof (void*);
  }
  else if (q->front < q->back) { /* no wraparound */
    if (q->back == q->size) {
      if (q->front)
        q->back = 0; /* perform wraparound */
      else { /* extend contents buffer */
        q->contents = dynalloc(q->contents, &q->allocated,
                               q->allocated + sizeof (void*));
        q->size = q->allocated / sizeof (void*);
      }
    }
  }
  else if (q->back && (q->back == q->front)) { /* wraparound */
    q->contents = dynalloc(q->contents, &q->allocated,
                           q->allocated + q->front * sizeof (void*));
    memcpy(q->contents + q->size, q->contents, q->front * sizeof (void*));
    /* dynalloc() always doubles the already allocated memory region, which
       means we always have some additional space after the copied memory region
       left to set the back pointer to (otherwise we would have to reset the
       back pointer to 0).
     */
    assert(q->front + q->size < q->allocated / sizeof (void*));
    q->back = q->front + q->size;
    q->size = q->allocated / sizeof (void*);
  }
}

void queue_add(Queue *q, void *contents)
{
  assert(q);
  check_space(q);
  q->contents[q->back++] = contents;
}

void* queue_get(Queue *q)
{
  void *contents;
  assert(q && queue_size(q));
  /* get contents */
  contents = q->contents[q->front++];
  /* adjust indices */
  if (q->front == q->back)
    q->front = q->back = 0; /* reset */
  else if (q->front == q->size)
    q->front = 0; /* wraparound */
  return contents;
}

void* queue_head(Queue *q)
{
  assert(q && queue_size(q));
  return q->contents[q->front];
}

int queue_iterate(Queue *q, QueueProcessor queueprocessor, void *info,
                  Error *err)
{
  unsigned long i;
  error_check(err);
  assert(q && queueprocessor);
  if (queue_size(q)) {
    if (q->front < q->back) { /* no wraparound */
      for (i = q->front; i < q->back; i++) {
        if (queueprocessor(q->contents[i], info, err))
          return -1;
      }
    }
    else { /* wraparound */
      for (i = q->front; i < q->size; i++) {
        if (queueprocessor(q->contents[i], info, err))
          return -1;
      }
      for (i = 0; i < q->back; i++) {
        if (queueprocessor(q->contents[i], info, err)) return -1;
      }
    }
  }
  return 0;
}

unsigned long queue_size(const Queue *q)
{
  assert(q);
  if ((q->front < q->back) || ((q->front == 0) && (q->back == 0)))
    return q->back - q->front; /* no wraparound */
  return q->size - (q->front - q->back); /* wraparound */
}

static int check_queue(void *elem, void *info, Error *err)
{
  unsigned long *check_counter = info;
  int had_err = 0;
  error_check(err);
  assert(check_counter);
  ensure(had_err, *check_counter == (unsigned long) elem);
  if (!had_err)
    (*check_counter)++;
  return had_err;
}

static int fail_func(UNUSED void *elem, UNUSED void *info, UNUSED Error *err)
{
  return -1;
}

int queue_unit_test(Error *err)
{
  unsigned long check_counter = 0;
  unsigned long i;
  int had_err = 0;
  Queue *q;

  error_check(err);

  /* without wraparound */
  q = queue_new();
  ensure(had_err, !queue_size(q));
  for (i = 0; !had_err && i < 1024; i++) {
    queue_add(q, (void*) i);
    ensure(had_err, queue_size(q) == i + 1);
  }
  if (!had_err)
    had_err = queue_iterate(q, check_queue, &check_counter, err);
  if (!had_err)
    ensure(had_err, queue_iterate(q, fail_func, NULL, NULL));
  for (i = 0; !had_err && i < 1024; i++) {
    ensure(had_err, queue_head(q) == (void*) i);
    ensure(had_err, queue_get(q) == (void*) i);
    ensure(had_err, queue_size(q) == 1024 - i - 1);
  }
  ensure(had_err, !queue_size(q));
  queue_delete(q);

  /* with wraparound (without full queue) */
  if (!had_err) {
    q = queue_new();
    ensure(had_err, !queue_size(q));
    for (i = 0; !had_err && i < 1024; i++) {
      queue_add(q, (void*) i);
      ensure(had_err, queue_size(q) == i + 1);
    }
    check_counter = 0;
    if (!had_err)
      had_err = queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err)
      ensure(had_err, queue_iterate(q, fail_func, NULL, NULL));
    for (i = 0; !had_err && i < 512; i++) {
      ensure(had_err, queue_head(q) == (void*) i);
      ensure(had_err, queue_get(q) == (void*) i);
      ensure(had_err, queue_size(q) == 1024 - i - 1);
    }
    for (i = 0; !had_err && i < 512; i++) {
      queue_add(q, (void*) (i + 1024));
      ensure(had_err, queue_size(q) == 512 + i + 1);
    }
    check_counter = 512;
    if (!had_err)
      had_err = queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err)
      ensure(had_err, queue_iterate(q, fail_func, NULL, NULL));
    for (i = 0; !had_err && i < 1024; i++) {
      ensure(had_err, queue_head(q) == (void*) (512 + i));
      ensure(had_err, queue_get(q) == (void*) (512 + i));
      ensure(had_err, queue_size(q) == 1024 - i - 1);
    }
    ensure(had_err, !queue_size(q));
    queue_delete(q);
  }

  /* with wraparound (with full queue) */
  if (!had_err) {
    q = queue_new();
    ensure(had_err, !queue_size(q));
    for (i = 0; !had_err && i < 1024; i++) {
      queue_add(q, (void*) i);
      ensure(had_err, queue_size(q) == i + 1);
    }
    check_counter = 0;
    if (!had_err)
      had_err = queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err)
      ensure(had_err, queue_iterate(q, fail_func, NULL, NULL));
    for (i = 0; !had_err && i < 512; i++) {
      ensure(had_err, queue_head(q) == (void*) i);
      ensure(had_err, queue_get(q) == (void*) i);
      ensure(had_err, queue_size(q) == 1024 - i - 1);
    }
    for (i = 0; !had_err && i < 1024; i++) {
      queue_add(q, (void*) (i + 1024));
      ensure(had_err, queue_size(q) == 512 + i + 1);
    }
    check_counter = 512;
    if (!had_err)
      had_err = queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err)
      ensure(had_err, queue_iterate(q, fail_func, NULL, NULL));
    for (i = 0; !had_err && i < 1536; i++) {
      ensure(had_err, queue_head(q) == (void*) (512 + i));
      ensure(had_err, queue_get(q) == (void*) (512 + i));
      ensure(had_err, queue_size(q) == 1536 - i - 1);
    }
    ensure(had_err, !queue_size(q));
    queue_delete(q);
  }

  /* test a corner case */
  if (!had_err) {
    q = queue_new();
    queue_add(q, (void*) 1);
    ensure(had_err, queue_size(q) == 1);
    if (!had_err)
      queue_add(q, (void*) 1);
    ensure(had_err, queue_size(q) == 2);
    ensure(had_err, queue_get(q));
    ensure(had_err, queue_size(q) == 1);
    if (!had_err)
      queue_add(q, (void*) 1);
    ensure(had_err, queue_size(q) == 2);
    ensure(had_err, queue_get(q));
    ensure(had_err, queue_size(q) == 1);
    if (!had_err)
      queue_add(q, (void*) 1);
    ensure(had_err, queue_size(q) == 2);
    ensure(had_err, queue_get(q));
    ensure(had_err, queue_size(q) == 1);
    ensure(had_err, queue_get(q));
    ensure(had_err, queue_size(q) == 0);
    if (!had_err)
      queue_add(q, (void*) 1);
    ensure(had_err, queue_size(q) == 1);
    ensure(had_err, queue_get(q));
    ensure(had_err, queue_size(q) == 0);
    queue_delete(q);
  }

  /* delete with contents */
  if (!had_err) {
    q = queue_new();
    ensure(had_err, !queue_size(q));
    if (!had_err)
      queue_add(q, ma_calloc(1, 16));
    ensure(had_err, queue_size(q) == 1);
    if (!had_err)
      queue_add(q, ma_calloc(1, 32));
    ensure(had_err, queue_size(q) == 2);
    queue_delete_with_contents(q);
  }

  return had_err;
}
