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
#include "core/dynalloc.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/queue.h"
#include "core/unused_api.h"

struct GtQueue {
  void **contents;
  long front, /* f */
       back,  /* b */
       size;
  size_t allocated;
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

GtQueue* gt_queue_new(void)
{
  return gt_calloc(1, sizeof (GtQueue));
}

static void check_space(GtQueue *q)
{
  if (!q->allocated) { /* empty queue without allocated memory */
    q->contents = gt_dynalloc(q->contents, &q->allocated, sizeof (void*));
    q->size = q->allocated / sizeof (void*);
  }
  else if (q->front < q->back) { /* no wraparound */
    if (q->back == q->size) {
      if (q->front)
        q->back = 0; /* perform wraparound */
      else { /* extend contents buffer */
        q->contents = gt_dynalloc(q->contents, &q->allocated,
                                  q->allocated + sizeof (void*));
        q->size = q->allocated / sizeof (void*);
      }
    }
  }
  else if (q->back && (q->back == q->front)) { /* wraparound */
    q->contents = gt_dynalloc(q->contents, &q->allocated,
                              q->allocated + q->front * sizeof (void*));
    memcpy(q->contents + q->size, q->contents, q->front * sizeof (void*));
    /* dynalloc() always doubles the already allocated memory region, which
       means we always have some additional space after the copied memory region
       left to set the back pointer to (otherwise we would have to reset the
       back pointer to 0).
     */
    gt_assert((size_t) q->front + q->size < q->allocated / sizeof (void*));
    q->back = q->front + q->size;
    q->size = q->allocated / sizeof (void*);
  }
}

void gt_queue_add(GtQueue *q, void *contents)
{
  gt_assert(q);
  check_space(q);
  q->contents[q->back++] = contents;
}

void* gt_queue_get(GtQueue *q)
{
  void *contents;
  gt_assert(q && gt_queue_size(q));
  /* get contents */
  contents = q->contents[q->front++];
  /* adjust indices */
  if (q->front == q->back)
    q->front = q->back = 0; /* reset */
  else if (q->front == q->size)
    q->front = 0; /* wraparound */
  return contents;
}

void* gt_queue_head(GtQueue *q)
{
  gt_assert(q && gt_queue_size(q));
  return q->contents[q->front];
}

void gt_queue_remove(GtQueue *q, void *elem)
{
  long i, elemidx;
  gt_assert(q &&  gt_queue_size(q));
  if (q->front < q->back) { /* no wraparound */
    for (i = q->back-1; i >= q->front; i--) {
      if (q->contents[i] == elem)
        break;
    }
    elemidx = i;
    gt_assert(elemidx >= q->front); /* valid element found */
    for (i = elemidx+1; i < q->back; i++)
      q->contents[i-1] = q->contents[i];
    q->contents[q->back-1] = NULL;
    q->back--;
    if (q->front == q->back)
      q->front = q->back = 0; /* reset */
  }
  else { /* wraparound */
    for (i = q->back-1; i >= 0; i--) {
      if (q->contents[i] == elem)
        break;
    }
    elemidx = i;
    if (elemidx >= 0) { /* element found */
      for (i = elemidx+1; i < q->back; i++)
        q->contents[i-1] = q->contents[i];
      q->contents[q->back-1] = NULL;
      q->back--;
      if (q->back == 0) q->back = q->size;
      return;
    }
    for (i = q->size-1; i >= q->front; i--) {
      if (q->contents[i] == elem)
        break;
    }
    elemidx = i;
    gt_assert(elemidx >= q->front); /* valid element found */
    for (i = elemidx+1; i < q->size; i++)
      q->contents[i-1] = q->contents[i];
    q->contents[q->size-1] = q->contents[0];
    for (i = 1; i < q->back; i++)
      q->contents[i-1] = q->contents[i];
    q->contents[q->back-1] = NULL;
    q->back--;
    if (q->back == 0)
      q->back = q->size;
  }
}

int gt_queue_iterate(GtQueue *q, GtQueueProcessor gt_queue_processor,
                     void *info, GtError *err)
{
  long i;
  int rval;
  gt_error_check(err);
  gt_assert(q && gt_queue_processor);
  if (gt_queue_size(q)) {
    if (q->front < q->back) { /* no wraparound */
      for (i = q->front; i < q->back; i++) {
        if ((rval = gt_queue_processor(q->contents + i, info, err)))
          return rval;
      }
    }
    else { /* wraparound */
      for (i = q->front; i < q->size; i++) {
        if ((rval = gt_queue_processor(q->contents + i, info, err)))
          return rval;
      }
      for (i = 0; i < q->back; i++) {
        if ((rval = gt_queue_processor(q->contents+i, info, err))) return rval;
      }
    }
  }
  return 0;
}

int gt_queue_iterate_reverse(GtQueue *q, GtQueueProcessor gt_queue_processor,
                             void *info, GtError *err)
{
  long i;
  int rval;
  gt_error_check(err);
  gt_assert(q && gt_queue_processor);
  if (gt_queue_size(q)) {
    if (q->front < q->back) { /* no wraparound */
      for (i = q->back-1; i >= q->front; i--) {
        if ((rval = gt_queue_processor(q->contents + i, info, err)))
          return rval;
      }
    }
    else { /* wraparound */
      for (i = q->back-1; i >= 0; i--) {
        if ((rval = gt_queue_processor(q->contents + i, info, err)))
          return rval;
      }
      for (i = q->size-1; i >= q->front; i--) {
        if ((rval = gt_queue_processor(q->contents +i, info, err))) return rval;
      }
    }
  }
  return 0;
}

unsigned long gt_queue_size(const GtQueue *q)
{
  gt_assert(q);
  if ((q->front < q->back) || ((q->front == 0) && (q->back == 0)))
    return q->back - q->front; /* no wraparound */
  return q->size - (q->front - q->back); /* wraparound */
}

static void queue_wrap(GtQueue *q, bool free_contents)
{
  gt_assert(q);
  if (free_contents) {
    while (gt_queue_size(q))
      gt_free(gt_queue_get(q));
  }
  gt_free(q->contents);
  gt_free(q);
}

void gt_queue_delete(GtQueue *q)
{
  if (!q) return;
  queue_wrap(q, false);
}

void gt_queue_delete_with_contents(GtQueue *q)
{
  if (!q) return;
  queue_wrap(q, true);
}

static int check_queue(void **elem, void *info, GtError *err)
{
  long *check_counter = info;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(check_counter);
  gt_ensure(had_err, *check_counter == *(long*) elem);
  if (!had_err)
    (*check_counter)++;
  return had_err;
}

static int check_queue_reverse(void **elem, void *info, GtError *err)
{
  long *check_counter_reverse = info;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(check_counter_reverse);
  gt_ensure(had_err, *check_counter_reverse == *(long*) elem);
  if (!had_err)
    (*check_counter_reverse)--;
  return had_err;
}

static int fail_func(GT_UNUSED void **elem, GT_UNUSED void *info,
                     GT_UNUSED GtError *err)
{
  return -1;
}

int gt_queue_unit_test(GtError *err)
{
  long check_counter = 0, check_counter_reverse = 1023;
  unsigned long i;
  int had_err = 0;
  GtQueue *q;

  gt_error_check(err);

  /* without wraparound */
  q = gt_queue_new();
  gt_ensure(had_err, !gt_queue_size(q));
  for (i = 0; !had_err && i < 1024; i++) {
    gt_queue_add(q, (void*) i);
    gt_ensure(had_err, gt_queue_size(q) == i + 1);
  }
  if (!had_err)
    had_err = gt_queue_iterate(q, check_queue, &check_counter, err);
  if (!had_err) {
    had_err = gt_queue_iterate_reverse(q, check_queue_reverse,
                                    &check_counter_reverse, err);
  }
  gt_ensure(had_err, gt_queue_iterate(q, fail_func, NULL, NULL));
  gt_ensure(had_err, gt_queue_iterate_reverse(q, fail_func, NULL, NULL));
  if (!had_err) {
    gt_queue_remove(q, (void*) 0);
    gt_ensure(had_err, gt_queue_size(q) == 1023);
  }
  for (i = 1; !had_err && i < 1024; i++) {
    gt_ensure(had_err, gt_queue_head(q) == (void*) i);
    gt_ensure(had_err, gt_queue_get(q) == (void*) i);
    gt_ensure(had_err, gt_queue_size(q) == 1024 - i - 1);
  }
  gt_ensure(had_err, !gt_queue_size(q));
  gt_queue_delete(q);

  /* with wraparound (without full queue) */
  if (!had_err) {
    q = gt_queue_new();
    gt_ensure(had_err, !gt_queue_size(q));
    for (i = 0; !had_err && i < 1024; i++) {
      gt_queue_add(q, (void*) i);
      gt_ensure(had_err, gt_queue_size(q) == i + 1);
    }
    check_counter = 0;
    check_counter_reverse = 1023;
    if (!had_err)
      had_err = gt_queue_iterate(q, check_queue, &check_counter, err);
    gt_ensure(had_err, gt_queue_iterate(q, fail_func, NULL, NULL));
    gt_ensure(had_err, gt_queue_iterate_reverse(q, fail_func, NULL, NULL));
    if (!had_err) {
      had_err = gt_queue_iterate_reverse(q, check_queue_reverse,
                                         &check_counter_reverse, err);
    }
    for (i = 0; !had_err && i < 512; i++) {
      gt_ensure(had_err, gt_queue_head(q) == (void*) i);
      gt_ensure(had_err, gt_queue_get(q) == (void*) i);
      gt_ensure(had_err, gt_queue_size(q) == 1024 - i - 1);
    }
    for (i = 0; !had_err && i < 512; i++) {
      gt_queue_add(q, (void*) (i + 1024));
      gt_ensure(had_err, gt_queue_size(q) == 512 + i + 1);
    }
    check_counter = 512;
    check_counter_reverse = 1535;
    if (!had_err)
      had_err = gt_queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err) {
      had_err = gt_queue_iterate_reverse(q, check_queue_reverse,
                                         &check_counter_reverse, err);
    }
    gt_ensure(had_err, gt_queue_iterate(q, fail_func, NULL, NULL));
    gt_ensure(had_err, gt_queue_iterate_reverse(q, fail_func, NULL, NULL));
    if (!had_err) {
      gt_queue_remove(q, (void*) 512);
      gt_ensure(had_err, gt_queue_size(q) == 1023);
    }
    for (i = 1; !had_err && i < 1024; i++) {
      gt_ensure(had_err, gt_queue_head(q) == (void*) (512 + i));
      gt_ensure(had_err, gt_queue_get(q) == (void*) (512 + i));
      gt_ensure(had_err, gt_queue_size(q) == 1024 - i - 1);
    }
    gt_ensure(had_err, !gt_queue_size(q));
    gt_queue_delete(q);
  }

  /* with wraparound (with full queue) */
  if (!had_err) {
    q = gt_queue_new();
    gt_ensure(had_err, !gt_queue_size(q));
    for (i = 0; !had_err && i < 1024; i++) {
      gt_queue_add(q, (void*) i);
      gt_ensure(had_err, gt_queue_size(q) == i + 1);
    }
    check_counter = 0;
    check_counter_reverse = 1023;
    if (!had_err)
      had_err = gt_queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err) {
      had_err = gt_queue_iterate_reverse(q, check_queue_reverse,
                                      &check_counter_reverse, err);
    }
    gt_ensure(had_err, gt_queue_iterate(q, fail_func, NULL, NULL));
    gt_ensure(had_err, gt_queue_iterate_reverse(q, fail_func, NULL, NULL));
    for (i = 0; !had_err && i < 512; i++) {
      gt_ensure(had_err, gt_queue_head(q) == (void*) i);
      gt_ensure(had_err, gt_queue_get(q) == (void*) i);
      gt_ensure(had_err, gt_queue_size(q) == 1024 - i - 1);
    }
    for (i = 0; !had_err && i < 1024; i++) {
      gt_queue_add(q, (void*) (i + 1024));
      gt_ensure(had_err, gt_queue_size(q) == 512 + i + 1);
    }
    check_counter = 512;
    check_counter_reverse = 2047;
    if (!had_err)
      had_err = gt_queue_iterate(q, check_queue, &check_counter, err);
    if (!had_err) {
      had_err = gt_queue_iterate_reverse(q, check_queue_reverse,
                                      &check_counter_reverse, err);
    }
    gt_ensure(had_err, gt_queue_iterate(q, fail_func, NULL, NULL));
    gt_ensure(had_err, gt_queue_iterate_reverse(q, fail_func, NULL, NULL));
    if (!had_err) {
      gt_queue_remove(q, (void*) 512);
      gt_ensure(had_err, gt_queue_size(q) == 1535);
    }
    for (i = 1; !had_err && i < 1536; i++) {
      gt_ensure(had_err, gt_queue_head(q) == (void*) (512 + i));
      gt_ensure(had_err, gt_queue_get(q) == (void*) (512 + i));
      gt_ensure(had_err, gt_queue_size(q) == 1536 - i - 1);
    }
    gt_ensure(had_err, !gt_queue_size(q));
    gt_queue_delete(q);
  }

  /* test a corner case */
  if (!had_err) {
    q = gt_queue_new();
    gt_queue_add(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 1);
    if (!had_err)
      gt_queue_add(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 2);
    gt_ensure(had_err, gt_queue_get(q));
    gt_ensure(had_err, gt_queue_size(q) == 1);
    if (!had_err)
      gt_queue_add(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 2);
    gt_ensure(had_err, gt_queue_get(q));
    gt_ensure(had_err, gt_queue_size(q) == 1);
    if (!had_err)
      gt_queue_add(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 2);
    gt_ensure(had_err, gt_queue_get(q));
    gt_ensure(had_err, gt_queue_size(q) == 1);
    gt_ensure(had_err, gt_queue_get(q));
    gt_ensure(had_err, gt_queue_size(q) == 0);
    if (!had_err)
      gt_queue_add(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 1);
    gt_ensure(had_err, gt_queue_get(q));
    gt_ensure(had_err, gt_queue_size(q) == 0);
    gt_queue_delete(q);
  }

  /* gt_queue_remove() corner case */
  if (!had_err) {
    q = gt_queue_new();
    gt_queue_add(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 1);
    gt_queue_remove(q, (void*) 1);
    gt_ensure(had_err, gt_queue_size(q) == 0);
    gt_queue_delete(q);
  }

  /* gt_queue_remove() corner case */
  if (!had_err) {
    q = gt_queue_new();
    gt_queue_add(q, (void*) 0);
    gt_queue_add(q, (void*) 1);
    gt_queue_add(q, (void*) 2);
    gt_queue_add(q, (void*) 3);
    gt_ensure(had_err, gt_queue_get(q) == (void*) 0);
    gt_ensure(had_err, gt_queue_get(q) == (void*) 1);
    gt_queue_add(q, (void*) 4);
    gt_queue_add(q, (void*) 5);
    gt_queue_remove(q, (void*) 4);
    gt_queue_remove(q, (void*) 2);
    gt_queue_remove(q, (void*) 5);
    gt_queue_remove(q, (void*) 3);
    gt_ensure(had_err, gt_queue_size(q) == 0);
    gt_queue_delete(q);
  }

  /* delete with contents */
  if (!had_err) {
    q = gt_queue_new();
    gt_ensure(had_err, !gt_queue_size(q));
    if (!had_err)
      gt_queue_add(q, gt_calloc(1, 16));
    gt_ensure(had_err, gt_queue_size(q) == 1);
    if (!had_err)
      gt_queue_add(q, gt_calloc(1, 32));
    gt_ensure(had_err, gt_queue_size(q) == 2);
    gt_queue_delete_with_contents(q);
  }

  return had_err;
}
