/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/thread.h"
#include "core/unused_api.h"

#ifdef GT_THREADS_ENABLED

#include <pthread.h>

GtThread* gt_thread_new(GtThreadFunc function, void *data, GtError *err)
{
  GtThread *thread;
  int rval;
  gt_error_check(err);
  gt_assert(function);
  thread = gt_malloc(sizeof (pthread_t));
  rval = pthread_create((pthread_t*) thread, NULL, function, data);
  if (rval) {
    gt_error_set(err, "cannot create thread: %s\n", strerror(rval));
    gt_free(thread);
    return NULL;
  }
  return thread;
}

GtMutex* gt_mutex_new(void)
{
  GtMutex* mutex;
  int rval;
  mutex = gt_malloc(sizeof (pthread_mutex_t));
  /* initialize mutex with default attributes */
  rval = pthread_mutex_init((pthread_mutex_t*) mutex, NULL);
  gt_assert(!rval);
  return mutex;
}

void gt_mutex_delete(GtMutex *mutex)
{
  int rval;
  if (!mutex) return;
  rval = pthread_mutex_destroy((pthread_mutex_t*) mutex);
  gt_assert(!rval);
  gt_free(mutex);
}

void gt_mutex_lock_func(GtMutex *mutex)
{
 int rval = pthread_mutex_lock((pthread_mutex_t*) mutex);
 gt_assert(!rval);
}

void gt_mutex_unlock_func(GtMutex *mutex)
{
  int rval = pthread_mutex_unlock((pthread_mutex_t*) mutex);
  gt_assert(!rval);
}

#else

GtThread* gt_thread_new(GtThreadFunc function, void *data, GtError *err)
{
  GtThread *thread;
  gt_error_check(err);
  gt_assert(function);
  thread = gt_malloc(sizeof (void*));
  function(data);
  return thread;
}

GtMutex* gt_mutex_new(void)
{
  return NULL;
}

void gt_mutex_delete(GT_UNUSED GtMutex *mutex)
{
  return;
}

#endif

void gt_thread_delete(GtThread *thread)
{
  if (!thread);
  gt_free(thread);
}
