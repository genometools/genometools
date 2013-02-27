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

#include <errno.h>
#include <string.h>
#include "core/array_api.h"
#include "core/ma_api.h"
#include "core/thread_api.h"
#include "core/unused_api.h"

unsigned int gt_jobs = 1;

#ifdef GT_THREADS_ENABLED

#include <pthread.h>

int gt_multithread(GtThreadFunc function, void *data, GtError *err)
{
  GtArray *threads;
  GtThread *thread;
  unsigned int i, j;

  gt_error_check(err);
  gt_assert(function);

  threads = gt_array_new(sizeof (GtThread*));

  /* start all other threads and store them */
  for (i = 1; i < gt_jobs; i++) {
    if (!(thread = gt_thread_new(function, data, err))) {
      for (j = 0; j < gt_array_size(threads); j++)
        gt_thread_delete(*(GtThread**) gt_array_get(threads, j));
      gt_array_delete(threads);
      return -1;
    }
    gt_array_add(threads, thread);
  }

  function(data); /* execute function in main thread, too */

  /* wait until all other threads are finished */
  for (i = 0; i < gt_array_size(threads); i++) {
    thread = *(GtThread**) gt_array_get(threads, i);
    gt_thread_join(thread);
    gt_thread_delete(thread);
  }
  gt_array_delete(threads);

  return 0;
}

GtThread* gt_thread_new(GtThreadFunc function, void *data, GtError *err)
{
  GtThread *thread;
  int rval;
  gt_error_check(err);
  gt_assert(function);
  thread = gt_malloc(sizeof (pthread_t));
  rval = pthread_create((pthread_t*) thread, NULL, function, data);
  if (rval) {
    gt_error_set(err, "cannot create thread: %s", strerror(rval));
    gt_free(thread);
    return NULL;
  }
  return thread;
}

void gt_thread_join(GtThread *thread)
{
  GT_UNUSED int rval;
  gt_assert(thread);
  rval = pthread_join(*(pthread_t*) thread, NULL);
  gt_assert(!rval); /* XXX */
}

static void* thread_xmalloc(size_t size, const char *filename, int line)
{
  void *p;
  if ((p = malloc(size)) == NULL) {
    fprintf(stderr, "cannot malloc(%zu) memory: %s\n", size, strerror(errno));
    fprintf(stderr, "attempted on line %d in file \"%s\"\n", line, filename);
    exit(EXIT_FAILURE);
  }
  return p;
}

GtRWLock* gt_rwlock_new(void)
{
  GtRWLock *rwlock;
  GT_UNUSED int rval;
  /* XXX: can we use gt_malloc() here? */
  rwlock = thread_xmalloc(sizeof (pthread_rwlock_t), __FILE__, __LINE__);
  /* initialize read/write lock with default attributes */
  rval = pthread_rwlock_init((pthread_rwlock_t*) rwlock, NULL);
  gt_assert(!rval);
  return rwlock;
}

void gt_rwlock_delete(GtRWLock *rwlock)
{
  GT_UNUSED int rval;
  if (!rwlock) return;
  rval = pthread_rwlock_destroy((pthread_rwlock_t*) rwlock);
  gt_assert(!rval);
  free(rwlock);
}

void gt_rwlock_rdlock_func(GtRWLock *rwlock)
{
  GT_UNUSED int rval;
  gt_assert(rwlock);
  rval = pthread_rwlock_rdlock((pthread_rwlock_t*) rwlock);
  gt_assert(!rval);
}

void gt_rwlock_wrlock_func(GtRWLock *rwlock)
{
  GT_UNUSED int rval;
  gt_assert(rwlock);
  rval = pthread_rwlock_wrlock((pthread_rwlock_t*) rwlock);
  gt_assert(!rval);
}

void gt_rwlock_unlock_func(GtRWLock *rwlock)
{
  GT_UNUSED int rval;
  gt_assert(rwlock);
  rval = pthread_rwlock_unlock((pthread_rwlock_t*) rwlock);
  gt_assert(!rval);
}

GtMutex* gt_mutex_new(void)
{
  GtMutex *mutex;
  GT_UNUSED int rval;
  /* XXX: can we use gt_malloc() here? */
  mutex = thread_xmalloc(sizeof (pthread_mutex_t), __FILE__, __LINE__);
  /* initialize mutex with default attributes */
  rval = pthread_mutex_init((pthread_mutex_t*) mutex, NULL);
  gt_assert(!rval);
  return mutex;
}

void gt_mutex_delete(GtMutex *mutex)
{
  GT_UNUSED int rval;
  if (!mutex) return;
  rval = pthread_mutex_destroy((pthread_mutex_t*) mutex);
  gt_assert(!rval);
  free(mutex);
}

void gt_mutex_lock_func(GtMutex *mutex)
{
  GT_UNUSED int rval;
  gt_assert(mutex);
  rval = pthread_mutex_lock((pthread_mutex_t*) mutex);
  gt_assert(!rval);
}

void gt_mutex_unlock_func(GtMutex *mutex)
{
  GT_UNUSED int rval;
  gt_assert(mutex);
  rval = pthread_mutex_unlock((pthread_mutex_t*) mutex);
  gt_assert(!rval);
}

#else

int gt_multithread(GtThreadFunc function, void *data, GT_UNUSED GtError *err)
{
  unsigned int i;
  gt_error_check(err);
  gt_assert(function);
  for (i = 0; i < gt_jobs; i++)
    function(data);
  return 0;
}

GtThread* gt_thread_new(GtThreadFunc function, void *data,
                        GT_UNUSED GtError *err)
{
  GtThread *thread;
  gt_error_check(err);
  gt_assert(function);
  thread = gt_malloc(sizeof (void*));
  function(data);
  return thread;
}

GtRWLock* gt_rwlock_new(void)
{
  return NULL;
}

void gt_rwlock_delete(GT_UNUSED GtRWLock *rwlock)
{
  return;
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
  if (!thread) return;
  gt_free(thread);
}
