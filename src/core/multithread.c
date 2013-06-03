/*
  Copyright (c) 2010, 2013 Gordon Gremme <gordon@gremme.org>

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

#include "core/array_api.h"
#include "core/multithread_api.h"
#include "core/unused_api.h"

#ifdef GT_THREADS_ENABLED

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

#endif
