/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include "core/thread_api.h"

static GtMutex* gt_class_alloc_lock = NULL;

void gt_class_alloc_lock_init(void)
{
  gt_class_alloc_lock = gt_mutex_new();
}

void gt_class_alloc_lock_clean(void)
{
  gt_mutex_delete(gt_class_alloc_lock);
}

void gt_class_alloc_lock_enter_func(void)
{
  gt_assert(gt_class_alloc_lock);
  gt_mutex_lock(gt_class_alloc_lock);
}

void gt_class_alloc_lock_leave_func(void)
{
  gt_assert(gt_class_alloc_lock);
  gt_mutex_unlock(gt_class_alloc_lock);
}
