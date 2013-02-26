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

#ifndef CLASS_ALLOC_LOCK_H
#define CLASS_ALLOC_LOCK_H

/* Initializes the big lock for runtime class instatiation. */
void    gt_class_alloc_lock_init(void);
/* Cleans static resources required for the runtime class instatiation lock. */
void    gt_class_alloc_lock_clean(void);

/* Marks the beginning of the critical section when instantiating class
   objects. */
#ifdef GT_THREADS_ENABLED
#define gt_class_alloc_lock_enter() \
        gt_class_alloc_lock_enter_func()
void    gt_class_alloc_lock_enter_func(void);
#else
#define gt_class_alloc_lock_enter() \
        ((void) 0)
#endif

/* Marks the end of the critical section when instantiating class objects. */
#ifdef GT_THREADS_ENABLED
#define gt_class_alloc_lock_leave() \
        gt_class_alloc_lock_leave_func()
void    gt_class_alloc_lock_leave_func(void);
#else
#define gt_class_alloc_lock_leave() \
        ((void) 0)
#endif

#endif
