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

#ifndef QUEUE_H
#define QUEUE_H

#include <stdio.h>
#include "core/error.h"

typedef struct GT_Queue GT_Queue;

typedef int (*GT_QueueProcessor)(void **elem, void *info, GT_Error*);

GT_Queue*     gt_queue_new(void);
void          gt_queue_delete(GT_Queue*);
void          gt_queue_delete_with_contents(GT_Queue*);
void          gt_queue_add(GT_Queue*, void*);
void*         gt_queue_get(GT_Queue*);
void*         gt_queue_head(GT_Queue*);
/* Remove <elem> from <queue> (<elem> has to be in <queue>).
   Thereby <queue> is traversed in reverse order, leading to O(gt_queue_size)
   worst-case running time. */
void          gt_queue_remove(GT_Queue *queue, void *elem);
/* Iterate over all elements in <queue> and call <gt_queue_processor> with them.
   <info> and <err> are passed to <queue_processor>.
   If <queue_processor> returns a value != 0, the iteration is stopped and the
   return value of <queue_processor> is returned. */
int           gt_queue_iterate(GT_Queue *queue,
                               GT_QueueProcessor queue_processor,
                               void *info, GT_Error *err);
/* Similar to <gt_queue_iterate>, except that the <queue> is traversed in
   reverse order. */
int           gt_queue_iterate_reverse(GT_Queue *queue,
                                       GT_QueueProcessor queue_processor,
                                       void *info, GT_Error *err);
unsigned long gt_queue_size(const GT_Queue*);
int           gt_queue_unit_test(GT_Error*);

#endif
