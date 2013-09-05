/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef PRIORITY_QUEUE_H
#define PRIORITY_QUEUE_H

#include "core/fptr_api.h"

typedef struct GtPriorityQueue GtPriorityQueue;

GtPriorityQueue* gt_priority_queue_new(GtCompare cmpfun,
                                       GtUword maxnumofelements);
void             gt_priority_queue_add(GtPriorityQueue *pq,
                                       void *value);
void*            gt_priority_queue_extract_min(GtPriorityQueue *pq);
const void*      gt_priority_queue_find_min(const GtPriorityQueue *pq);
bool             gt_priority_queue_is_empty(const GtPriorityQueue *pq);
bool             gt_priority_queue_is_full(const GtPriorityQueue *pq);
void             gt_priority_queue_delete(GtPriorityQueue *pq);
int              gt_priority_queue_unit_test(GtError *err);

#endif
