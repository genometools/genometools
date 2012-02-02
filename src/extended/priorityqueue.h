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

#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

typedef struct GtPriorityQueue GtPriorityQueue;
typedef long GtPQelementtype;

GtPriorityQueue *priorityqueue_new(unsigned long maxnumofelements);
void priorityqueue_add(GtPriorityQueue *pq, GtPQelementtype value);
GtPQelementtype priorityqueue_delete_min(GtPriorityQueue *pq);
GtPQelementtype priorityqueue_find_min(const GtPriorityQueue *pq);
bool priorityqueue_is_empty(const GtPriorityQueue *pq);
bool priorityqueue_is_full(const GtPriorityQueue *pq);
void priorityqueue_delete(GtPriorityQueue *pq);

#endif
