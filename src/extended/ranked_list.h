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

#ifndef RANKED_LIST_H
#define RANKED_LIST_H

#include "extended/rbtree.h"

/* The <GtRankedListIter> class implements an iterator over the elements of
   a <GtRankedList>. */
typedef GtRBTreeIter GtRankedListIter;

/* The <GtRankedList> class holds a dynamic sorted collection of <n> items
   such that the <n> items with the highest rank are kept. */
typedef struct GtRankedList GtRankedList;

/* Returns a new <GtRankedList> object with maximumsize <maxsize>. The
   comparator function <comparefunction> is used to define an order on the
   inserted elements, with <compareinfo> to be used as additional external
   data. */
GtRankedList* gt_ranked_list_new(unsigned long maxsize,
                                 GtCompareWithData comparefunction,
                                 void *compareinfo);

/* Inserts <elem> into <ranked_list>. Does not take ownership of <elem>. */
void          gt_ranked_list_insert(GtRankedList *ranked_list, void *elem);

/* Returns the element in <ranked_list> with the highest rank. */
void*         gt_ranked_list_first(const GtRankedList *ranked_list);

/* Returns the element in <ranked_list> with the lowest rank. */
void*         gt_ranked_list_last(const GtRankedList *ranked_list);

/* Returns the number of elements currently stored in <ranked_list>. */
unsigned long gt_ranked_list_size(const GtRankedList *ranked_list);

/* Deletes <ranked_list> and frees all associated memory. */
void          gt_ranked_list_delete(GtRankedList *ranked_list);

int           gt_ranked_list_unit_test(GtError *err);

/* Returns a new <GtRankedListIter>, initialized to the element in <ranked_list>
   with the highest rank. The object is deleted along with <ranked_list> and
   thus does not need to be freed. */
GtRankedListIter* gt_ranked_list_iter_new_from_first(GtRankedList *ranked_list);

/* Returns a new <GtRankedListIter>, initialized to the element in <ranked_list>
   with the lowest rank. The object is deleted along with <ranked_list> and
   thus does not need to be freed */
GtRankedListIter* gt_ranked_list_iter_new_from_last(GtRankedList *ranked_list);

/* Returns the next element according to <ranked_list_iter>. */
void*             gt_ranked_list_iter_next(GtRankedListIter *ranked_list_iter);

/* Returns the previous element according to <ranked_list_iter>. */
void*             gt_ranked_list_iter_prev(GtRankedListIter *ranked_list_iter);

#endif
